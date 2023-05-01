#include "ProgressiveEmbedder.hh"
#include <random>

namespace OpenVolumeMesh{



ProgressiveEmbedder::ProgressiveEmbedder(TetrahedralMesh& mesh,
                                         const std::string& mesh_name,
                                         const std::string& output_file_path) :
    mesh_(mesh),
    mesh_name_(mesh_name),
    output_file_path_(output_file_path)
    /*cluster_prop_(mesh_.request_vertex_property<int>()),
    to_shrink_prop_(mesh_.request_edge_property<bool>()),
    shrunk_prop_(mesh_.request_edge_property<bool>()),
    zero_volume_cells_prop_(mesh_.request_cell_property<bool>()),
    changed_neighboring_cluster_at_last_iteration_prop_(mesh_.request_vertex_property<bool>()),
    total_negative_volume_in_1_ring_prop_(mesh_.request_edge_property<double>()),
    computed_1_ring_volume_at_least_once_prop_(mesh_.request_edge_property<bool>()),
    best_position_prop_(mesh_.request_edge_property<VertexPosition>()),
    vertex_position_prop_(mesh_.request_vertex_property<VertexPosition>()),
    n_initial_interior_edges_(count_interior_edges())*/{

}


int ProgressiveEmbedder::shrinkAndExpand(TetrahedralMesh& mesh,
                                         const std::string& mesh_name,
                                         const std::string& output_file_path,
                                         int boundary_mapping_method,
                                         int shrinkage_method){

    ProgressiveEmbedder embedder(mesh, mesh_name, output_file_path);
    return embedder.shrink_and_expand(boundary_mapping_method, shrinkage_method);

}



int ProgressiveEmbedder::shrink_and_expand(int boundary_mapping_method,
                                           int shrinkage_method){

    if(PEHelpers::BCI_edges_count(mesh_)){
        std::cout<<" ERROR - cannot Shrink-and-Expand mesh with boundary-connecting interior edges"<<std::endl;
        return -1;
    }

    TetrahedralMesh input_mesh = mesh_;

   if(!Expander::at_least_one_expanded_vertex(mesh_) && boundary_mapping_method >= 0){

        SphereMapper::center_mesh(mesh_);
        SphereMapper::map_interior_to_centroid(mesh_);

        for(auto f: mesh_.faces()){
            if(mesh_.is_boundary(f)){
                auto f_vertices = mesh_.get_halfface_vertices(mesh_.halfface_handle(f, 0));
                auto tri = CGAL_Triangle3(OVMvec3ToCGALPoint3(mesh_.vertex(f_vertices[0])),
                        OVMvec3ToCGALPoint3(mesh_.vertex(f_vertices[1])),
                        OVMvec3ToCGALPoint3(mesh_.vertex(f_vertices[2])));

                if(CGAL::collinear(OVMvec3ToCGALPoint3(mesh_.vertex(f_vertices[0])),
                                   OVMvec3ToCGALPoint3(mesh_.vertex(f_vertices[1])),
                                   OVMvec3ToCGALPoint3(mesh_.vertex(f_vertices[2])))){
                    std::cout<<" --> input boundary face "<<f<<": "<<f_vertices<<" is collinear "<<std::endl;
                    return -1;

                }

                /*if(tri.is_degenerate()){
                    std::cout<<" --> input boundary face "<<f<<": "<<f_vertices<<" is degenerate "<<std::endl;
                    return 12123;
                }*/
            }
        }

        switch(boundary_mapping_method){

        case TET_ONLY:{
            //TetMapper::map_boundary_to_tet(mesh_);
            //SphereMapper::map_interior_to_centroid(mesh_);
            map_to_unit_tet(mesh_);
            break;
        }
        case STIFF_TET:{
            //TetMapper::map_boundary_to_tet(mesh_);
            //SphereMapper::map_interior_to_centroid(mesh_);
            map_to_stiff_unit_tet(mesh_);
            break;
        }
        case TET_TO_SPHERE:{
            /*TetMapper::map_boundary_to_tet(mesh_);
            SphereMapper::map_interior_to_centroid(mesh_);
            SphereMapper::project_boundary_to_sphere(mesh_,1.);*/
            map_to_unit_ball_using_tet(mesh_);
            break;
        }
        case TET_TO_RANDOM_STAR_SHAPE:{

            map_to_random_star_shape_using_tet(mesh_);
            break;
        }
        default:{
            std::cout<<" ERROR - unhandled boundary mapping case "<<boundary_mapping_method<<std::endl;
            return -1;
        }
        }

        if(!SphereMapper::is_boundary_star_shaped_and_non_degenerate(mesh_)){
            std::cout<<" ERROR - Input mesh contains degenerate boundary triangles or is not star-shaped after boundary mapping"<<std::endl;
            return -1;
        }

    }else{
        std::cout<<" ==> MESH WAS ALREADY PARTIALLY EXPANDED OR BOUNDARY-MAPPING EXPLICITLY SKIPPED, SKIPPING BOUNDARY-MAPPING"<<std::endl;
   }


    auto boundary_area = SphereMapper::boundary_area(mesh_);
    std::cout<<" WARNING - set centroid to origin"<<std::endl;
    Vec3d centroid = {0,0,0};
    //auto result = mesh_cone.find_max_min_volume_center(centroid, 1, boundary_area);

    //centroid = {0,0,0};

    for(auto v: mesh_.vertices()){
        if(!mesh_.is_boundary(v)){
            mesh_.set_vertex(v, centroid);
        }
    }

    std::cout<<" - mesh centroid = "<<vec2vec(centroid)<<std::endl;

    if(BadTetFinder::meshContainsFlippedTets(mesh_)){

        std::cout<<" ERROR - mesh contains flipped tets with centroid = "<<centroid<<std::endl;

        auto init_flipped_tets_count = BadTetFinder::findBadTets(mesh_).second.size();
        std::cout<<" - initial flipped tets with min-max volume centroid: "<<init_flipped_tets_count<<std::endl;

        IO::FileManager fm;
        fm.writeFile("a_bad_centroid.ovm", mesh_);
        return -1;
    }


    int SAE_result(-1);

    //prepare the data exporter
    std::stringstream output_stream;
    JsonExporter exporter(output_stream);
    setup_json_exporter(exporter);

    export_shrinkage_data_to_json(exporter, 0, 0);

    auto expansion_start_time = std::chrono::high_resolution_clock::now();
    double expansion_time_s(0);
    Expander expander(mesh_, input_mesh, mesh_name_, false);
    SAE_result = expander.fully_expand_mesh("expansion_details_"+output_file_path_,
                                            DEFAULT_EXPANSION_SCHEME,
                                            expansion_time_s);

    auto expansion_end_time = std::chrono::high_resolution_clock::now();
    float expansion_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(expansion_end_time - expansion_start_time).count() /1000000;

    std::cout<<" TEMP - initial centroid = "<<vec2vec(centroid)<<std::endl;

    std::cout<<" done with local expansions in "<<expansion_duration_s<<" seconds"<<std::endl;


    export_expansion_data_to_json(exporter,
                                  expansion_time_s,
                                  SAE_result,
                                  expander.left_to_expand_count(),
                                  expander.iteration_count(),
                                  expander.star_shapification_count(),
                                  expander.cluster_expansion_count(),
                                  expander.cluster_star_shapification_count());

    exporter.close();

    std::ofstream output_data_file;
    output_data_file.open(output_file_path_);
    if(!output_data_file.is_open()){
        std::cerr<<" ERROR - couldn't open file "<<output_file_path_<<std::endl;
    }/*else{
            std::cout<<" JSON file: "<<output_file_path_<<std::endl;
        }*/
    //exporter.close();
    //output_data_file << std::endl <<"}"<<std::endl;

    output_data_file << exporter.stream().rdbuf();

    output_data_file.close();


    if(SAE_result == -1){
        return -2;
    }else{
        return SAE_result;
    }
}




void ProgressiveEmbedder::interior_uniform_smoothing(TetrahedralMesh& mesh_,
                                                     const double delta_eps,
                                                     const int max_iterations){

    std::cout<<" ================================================================ "<<std::endl;
    std::cout<<" UNIFORMLY SMOOTHING INTERIOR WITH DELTA EPSILON = "<<delta_eps<<std::endl;

    auto next_pos_prop = mesh_.request_vertex_property<TetrahedralMesh::PointT>();

    int i(0);
    double max_delta(delta_eps + 1);

    while(i < max_iterations && max_delta > delta_eps){

        max_delta = 0;

        for(auto v: mesh_.vertices()){
            if(!mesh_.is_boundary(v)){
                next_pos_prop[v] = {0,0,0};
                double n(0);
                for(auto vv_it = mesh_.vv_iter(v); vv_it.valid(); vv_it++){
                    next_pos_prop[v] += mesh_.vertex(*vv_it);
                    n++;
                }
                next_pos_prop[v] /= n;
                auto delta = (next_pos_prop[v] - mesh_.vertex(v)).norm();
                if(delta > max_delta){
                    max_delta = delta;
                }
            }
        }

        for(auto v: mesh_.vertices()){
            if(!mesh_.is_boundary(v)){
                mesh_.set_vertex(v, next_pos_prop[v]);
            }
        }

        i++;
    }


    std::cout<<" ... done smoothing after "<<i<<" iterations "<<std::endl;
    std::cout<<" ================================================================ "<<std::endl;

}


int ProgressiveEmbedder::map_to_unit_tet(TetrahedralMesh& mesh){

    SphereMapper::move_centroid_and_scale_mesh_boundary_surface(mesh, {0,0,0}, 1);
    SphereMapper::map_interior_to_centroid(mesh);
    TetMapper::map_boundary_to_tet(mesh);
    SphereMapper::map_interior_to_centroid(mesh);

    return mesh_contains_flipped_tets_or_degenerate_boundary_triangles(mesh);
}


int ProgressiveEmbedder::map_to_stiff_unit_tet(TetrahedralMesh& mesh){

    std::cout<<" mapping to stiff tet"<<std::endl;
    SphereMapper::move_centroid_and_scale_mesh_boundary_surface(mesh, {0,0,0}, 1);
    SphereMapper::map_interior_to_centroid(mesh);
    TetMapper::map_boundary_to_tet(mesh, true);
    SphereMapper::map_interior_to_centroid(mesh);

    return mesh_contains_flipped_tets_or_degenerate_boundary_triangles(mesh);
}


int ProgressiveEmbedder::map_to_unit_ball_using_tet(TetrahedralMesh& mesh){
    auto tet_mapping_result = map_to_unit_tet(mesh);
    if(tet_mapping_result){
        std::cout<<" couldn't map to unit tet, not mapping to unit ball"<<std::endl;
        return 1;
    }

    std::cout<<" - projecting boundary to sphere..."<<std::endl;
    SphereMapper::project_boundary_to_sphere(mesh, 1);

    return mesh_contains_flipped_tets_or_degenerate_boundary_triangles(mesh);
}


int ProgressiveEmbedder::map_to_random_star_shape_using_tet(TetrahedralMesh& mesh){
    auto sphere_mapping_result = map_to_unit_ball_using_tet(mesh);
    if(sphere_mapping_result){
        std::cout<<" couldn't map to unit sphere, not mapping to random star shape"<<std::endl;
        return 1;
    }


    const double max_variation_amplitude(10);
    const int max_variation(10000);
    const int seed(1234);
    std::mt19937 rng(seed);

    std::cout<<" - applying random scaling to boundary vertices. Max scaling = "<<max_variation_amplitude<<std::endl;

    for(auto v: mesh.vertices()){
        if(mesh.is_boundary(v)){
            double amplitude = 1+(rng() % max_variation) * (max_variation_amplitude-1) / (double)max_variation;
            mesh.set_vertex(v, mesh.vertex(v) * amplitude);
        }
    }

    std::cout<<" centroid = "<<SphereMapper::boundary_centroid(mesh)<<std::endl;

    return mesh_contains_flipped_tets_or_degenerate_boundary_triangles(mesh);
}


int ProgressiveEmbedder::mesh_contains_flipped_tets_or_degenerate_boundary_triangles(TetrahedralMesh& mesh){

    auto flipped_tets = BadTetFinder::findBadTets(mesh).second;
    //std::cout<<" flipped tets: "<<flipped_tets.size()<<std::endl;
    if(!flipped_tets.empty()){
        std::cout<<" --> mesh contains flipped tets"<<std::endl;
        return 1;
    }

    int deg_boundary_face_count = PEHelpers::degenerate_boundary_faces_count(mesh);
    //std::cout<<" deg boundary tris : "<<deg_boundary_face_count<<std::endl;
    if(deg_boundary_face_count){
        std::cout<<" --> degenerate triangles on tet boundary"<<std::endl;
        return 1;
    }

    return 0;
}




int ProgressiveEmbedder::count_interior_vertices() const{
    return mesh_.n_vertices() - count_boundary_vertices();
}

int ProgressiveEmbedder::count_boundary_vertices() const{

    int count(0);
    for(auto v: mesh_.vertices()){
        count += mesh_.is_boundary(v);
    }
    return count;
}


int ProgressiveEmbedder::count_interior_edges() const{
    int count(0);
    for(auto e: mesh_.edges()){
        count += is_interior(e);
    }
    return count;
}


int ProgressiveEmbedder::count_boundary_edges() const{
    int count(0);
    for(auto e: mesh_.edges()){
        count += mesh_.is_boundary(e);
    }
    return count;
}

bool ProgressiveEmbedder::is_interior(const EdgeHandle& eh) const{
    return !mesh_.is_boundary(mesh_.edge(eh).from_vertex()) &&
           !mesh_.is_boundary(mesh_.edge(eh).to_vertex());
}



void ProgressiveEmbedder::setup_json_exporter(JsonExporter& exporter){

    //11 vertex valence histogram
    const int valence_cap(20);
    std::vector<int> vertex_valence_histogram(valence_cap+1);
    for(auto v: mesh_.vertices()){
        auto val = mesh_.valence(v);
        if(val > valence_cap){
            vertex_valence_histogram[valence_cap]++;
        }else{
            vertex_valence_histogram[val]++;
        }
    }

    //12 edge valence histogram
    std::vector<int> edge_valence_histogram(valence_cap+1);
    for(auto e: mesh_.edges()){
        auto val = mesh_.valence(e);
        if(val > valence_cap){
            edge_valence_histogram[valence_cap]++;
        }else{
            edge_valence_histogram[val]++;
        }
    }

    double kappa, stdev;
    //TopoHelper::connectivity_matrix_condition_number(mesh_, kappa);

    double avg_val_grad;
    auto val_grad_prop = mesh_.request_vertex_property<double>();
    TopoHelper::compute_valence_gradient(mesh_, val_grad_prop, avg_val_grad);

    exporter.write("mesh_name", mesh_name_);
    exporter.write_vector("mesh", std::vector<size_t>({mesh_.n_vertices(), mesh_.n_edges(), mesh_.n_faces(), mesh_.n_cells()}));
    exporter.write("n_interior_vertices", count_interior_vertices());
    exporter.write("n_interior_edges", count_interior_edges());
    exporter.write("n_boundary_edges", count_boundary_edges());
    exporter.write_vector("vertex_valence_histogram", vertex_valence_histogram);
    exporter.write_vector("edge_valence_histogram", edge_valence_histogram);
    exporter.write("n_initial_flipped_tets", BadTetFinder::findBadTets(mesh_).second.size());
    exporter.write("initial_avg_val_grad", avg_val_grad);
}



void ProgressiveEmbedder::export_shrinkage_data_to_json(JsonExporter& exporter,
                                                        const float shrinkage_time_s,
                                                        const int shrinkage_result){

    std::vector<int> cluster_sizes;

    int shrunk_edges_count(0);
    for(auto e: mesh_.edges()){
        if(mesh_.vertex(mesh_.edge(e).from_vertex()) == mesh_.vertex(mesh_.edge(e).to_vertex())){
            shrunk_edges_count++;
        }
    }

    exporter.write("shrinkage_time_s", shrinkage_time_s);
    exporter.write("shrinkage_result", shrinkage_result);
    exporter.write("full_shrinkage_necessary", -1);

    std::cout<<" cluster sizes: ";
    for(auto c: cluster_sizes){
        std::cout<<" "<<c;
    }
    std::cout<<std::endl;

    if(cluster_sizes.empty()){
        exporter.write_vector("cluster_sizes", std::vector<int>({count_interior_vertices()}));
    }else{
        exporter.write_vector("cluster_sizes", cluster_sizes);
    }
    exporter.write("n_shrunk_edges", shrunk_edges_count);
    exporter.write("n_degenerate_tets_after_shrinkage", BadTetFinder::findBadTets(mesh_).first.size());
    exporter.write("n_LPs_solved", -1);
    exporter.write("n_back_up_positions", -1);
    exporter.write("n_skipped_flipped_tet_edges", -1);

}



void ProgressiveEmbedder::export_expansion_data_to_json(JsonExporter& exporter,
                                                        const float expansion_time_s,
                                                        const int expansion_result,
                                                        const int left_to_expand_count,
                                                        const int expansion_iteration_count,
                                                        const int star_shapifications_count,
                                                        const int cluster_expansions_count,
                                                        const int cluster_star_shapifications_count){

    auto bad_tets = BadTetFinder::findBadTets(mesh_);


    double kappa, stdev;
    //TopoHelper::connectivity_matrix_condition_number(mesh_, kappa);

    double avg_val_grad;
    auto val_grad_prop = mesh_.request_vertex_property<double>();
    TopoHelper::compute_valence_gradient(mesh_, val_grad_prop, avg_val_grad);

    exporter.write("expansion_time_s", expansion_time_s);
    exporter.write("expansion_result", expansion_result);
    exporter.write("n_vertices_after_expansion", mesh_.n_logical_vertices());
    exporter.write("n_edges_after_expansion", mesh_.n_logical_edges());
    exporter.write("n_faces_after_expansion", mesh_.n_logical_faces());
    exporter.write("n_cells_after_expansion", mesh_.n_logical_cells());
    exporter.write("n_remaining_unexpanded_vertices", left_to_expand_count);
    exporter.write("n_degenerate_tets_after_expansion", bad_tets.first.size());
    exporter.write("n_flipped_tets_after_expansion", bad_tets.second.size());
    exporter.write("n_expansion_iterations", expansion_iteration_count);
    exporter.write("n_star_shapifications", star_shapifications_count);
    exporter.write("n_cluster_expansions", cluster_expansions_count);
    exporter.write("n_cluster_star_shapifications", cluster_star_shapifications_count);
    exporter.write("final_avg_val_grad", avg_val_grad);
    std::cout<<" EXPORTED EXPANSION DATA"<<std::endl;
}
}


