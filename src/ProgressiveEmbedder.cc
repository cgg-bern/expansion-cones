#include "ProgressiveEmbedder.hh"
#include <random>

namespace OpenVolumeMesh{



ProgressiveEmbedder::ProgressiveEmbedder(TetrahedralMesh& domain_mesh,
                                         TetrahedralMesh& codomain_mesh,
                                         const std::string& mesh_name,
                                         const std::string& output_file_path):
        domain_mesh_(domain_mesh),
        codomain_mesh_(codomain_mesh),
        mesh_name_(mesh_name),
        output_file_path_(output_file_path){ }


int ProgressiveEmbedder::shrinkAndExpand(TetrahedralMesh& domain_mesh,
                                         TetrahedralMesh& codomain_mesh,
                                         const std::string& mesh_name,
                                         const std::string& output_json_path,
                                         int debug_expander){

    ProgressiveEmbedder embedder(domain_mesh, codomain_mesh, mesh_name, output_json_path);
    return embedder.shrink_and_expand(debug_expander);

}

bool ProgressiveEmbedder::is_codomain_boundary_valid() const {
    bool boundary_ok(true);
    for(auto f: codomain_mesh_.faces()){
        if(codomain_mesh_.is_boundary(f)){
            auto f_vertices = codomain_mesh_.get_halfface_vertices(codomain_mesh_.halfface_handle(f, 0));
            auto tri = CGAL_Triangle3(OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[0])),
                                      OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[1])),
                                      OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[2])));

            if(CGAL::collinear(OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[0])),
                               OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[1])),
                               OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[2])))){
                std::cout<<" --> input boundary face "<<f<<": "<<f_vertices<<" is collinear "<<std::endl;
                //return -1;
                boundary_ok = false;
            }

            auto op_v = codomain_mesh_.halfface_opposite_vertex(codomain_mesh_.halfface_handle(f, 0));
            if(!op_v.is_valid()){
                op_v = codomain_mesh_.halfface_opposite_vertex(codomain_mesh_.halfface_handle(f, 1));
            }
            if(CGAL::coplanar(OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[0])),
                              OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[1])),
                              OVMvec3ToCGALPoint3(codomain_mesh_.vertex(f_vertices[2])),
                              OVMvec3ToCGALPoint3(codomain_mesh_.vertex(op_v)))){
                std::cout<<" --> input boundary face "<<f<<": "<<f_vertices<<" is coplanar with cluster "<<std::endl;
                boundary_ok = false;
                //return -1;

            }

            /*if(tri.is_degenerate()){
                std::cout<<" --> input boundary face "<<f<<": "<<f_vertices<<" is degenerate "<<std::endl;
                return 12123;
            }*/
        }
    }
    return boundary_ok;
}


int ProgressiveEmbedder::shrink_and_expand(int debug_expander){

    if(PEHelpers::BCI_edges_count(codomain_mesh_)){
        std::cout<<" ERROR - cannot Shrink-and-Expand mesh with boundary-connecting interior edges"<<std::endl;
        return -1;
    }

    Vec3d centroid = {0,0,0};

    if(!Expander::at_least_one_expanded_vertex(codomain_mesh_)){


        std::cout<<" WARNING - USING ORIGIN AS CENTROID"<<std::endl;

        for(auto v: codomain_mesh_.vertices()){
            if(!codomain_mesh_.is_boundary(v)){
                codomain_mesh_.set_vertex(v, centroid);
            }
        }
        //SphereMapper::center_mesh(codomain_mesh_);
        //SphereMapper::map_interior_to_centroid(codomain_mesh_);


        if(!is_codomain_boundary_valid()){
            std::cout<<" WARNING - found some collinear boundary faces or some coplanar boundary tets"<<std::endl;
            std::cout<<" WARNING - trying Chebyshev center as alternative for cluster initialization"<<std::endl;
            std::cout<<" WARNING - Performances for this initialization have NOT been tested"<<std::endl;

            ExpansionCone cone = codomain_mesh_;
            VertexPosition cheb_centroid;
            auto exp = cone.is_geo_expandable(cheb_centroid);
            if(exp){
                std::cout<<" --> Codomain boundary is NOT star-shaped"<<std::endl;
                return -1;
            }

            centroid = vec2vec(cheb_centroid);
            for(auto v: codomain_mesh_.vertices()){
                if(!codomain_mesh_.is_boundary(v)){
                    codomain_mesh_.set_vertex(v, centroid);
                }
            }

            if(!is_codomain_boundary_valid()) {
                std::cout<<" ERROR - codomain boundary is still not valid, even with Chebyshev center."<<std::endl;
                std::cout<<" This can typically happen if your kernel is too small and the truncation to double failed"<<std::endl;
                std::cout<<" My suggestion is to use boundary conditions for which the origin is inside the kernel."<<std::endl;
                std::cout<<" If you have further questions, please reach out to yours truly (Valentin Nigolian)"<<std::endl;
                return -1;
            }
        }

        if(!SphereMapper::is_boundary_star_shaped_and_non_degenerate(codomain_mesh_)){
            std::cout<<" ERROR - Input mesh contains degenerate boundary triangles or is not star-shaped after boundary mapping"<<std::endl;
            return -1;
        }

    }else{
        std::cout<<" ==> MESH WAS ALREADY PARTIALLY EXPANDED, SKIPPING BOUNDARY CHECK"<<std::endl;
   }


    if(BadTetFinder::meshContainsFlippedTets(codomain_mesh_)){

        std::cout<<" ERROR - mesh contains flipped tets with centroid = "<<centroid<<std::endl;

        auto init_flipped_tets_count = BadTetFinder::findBadTets(codomain_mesh_).second.size();
        std::cout<<" - initial flipped tets: "<<init_flipped_tets_count<<std::endl;

        //IO::FileManager fm;
        //fm.writeFile("a_bad_centroid.ovm", codomain_mesh_);
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
    Expander expander(codomain_mesh_, domain_mesh_, mesh_name_, !debug_expander);
    SAE_result = expander.fully_expand_mesh(output_file_path_,
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
    return codomain_mesh_.n_vertices() - count_boundary_vertices();
}

int ProgressiveEmbedder::count_boundary_vertices() const{

    int count(0);
    for(auto v: codomain_mesh_.vertices()){
        count += codomain_mesh_.is_boundary(v);
    }
    return count;
}


int ProgressiveEmbedder::count_interior_edges() const{
    int count(0);
    for(auto e: codomain_mesh_.edges()){
        count += is_interior(e);
    }
    return count;
}


int ProgressiveEmbedder::count_boundary_edges() const{
    int count(0);
    for(auto e: codomain_mesh_.edges()){
        count += codomain_mesh_.is_boundary(e);
    }
    return count;
}

bool ProgressiveEmbedder::is_interior(const EdgeHandle& eh) const{
    return !codomain_mesh_.is_boundary(codomain_mesh_.edge(eh).from_vertex()) &&
           !codomain_mesh_.is_boundary(codomain_mesh_.edge(eh).to_vertex());
}



void ProgressiveEmbedder::setup_json_exporter(JsonExporter& exporter){

    //11 vertex valence histogram
    const int valence_cap(20);
    std::vector<int> vertex_valence_histogram(valence_cap+1);
    for(auto v: codomain_mesh_.vertices()){
        auto val = codomain_mesh_.valence(v);
        if(val > valence_cap){
            vertex_valence_histogram[valence_cap]++;
        }else{
            vertex_valence_histogram[val]++;
        }
    }

    //12 edge valence histogram
    std::vector<int> edge_valence_histogram(valence_cap+1);
    for(auto e: codomain_mesh_.edges()){
        auto val = codomain_mesh_.valence(e);
        if(val > valence_cap){
            edge_valence_histogram[valence_cap]++;
        }else{
            edge_valence_histogram[val]++;
        }
    }

    double kappa, stdev;
    //TopoHelper::connectivity_matrix_condition_number(codomain_mesh_, kappa);

    double avg_val_grad;
    auto val_grad_prop = codomain_mesh_.request_vertex_property<double>();
    TopoHelper::compute_valence_gradient(codomain_mesh_, val_grad_prop, avg_val_grad);

    exporter.write("mesh_name", mesh_name_);
    exporter.write_vector("mesh", std::vector<size_t>({codomain_mesh_.n_vertices(), codomain_mesh_.n_edges(), codomain_mesh_.n_faces(), codomain_mesh_.n_cells()}));
    exporter.write("n_interior_vertices", count_interior_vertices());
    exporter.write("n_interior_edges", count_interior_edges());
    exporter.write("n_boundary_edges", count_boundary_edges());
    exporter.write_vector("vertex_valence_histogram", vertex_valence_histogram);
    exporter.write_vector("edge_valence_histogram", edge_valence_histogram);
    exporter.write("n_initial_flipped_tets", BadTetFinder::findBadTets(codomain_mesh_).second.size());
    exporter.write("initial_avg_val_grad", avg_val_grad);
}



void ProgressiveEmbedder::export_shrinkage_data_to_json(JsonExporter& exporter,
                                                        const float shrinkage_time_s,
                                                        const int shrinkage_result){

    std::vector<int> cluster_sizes;

    int shrunk_edges_count(0);
    for(auto e: codomain_mesh_.edges()){
        if(codomain_mesh_.vertex(codomain_mesh_.edge(e).from_vertex()) == codomain_mesh_.vertex(codomain_mesh_.edge(e).to_vertex())){
            shrunk_edges_count++;
        }
    }

    exporter.write("shrinkage_time_s", shrinkage_time_s);
    exporter.write("shrinkage_result", shrinkage_result);
    exporter.write("full_shrinkage_necessary", -1);

    if(cluster_sizes.empty()){
        exporter.write_vector("cluster_sizes", std::vector<int>({count_interior_vertices()}));
    }else{
        exporter.write_vector("cluster_sizes", cluster_sizes);
    }
    exporter.write("n_shrunk_edges", shrunk_edges_count);
    exporter.write("n_degenerate_tets_after_shrinkage", BadTetFinder::findBadTets(codomain_mesh_).first.size());
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

    auto bad_tets = BadTetFinder::findBadTets(codomain_mesh_);


    double kappa, stdev;
    //TopoHelper::connectivity_matrix_condition_number(codomain_mesh_, kappa);

    double avg_val_grad;
    auto val_grad_prop = codomain_mesh_.request_vertex_property<double>();
    TopoHelper::compute_valence_gradient(codomain_mesh_, val_grad_prop, avg_val_grad);

    exporter.write("expansion_time_s", expansion_time_s);
    exporter.write("expansion_result", expansion_result);
    exporter.write("n_vertices_after_expansion", codomain_mesh_.n_logical_vertices());
    exporter.write("n_edges_after_expansion", codomain_mesh_.n_logical_edges());
    exporter.write("n_faces_after_expansion", codomain_mesh_.n_logical_faces());
    exporter.write("n_cells_after_expansion", codomain_mesh_.n_logical_cells());
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


