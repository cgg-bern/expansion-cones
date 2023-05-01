#include "ProgEmbeddingHelpers.hh"

namespace OpenVolumeMesh{


void PEHelpers::compute_and_output_quality_statistics(TetrahedralMesh& mesh, const std::string& output_file_path){

    auto str_domain_exact_pos_prop = mesh.request_vertex_property<std::string>("str_domain_exact_position");
    auto str_exact_pos_prop = mesh.request_vertex_property<std::string>("str_exact_position");

    auto   domain_exact_pos_prop = mesh.request_vertex_property<VertexPosition>();
    auto codomain_exact_pos_prop = mesh.request_vertex_property<VertexPosition>();

    for(auto v: mesh.vertices()){

        auto   domain_exact_pos = string_to_position(str_domain_exact_pos_prop[v]);
        auto codomain_exact_pos = string_to_position(str_exact_pos_prop[v]);

        auto diff_norm = (vec2vec(codomain_exact_pos) - mesh.vertex(v)).norm();

        if(diff_norm > 0){
            std::cout<<" ERROR - exact pos for vertex "<<v<<" is "<<codomain_exact_pos<<" = "<<vec2vec(codomain_exact_pos)<<
                       " and mesh pos is "<<mesh.vertex(v)<<
                       ", diff = "<<(vec2vec(codomain_exact_pos)-mesh.vertex(v)).norm()<<std::endl;
            return;
        }

          domain_exact_pos_prop[v] =   domain_exact_pos;
        codomain_exact_pos_prop[v] = codomain_exact_pos;
    }

    std::ofstream output_file(output_file_path/*, std::ios_base::app*/);
    if(!output_file.is_open()){
        std::cout<<" ERROR - couldn't open stat file "<<output_file_path<<std::endl;
        return;
    }


    ExactType total_domain_volume(0), total_codomain_volume(0);
    for(auto c: mesh.cells()){
        total_domain_volume   += OVMtetToCGALtet(mesh,   domain_exact_pos_prop, c).volume();
        total_codomain_volume += OVMtetToCGALtet(mesh, codomain_exact_pos_prop, c).volume();
    }

    for(auto c: mesh.cells()){
        auto cell_vertices = mesh.get_cell_vertices(c);

        std::vector<VertexPosition> domain_positions;
        std::vector<VertexPosition> codomain_positions;
        /*std::cout<<" ---------------------------------------------------------------"<<std::endl;
        std::cout<<" - computing quality for cell "<<cell_vertices<<std::endl;
        std::cout<<" -   domain volume = "<<OVMtetToCGALtet(mesh, domain_exact_pos_prop, c).volume().to_double()<<std::endl;
        std::cout<<" - codomain volume = "<<OVMtetToCGALtet(mesh, codomain_exact_pos_prop, c).volume().to_double()<<std::endl;*/
        //std::cout<<" positions: "<<std::endl;
        for(auto v: cell_vertices){
            domain_positions.push_back(domain_exact_pos_prop[v]);
            codomain_positions.push_back(codomain_exact_pos_prop[v]);
            /*std::cout<<"   - "<<vec2vec(domain_exact_pos_prop[v])<<" -> "<<vec2vec(codomain_exact_pos_prop[v])<<std::endl;
            std::cout<<"   - "<<domain_exact_pos_prop[v]<<" -> "<<codomain_exact_pos_prop[v]<<std::endl;
            std::cout<<" ---------------------------------------------------------------"<<std::endl;*/
        }

        if(domain_positions.size() != 4){
            std::cout<<" ERROR - only "<<cell_vertices.size()<<" vertices in cell "<<c<<": "<<cell_vertices<<std::endl;
            return;
        }


        ExactType conformal, areal;
        conformal_and_areal_components(domain_positions,
                                       codomain_positions,
                                       total_domain_volume,
                                       total_codomain_volume,
                                       conformal,
                                       areal);

        output_file<<conformal.to_double()<<";"<<areal.to_double()<<std::endl;
        //output_file<<conformal<<";"<<areal<<std::endl;
        //std::cout<<" --> conformal = "<<conformal.to_double()<<", areal = "<<areal.to_double()<<std::endl;

    }


    output_file.close();
}


void PEHelpers::compute_and_output_quality_statistics(const TetrahedralMesh& domain_mesh,
                                                      const TetrahedralMesh& codomain_mesh,
                                                      const std::string& output_file_path){


    std::ofstream output_file(output_file_path/*, std::ios_base::app*/);
    if(!output_file.is_open()){
        std::cout<<" ERROR - couldn't open stat file "<<output_file_path<<std::endl;
        return;
    }
    double total_domain_volume(0), total_codomain_volume(0);
    for(auto c: domain_mesh.cells()){
        total_domain_volume   += OVMtetToCGALtet(  domain_mesh, c).volume();
        total_codomain_volume += OVMtetToCGALtet(codomain_mesh, c).volume();
    }

    std::cout<<" total   domain volume = "<<total_domain_volume<<std::endl;
    std::cout<<" total codomain volume = "<<total_codomain_volume<<std::endl;

    for(auto c: domain_mesh.cells()){
        auto cell_vertices = domain_mesh.get_cell_vertices(c);

        std::vector<Vec3d> domain_positions;
        std::vector<Vec3d> codomain_positions;

        for(auto v: cell_vertices){
            domain_positions.push_back(domain_mesh.vertex(v));
            codomain_positions.push_back(codomain_mesh.vertex(v));
        }

        double conformal(0), areal(0);
        conformal_and_areal_components(domain_positions,
                                       codomain_positions,
                                       total_domain_volume,
                                       total_codomain_volume,
                                       conformal,
                                       areal);

        /*auto d_vol = OVMtetToCGALtet(  domain_mesh, c).volume();
        auto cd_vol =OVMtetToCGALtet(codomain_mesh, c).volume();*/

        //output_file<<conformal.to_double()<<";"<<areal.to_double()<<std::endl;
        output_file<<conformal<<";"<<areal<<";"<<std::endl;


    }


    output_file.close();
}



VertexPosition PEHelpers::string_to_position(const std::string& str_exact_pos){

    //PRINT_IF_NOT_SILENT("----------------------- "<<std::endl);
    //PRINT_IF_NOT_SILENT(" - exact pos string: "<<str_exact_pos<<std::endl);

    auto v_pos_x_pos = str_exact_pos.find(';');
    auto v_pos_x     = str_exact_pos.substr(0, v_pos_x_pos);

    auto v_pos_y_pos = str_exact_pos.find(';', v_pos_x_pos + 1);
    auto v_pos_y     = str_exact_pos.substr(v_pos_x_pos + 1, (v_pos_y_pos - v_pos_x_pos - 1));

    auto v_pos_z     = str_exact_pos.substr(v_pos_y_pos + 1);

    //PRINT_IF_NOT_SILENT(" - exact pos x string: "<<v_pos_x<<std::endl);
    //PRINT_IF_NOT_SILENT(" - exact pos y string: "<<v_pos_y<<std::endl);
    //PRINT_IF_NOT_SILENT(" - exact pos z string: "<<v_pos_z<<std::endl);

    CGAL::Gmpq x,y,z;

    std::stringstream ss_x;
    ss_x<<v_pos_x;
    ss_x>>x;

    std::stringstream ss_y;
    ss_y<<v_pos_y;
    ss_y>>y;

    std::stringstream ss_z;
    ss_z<<v_pos_z;
    ss_z>>z;

    /*PRINT_IF_NOT_SILENT(" - exact pos x Gmpq: "<<x<<std::endl);
    PRINT_IF_NOT_SILENT(" - exact pos y Gmpq: "<<y<<std::endl);
    PRINT_IF_NOT_SILENT(" - exact pos z Gmpq: "<<z<<std::endl);*/

    return {x,y,z};
}




    int PEHelpers::DBCI_vertices_count(const TetrahedralMesh& mesh){

        int count(0);

        for(auto v: mesh.vertices()){
            if(!mesh.is_boundary(v)){
                int boundary_neighbors_count(0);
                for(auto out_he: mesh.outgoing_halfedges(v)){
                    if(boundary_neighbors_count >= 2){
                        break;
                    }
                    auto neighbor = mesh.to_vertex_handle(out_he);

                    if(mesh.is_boundary(neighbor)){
                        boundary_neighbors_count++;
                    }
                }

                if(boundary_neighbors_count >= 2){
                    count++;
                }
            }
        }

        return count;
    }

    int PEHelpers::BCI_edges_count(const TetrahedralMesh& mesh){

        int count(0);
        for(auto e: mesh.edges()){
            if(!mesh.is_boundary(e)){
                count += mesh.is_boundary(mesh.edge(e).from_vertex()) &&
                         mesh.is_boundary((mesh.edge(e).to_vertex()));
            }
        }
        return count;
    }




    int PEHelpers::degenerate_boundary_faces_count(const TetrahedralMesh& mesh){

        int count(0);

        for(auto fh: mesh.faces()){
            if(mesh.is_boundary(fh)){
                auto f_vertices = mesh.get_halfface_vertices(mesh.halfface_handle(fh, 0));
                CGAL_Triangle3 tri(OVMvec3ToCGALPoint3(mesh.vertex(f_vertices[0])),
                                   OVMvec3ToCGALPoint3(mesh.vertex(f_vertices[1])),
                                   OVMvec3ToCGALPoint3(mesh.vertex(f_vertices[2])));

                count += tri.is_degenerate();

                /*if(tri.is_degenerate()){
                    std::cout<<" -> face "<<fh<<": "<<f_vertices<<" is degenerate: "<<std::endl;
                    for(auto v: f_vertices){
                        std::cout<<"    -- "<<v<<" at "<<std::setprecision(30)<<mesh.vertex(v)<<" is boundary: "<<mesh.is_boundary(v)<<std::endl;
                    }
                }*/
            }
        }


        return count;
    }



    size_t PEHelpers::find_maximum_shift_factor(const VertexPosition& new_vertex_position){

        int byte_size = max_min_size(new_vertex_position);

        int min_shift_factor = 7 * byte_size;
        int max_shift_factor = 8 * byte_size;

        int shift_factor = max_shift_factor;

        //std::cout<<" initial min-max shift factor = "<<min_shift_factor<<"-"<<max_shift_factor<<std::endl;

        bool found_best_factor(false);

        while(!found_best_factor){


            auto shifted_x_denom = (new_vertex_position[0].denominator() >> shift_factor);
            auto shifted_y_denom = (new_vertex_position[1].denominator() >> shift_factor);
            auto shifted_z_denom = (new_vertex_position[2].denominator() >> shift_factor);

            /*std::cout<<" -----"<<std::endl;
            std::cout<<" - min-max shift = "<<min_shift_factor<<" - "<<max_shift_factor<<std::endl;
            std::cout<<" - shift factor = "<<shift_factor<<std::endl;
            std::cout<<" - shifted x,y,z = "<<
                       shifted_x_denom<<
                       ", "<<shifted_y_denom<<
                       ", "<<shifted_z_denom<<std::endl;*/


            //if the shift is too high, we use it as upper bound
            //NOTE: we use 10 as a minimum denominator so new positions
            //are much less likely to create slivers
            if(shifted_x_denom <= 10 ||
               shifted_y_denom <= 10 ||
               shifted_z_denom <= 10){

                max_shift_factor = shift_factor;
            }else{
                min_shift_factor = shift_factor;
            }

            shift_factor = (min_shift_factor + max_shift_factor) / 2;
            //std::cout<<" - updated shift factor = "<<shift_factor<<std::endl;

            if(shift_factor <= min_shift_factor){
                found_best_factor = true;
            }
        }

        //std::cout<<" final max shift factor = "<<shift_factor<<std::endl;

        return shift_factor;
    }



    void PEHelpers::lower_precision(unsigned long shift_factor,
                                    VertexPosition& pos){

        //std::cout<<" --- lowering precision with shift factor = "<<shift_factor<<"..."<<std::endl;
        for(auto i(0); i<3; i++){

            auto size = std::min(pos[i].numerator().bit_size(), pos[i].denominator().bit_size());

            //std::cout<<"        size "<<i<<" = "<<size<<std::endl;
            if(shift_factor < size){
                pos[i] = {pos[i].numerator()  >> shift_factor,
                          pos[i].denominator() >> shift_factor};
                //std::cout<<"         shifted"<<std::endl;
            }
            //size = pos[i].numerator().bit_size();
            //std::cout<<"        updated size "<<i<<" = "<<size<<std::endl;


        }
        //std::cout<<" --- done"<<std::endl;
    }



    int PEHelpers::count_entangled_edges(const TetrahedralMesh& mesh){

        int count(0);

        for(auto e: mesh.edges()){
            if(!mesh.is_boundary(e)){
                continue;
            }
            //std::cout<<" --------------------------------------"<<std::endl;
            //std::cout<<" - checking edge "<<mesh.edge(e)<<std::endl;
            auto from_v = mesh.edge(e).from_vertex();
            auto to_v = mesh.edge(e).to_vertex();

            double dihedral_angle_sum(0);
            for(auto ec_it = mesh.ec_iter(e); ec_it.valid(); ec_it++){
                //std::cout<<" -------- checking cell "<<mesh.get_cell_vertices(*ec_it)<<std::endl;
                std::vector<HalfFaceHandle> dihedral_faces;
                std::vector<VertexHandle> op_vertices;
                for(auto chf_it = mesh.chf_iter(*ec_it); chf_it.valid(); chf_it++){
                    for(auto hfe_it = mesh.hfe_iter(*chf_it); hfe_it.valid(); hfe_it++){
                        if(*hfe_it == e){
                            dihedral_faces.push_back(*chf_it);
                            for(auto v: mesh.get_halfface_vertices(*chf_it)){
                                if(v != from_v && v != to_v){
                                    op_vertices.push_back(v);
                                }
                            }
                            break;
                        }
                    }
                }
                if(dihedral_faces.size() != 2){
                    std::cout<<" error - couldn't find two dihedral halffaces"<<std::endl;
                    return -1;
                }

                auto from_v_pos = mesh.vertex(from_v);
                auto   to_v_pos = mesh.vertex(to_v);

                auto next_v_pos1 = mesh.vertex(op_vertices[0]);
                auto next_v_pos2 = mesh.vertex(op_vertices[1]);

                /*std::cout<<" he: "<<mesh.halfedge(he)<<std::endl;
                std::cout<<" next he1: "<<mesh.halfedge(mesh.next_halfedge_in_halfface(he, dihedral_faces[0]))<<std::endl;
                std::cout<<" next he2: "<<mesh.halfedge(mesh.next_halfedge_in_halfface(mesh.halfedge_handle(e,1), dihedral_faces[1]))<<std::endl;
*/
                /*std::cout<<" - c = "<<next_v_pos1<<std::endl;
                std::cout<<" - d = "<<next_v_pos2<<std::endl;

                std::cout<<" - b-a = "<<(to_v_pos - from_v_pos)<<std::endl;
                std::cout<<" - c-a = "<<(next_v_pos1 - from_v_pos)<<std::endl;
                std::cout<<" - d-a = "<<(next_v_pos2 - from_v_pos)<<std::endl;*/

                auto normal1 = (to_v_pos - from_v_pos).cross(next_v_pos1 - from_v_pos).normalize();
                auto normal2 = (to_v_pos - from_v_pos).cross(next_v_pos2 - from_v_pos).normalize();

                //std::cout<<" - normal for halfface "<<mesh.get_halfface_vertices(dihedral_faces[0])<<" = "<<normal1<<std::endl;
                //std::cout<<" - normal for halfface "<<mesh.get_halfface_vertices(dihedral_faces[1])<<" = "<<normal2<<std::endl;

                auto theta = acos(normal1.dot(normal2));

                //std::cout<<"  --> dihedral angle = "<<theta<<"rad = "<<(theta * 180 / M_PI)<<"°"<<std::endl;

                dihedral_angle_sum += theta;
            }
            if(dihedral_angle_sum > 2*M_PI){
                std::cout<<" --> dihedral angles sum = "<<dihedral_angle_sum<<", diff to 2pi = "<<std::abs(dihedral_angle_sum - 2*M_PI)<<std::endl;
                count++;
            }
        }

        return count;
    }




}//namespace OVM



