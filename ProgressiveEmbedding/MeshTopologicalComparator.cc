
#include "MeshTopologicalComparator.hh"


namespace OpenVolumeMesh{


bool mark_sub_tets_as_part_of_tet(const TetrahedralMesh& reference_mesh,
                                  TetrahedralMesh& subdivided_mesh,
                                  const SubTriangleMap& sub_triangle_map,
                                  const CellHandle& reference_cell_handle,
                                  CellPropertyT<CellHandle>& reference_tet_prop);

bool store_reference_neighborhood(const TetrahedralMesh& reference_mesh,
                                 const CellHandle& ch,
                                 CellPropertyT<std::set<CellHandle>>& reference_neighbor_tets_prop);

bool check_neighborhood(const TetrahedralMesh& subdivided_mesh,
                        const CellHandle sub_tet_handle,
                        const CellPropertyT<CellHandle>& reference_tet_prop,
                        const CellPropertyT<std::set<CellHandle>>& reference_neighbor_tets_prop);



bool is_subdivided_topologically_equivalent(const TetrahedralMesh& reference_mesh,
                                            const TetrahedralMesh& subdivided_mesh,
                                            const SubTriangleMap& sub_triangle_map){

    std::cout<<" ============================================================"<<std::endl;
    std::cout<<" CHECKING IF SUBDIVIDED MESH IS EQUIVALENT TO REFERENCE MESH"<<std::endl;

    TopoHelper reference_mesh_helper(reference_mesh);
    TopoHelper subdivided_mesh_helper(subdivided_mesh);


    auto subdivided_mesh_copy = subdivided_mesh;
    auto reference_tet_prop = subdivided_mesh_copy.request_cell_property<CellHandle>();
    auto reference_neighbor_tets_prop = subdivided_mesh_copy.request_cell_property<std::set<CellHandle>>();


    //mark all subtets as parts of their reference tet
    for(auto ch: reference_mesh.cells()){

        const auto& cell_vertices = reference_mesh.get_cell_vertices(ch);
        if(!mark_sub_tets_as_part_of_tet(reference_mesh,
                                         subdivided_mesh_copy,
                                         sub_triangle_map,
                                         ch,
                                         reference_tet_prop)){
            std::cout<<" --> reference cell "<<ch<<": ("<<cell_vertices<<")"<<
                       " has no ball-homeomorphic region equivalent"<<std::endl;
            return false;

        }

        store_reference_neighborhood(reference_mesh,
                                    ch,
                                    reference_neighbor_tets_prop);

    }

    int bad_subtet_count(0);
    for(auto ch: subdivided_mesh.cells()){
        if(reference_tet_prop[ch].idx() == -1){
            bad_subtet_count++;

            std::cout<<" - sub-tet "<<ch<<":"
                       " ("<<subdivided_mesh.get_cell_vertices(ch)<<")"
                       " is is no reference tet... checking its hfs..."<<std::endl;

            for(auto chf_it = subdivided_mesh.chf_iter(ch); chf_it.valid(); chf_it++){
                auto hf_vertices = subdivided_mesh.get_halfface_vertices(*chf_it);

                std::cout<<" -- looking for hf "<<*chf_it<<"("<<hf_vertices<<") in sub-triangle map..."<<std::endl;

                for(auto key_value: sub_triangle_map){

                    for(auto subface: key_value.second){
                        if(SubTriangleFace(hf_vertices) == subface){
                            std::cout<<" ---> found subface "<<subface<<" in face "<<key_value.first<<std::endl;
                            auto hf_ref_handle = reference_mesh.halfface(key_value.first.get_vertices());

                            std::cout<<" ---> hf handle for those vertices: "<<hf_ref_handle<<std::endl;
                            if(hf_ref_handle.idx() != -1){
                                auto cell_ref_handle = reference_mesh.incident_cell(hf_ref_handle);
                                std::cout<<" ----> cell handle for this hf: "<<cell_ref_handle<<":"
                                           " ("<<reference_mesh.get_cell_vertices(cell_ref_handle)<<")"<<std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout<<" sub-tets without a reference tet: "<<bad_subtet_count<<std::endl;


    //check that all sub-tets are indeed part of an reference tet
    //and that their neighbors are either
    //  (i)  part of the same reference tet
    //  (ii) part of a reference neighbor tet
    for(auto ch: subdivided_mesh.cells()){

        auto cell_vertices = subdivided_mesh.get_cell_vertices(ch);

        if(!check_neighborhood(subdivided_mesh,
                               ch,
                               reference_tet_prop,
                               reference_neighbor_tets_prop)){
            std::cout<<" --> cell "<<ch<<": ("<<cell_vertices<<")"<<
                       " of subdivided mesh has a neighborhood inconsistent with the reference mesh"<<std::endl;
            return false;
        }
    }




    std::cout<<"... DONE! IT IS INDEED EQUIVALENT"<<std::endl;
    std::cout<<" ============================================================"<<std::endl;

    return true;
}




bool mark_sub_tets_as_part_of_tet(const TetrahedralMesh& reference_mesh,
                                  TetrahedralMesh& subdivided_mesh,
                                  const SubTriangleMap& sub_triangle_map,
                                  const CellHandle& reference_cell_handle,
                                  CellPropertyT<CellHandle>& reference_tet_prop){

    /*std::cout<<" ---------"<<std::endl;
    std::cout<<" marking subtets of cell "<<reference_cell_handle<<":"<<
               " ("<<TopoHelper(reference_mesh).cell_vertices(reference_cell_handle)<<")..."<<std::endl;*/


    auto boundary_subface_prop = subdivided_mesh.request_halfface_property<bool>();

    auto reference_cell_in_subdivided_mesh = TopoHelper::cell_exists(subdivided_mesh,
                                                                     reference_mesh.get_cell_vertices(reference_cell_handle));

    std::queue<CellHandle> subtets_to_mark;

    //if the cell is present as-is in the subdivided mesh, we're done
    if(reference_cell_in_subdivided_mesh.idx() != -1){
        //std::cout<<" - reference cell is present as-is in the subdivided mesh"<<std::endl;
        reference_tet_prop[reference_cell_in_subdivided_mesh] = reference_cell_handle;
        return true;
    }else{
        //otherwise we need to mark all sub-tets

        for(auto chf_it = reference_mesh.chf_iter(reference_cell_handle); chf_it.valid(); chf_it++){
            //std::cout<<" - checking hf "<<*chf_it<<": ("<<reference_mesh.get_halfface_vertices(*chf_it)<<")"<<std::endl;

            auto reference_hf_vertices = reference_mesh.get_halfface_vertices(*chf_it);
            if(reference_hf_vertices.size() != 3){
                std::cerr<<" ERROR - got "<<reference_hf_vertices.size()<<" vertices for hf "<<*chf_it<<std::endl;
                return false;
            }
            SubTriangleFace reference_hf(reference_hf_vertices);
            //std::cout<<" - reference hf: "<<reference_hf<<std::endl;

            auto reference_hf_handle_in_subdivided_mesh = subdivided_mesh.halfface(reference_hf_vertices);

            //if the hf is present in the subdivided mesh, we mark its incident cell
            if(reference_hf_handle_in_subdivided_mesh.idx() != -1){

                /*std::cout<<" --> which exists in the subdivided mesh,"
                           " adding its incident cell "<<subdivided_mesh.incident_cell(reference_hf_handle_in_subdivided_mesh)<<":"<<
                           " ("<<subdivided_mesh.get_cell_vertices(subdivided_mesh.incident_cell(reference_hf_handle_in_subdivided_mesh))<<")"<<
                           " for marking"<<std::endl;*/

                boundary_subface_prop[reference_hf_handle_in_subdivided_mesh] = true;

                subtets_to_mark.push(subdivided_mesh.incident_cell(reference_hf_handle_in_subdivided_mesh));

            }else if(sub_triangle_map.contains(reference_hf)){
                //std::cout<<" -- which is present in the sub-triangle map, adding its subfaces for marking"<<std::endl;

                //else we need to mark all tets incident to its subfaces
                SubFaceSet subfaces;
                sub_triangle_map.extract_all_subfaces(reference_hf,
                                                      subfaces);

                for(auto subface: subfaces){
                    auto subface_handle = subdivided_mesh.halfface(subface.get_vertices());
                    if(subface_handle.idx() != -1){
                        boundary_subface_prop[subface_handle] = true;
                        subtets_to_mark.push(subdivided_mesh.incident_cell(subface_handle));
                    }else{
                        std::cout<<" --> subface "<<subface<<" is not present in the subdivided mesh"<<std::endl;
                        return false;
                    }
                }
            }else{
                std::cout<<" --> halfface "<<*chf_it<<":"<<
                           " ("<<reference_mesh.get_halfface_vertices(*chf_it)<<")"<<
                           "doesn't exist in the subdivided mesh or in the sub-triangle map"<<std::endl;
                return false;
            }
        }
    }

    //tet-mesh consisting of all the subtets contained in the reference cell
    TetrahedralMesh subtets_mesh;

    auto visited = subdivided_mesh.request_cell_property<bool>();

    //std::cout<<" subtets to mark: "<<subtets_to_mark.size()<<std::endl;

    while(!subtets_to_mark.empty()){

        auto current_subtet_to_mark = subtets_to_mark.front();
        subtets_to_mark.pop();

        //std::cout<<" - current cell : "<<current_subtet_to_mark<<std::endl;
        if(current_subtet_to_mark.idx() == -1 ||
                current_subtet_to_mark.idx() >= (int)subdivided_mesh.n_cells()){
            std::cerr<<" ERROR - sub-tets list contains invalid cell handle"<<std::endl;
            return false;
        }
        //std::cout<<" - vertices: "<<subdivided_mesh.get_cell_vertices(current_subtet_to_mark)<<std::endl;

        if(!visited[current_subtet_to_mark]){
            //mark the current subtet
            reference_tet_prop[current_subtet_to_mark] = reference_cell_handle;
            visited[current_subtet_to_mark] = true;

            //add the neighbors to the queue
            for(auto chf_it = subdivided_mesh.chf_iter(current_subtet_to_mark); chf_it.valid(); chf_it++){
                if(!boundary_subface_prop[*chf_it]){
                    auto opposite_hf = subdivided_mesh.opposite_halfface_handle(*chf_it);
                    if(opposite_hf.idx() != -1){
                        auto neighbor_cell = subdivided_mesh.incident_cell(opposite_hf);
                        if(neighbor_cell.idx() != -1 &&
                                !visited[neighbor_cell]){
                            subtets_to_mark.push(neighbor_cell);
                            //std::cout<<" added neighbor cell "<<neighbor_cell<<": "<<subdivided_mesh.get_cell_vertices(neighbor_cell)<<std::endl;
                        }
                    }
                }
            }
        }
    }


    //and check that all sub-tets make up a topological ball


    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ---------"<<std::endl;
    return true;
}



bool store_reference_neighborhood(const TetrahedralMesh& reference_mesh,
                                 const CellHandle&      ch,
                                 CellPropertyT<std::set<CellHandle>>& reference_neighbor_tets_prop){

    /*std::cout<<" ---------"<<std::endl;
    std::cout<<" storing neighbors to tet "<<ch<<
               " ("<<TopoHelper(reference_mesh).cell_vertices(ch)<<")"
               " in property..."<<std::endl;*/


    for(auto cc_it = reference_mesh.cc_iter(ch); cc_it.valid(); cc_it++){
        reference_neighbor_tets_prop[ch].insert(*cc_it);
    }

    //std::cout<<" neighbors to cell "<<ch<<": "<<reference_neighbor_tets_prop[ch].size()<<std::endl;

    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ---------"<<std::endl;
    return false;
}



bool check_neighborhood(const TetrahedralMesh& subdivided_mesh,
                        const CellHandle sub_tet_handle,
                        const CellPropertyT<CellHandle>& reference_tet_prop,
                        const CellPropertyT<std::set<CellHandle>>& reference_neighbor_tets_prop){

    /*std::cout<<" ---------"<<std::endl;
    std::cout<<" checking that neighborhood of sub-tet "<<sub_tet_handle<<": "<<
               "("<<TopoHelper(subdivided_mesh).cell_vertices(sub_tet_handle)<<")"
               " is consistent with reference mesh..."<<std::endl;

    std::cout<<" --> reference tet: "<<reference_tet_prop[sub_tet_handle]<<std::endl;
    std::cout<<" --> reference tet neighbors : "<<reference_neighbor_tets_prop[sub_tet_handle].size()<<std::endl;
    for(auto c: reference_neighbor_tets_prop[sub_tet_handle]){
        std::cout<<c<<" ";
    }
    std::cout<<std::endl;*/

    if(reference_tet_prop[sub_tet_handle].idx() != -1){

        auto sub_tet_reference_tet = reference_tet_prop[sub_tet_handle];

        auto& reference_neighborhood = reference_neighbor_tets_prop[sub_tet_reference_tet];

        //then for all neighbors of the subtet,
        //check that they're either in the same reference tet
        //or in a neighbor of the reference tet
        for(auto cc_it = subdivided_mesh.cc_iter(sub_tet_handle); cc_it.valid(); cc_it++){

            const auto& reference_tet_of_neighbor = reference_tet_prop[*cc_it];
            //std::cout<<" -- neighbor: "<<*cc_it<<", in ref tet: "<<reference_tet_of_neighbor<<std::endl;

            if(reference_tet_of_neighbor != sub_tet_reference_tet &&
                    reference_neighborhood.find(reference_tet_of_neighbor) == reference_neighborhood.end()){
                std::cout<<" --> neighbor cell "<<*cc_it<<
                           " is part of the reference cell "<<reference_tet_of_neighbor<<
                           " which is not a neighbor of the reference tet of the sub-tet"<<std::endl;
                return false;

            }
        }
    }else{
        std::cout<<" sub-tet "<<sub_tet_handle<<":"
                   " ("<<subdivided_mesh.get_cell_vertices(sub_tet_handle)<<")"
                   " is part of no reference tet"<<std::endl;
        return false;
    }



    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ---------"<<std::endl;
    return true;
}

}



