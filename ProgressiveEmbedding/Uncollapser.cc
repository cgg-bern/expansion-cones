#include "Uncollapser.hh"


namespace OpenVolumeMesh{

Uncollapser::Uncollapser(TetrahedralMesh& mesh,
                         const SubTriangleMap& sub_triangle_map) :
    mesh_(mesh),
    sub_triangle_map_(sub_triangle_map),
    topo_helper_(mesh){

}




const SubTriangleMap& Uncollapser::get_sub_triangle_map() const{
    return sub_triangle_map_;
}


bool Uncollapser::uncollapse_sequence(CollapseSequence collapse_sequence){

    std::cout<<" -------------------------------------- "<<std::endl;
    std::cout<<"  UNDOING COLLAPSE SEQUENCE "<<std::endl;

    /*int i(0);
    for(auto collapse: collapse_sequence){
        std::cout<<"  COLLAPSE N° "<<i<<
                   " : ("<<collapse.from_vertex<<"-"<<collapse.to_vertex<<
                   ") -> "<<collapse.replacement_vertex<<
                   (collapse.TEMP_after_uncollapse ? " -> after uncollapse":"")<<std::endl;
        i++;
    }*/


    //return false;

    auto start_time = std::chrono::high_resolution_clock::now();


    int uncollapse_count(0);
    int total_uncollapse_count(collapse_sequence.size());
    VertexHandle uncollapse_result(0);
    while(collapse_sequence.size() &&
          uncollapse_result.idx() != -1){

        //std::cout<<" .............................................................. PERFORMING uncollapse "<<uncollapse_count<<"/"<<total_uncollapse_count<<std::endl;


        uncollapse_result = uncollapse(collapse_sequence.back());


        collapse_sequence.pop_back();
        uncollapse_count++;

        if(uncollapse_count && !(uncollapse_count % 100)){
            std::cout<<"... done "<<uncollapse_count<<" uncollapses"<<std::endl;
        }
    }

    if(uncollapse_result.idx() == -1){
        std::cout<<" ERROR WHILE UNDOING COLLAPSE SEQUENCE "<<std::endl;


        return false;
    }else{

        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout<<" DONE! FULLY UNCOLLAPSED MESH AFTER "<<total_uncollapse_count<<
                   " IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" -------------------------------------- "<<std::endl;
        return true;
    }

}



bool Uncollapser::uncollapse_sequence(SACSequence split_and_collapse_sequence){


    std::cout<<" -------------------------------------- "<<std::endl;
    std::cout<<"  UNDOING SPLIT-AND-COLLAPSE SEQUENCE "<<std::endl;

    /*int i(0);
    for(auto collapse: collapse_sequence){
        std::cout<<"  COLLAPSE N° "<<i<<
                   " : ("<<collapse.from_vertex<<"-"<<collapse.to_vertex<<
                   ") -> "<<collapse.replacement_vertex<<
                   (collapse.TEMP_after_uncollapse ? " -> after uncollapse":"")<<std::endl;
        i++;
    }*/


    //return false;

    auto start_time = std::chrono::high_resolution_clock::now();


    int split_success_count(0);
    int uncollapse_count(0);
    int total_uncollapse_count(split_and_collapse_sequence.size());
    int total_split_count(0);
    VertexHandle uncollapse_result(0);
    while(split_and_collapse_sequence.size() &&
          uncollapse_result.idx() != -1){

        const auto& split    = split_and_collapse_sequence.back().split;
        const auto& collapse = split_and_collapse_sequence.back().collapse;

        //std::cout<<" .............................................................. PERFORMING uncollapse "<<uncollapse_count<<"/"<<total_uncollapse_count<<std::endl;
        uncollapse_result = uncollapse(collapse);


        //undo split

        if(is_valid(split)){
            auto unsplit_result = unsplit(split);

            if(unsplit_result == UNSPLIT_ERROR){
                std::cerr<<" ERROR WHILE UNDOING SPLIT "<<split<<std::endl;
                return false;
            }else if(unsplit_result == UNSPLIT_SUCCESS){
                split_success_count ++;
            }
            total_split_count++;
        }


        uncollapse_count++;
        split_and_collapse_sequence.pop_back();


        if(uncollapse_count && !(uncollapse_count % 100)){
            std::cout<<"... done "<<uncollapse_count<<" uncollapses"<<std::endl;
        }
    }

    std::cout<<" split count: "<<split_success_count<<"/"<<total_split_count<<std::endl;


    if(uncollapse_result.idx() == -1){
        std::cout<<" ERROR WHILE UNDOING SPLIT-AND-COLLAPSE SEQUENCE "<<std::endl;
        return false;
    }else{

        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout<<" DONE! FULLY UNCOLLAPSED MESH AFTER "<<total_uncollapse_count<<
                   " IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" -------------------------------------- "<<std::endl;
        return true;
    }
}



UNSPLIT_RESULT Uncollapser::unsplit(const EdgeSplit& split){

    //std::cout<<" -------------------------------------------------- "<<std::endl;
    //std::cout<<" =========================> UNDOING SPLIT "<<split<<std::endl;

    if(!is_valid(split)){
        std::cerr<<" ERROR - tried to unsplit an invalid split "<<split<<std::endl;
        return UNSPLIT_ERROR;
    }

    UNSPLIT_RESULT result = UNSPLIT_ERROR;

    const auto& from_vertex   = split.from_vertex;
    const auto& middle_vertex = split.middle_vertex;
    const auto& to_vertex     = split.to_vertex;

    const auto& first_half = topo_helper_.halfedge_exists(middle_vertex, from_vertex);
    const auto& second_half = topo_helper_.halfedge_exists(middle_vertex, to_vertex);

    if(first_half.idx() == -1 ||
            second_half.idx() == -1){

        //std::cout<<" at least one of the two parts of the original edge no longer exists"<<std::endl;
        return UNSPLIT_TOPO_FAILURE;
    }

    if(!is_unsplittable(split, first_half) ||
            !is_unsplittable(split, second_half)){
        //std::cout<<" at least one of the two parts of the original edge is not unsplittable"<<std::endl;
        return UNSPLIT_TOPO_FAILURE;
    }

    HalfEdgeHandle to_unsplit(-1);

    if(GeometryHelper::is_collapsible(mesh_, second_half)){
        to_unsplit = second_half;
        result = UNSPLIT_SUCCESS;

    }else if(GeometryHelper::is_collapsible(mesh_, second_half)){
        to_unsplit = first_half;
        result = UNSPLIT_SUCCESS;

    }else{
        //std::cout<<" at least one of the two parts of the original edge is geometrically not unsplittable"<<std::endl;
        return UNSPLIT_GEO_FAILURE;
    }


    mesh_.collapse_edge(to_unsplit);
    register_unsplit_in_subtriangle_map(split);

    //std::cout<<" =========== DONE UNSPLIT "<<split<<std::endl;
    //std::cout<<" -------------------------------------------------- "<<std::endl;

    return result;
}



bool Uncollapser::is_valid(const EdgeSplit& split) const{
    return split.from_vertex.idx() != -1 &&
            split.middle_vertex.idx() != -1 &&
            split.to_vertex.idx() != -1;
}


bool Uncollapser::is_unsplittable(const EdgeSplit& split,
                                  const HalfEdgeHandle& heh) const{

    //std::cout<<" -----------------"<<std::endl;
    //std::cout<<" checking wether halfedge "<<mesh_.halfedge(heh)<<" is unsplittable"<<std::endl;

    const auto& from_vertex = mesh_.from_vertex_handle(heh);
    const auto& to_vertex   = mesh_.to_vertex_handle(heh);

    for(auto ev: split.equatorial_vertices){
        SubTriangleFace spoke_triangle({from_vertex,
                                        to_vertex,
                                        ev});

        //std::cout<<" - checking spoke triangle "<<spoke_triangle<<std::endl;

        if(sub_triangle_map_.contains(spoke_triangle)){
            //std::cout<<" --> contained in the sub-triangle map and thus split"<<std::endl;
            return false;
        }
    }

    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" -----------------"<<std::endl;

    return true;
}



bool Uncollapser::register_unsplit_in_subtriangle_map(const EdgeSplit& split){

    //std::cout<<" -----------------"<<std::endl;
    //std::cout<<" registering unsplit in sub-triangle map"<<std::endl;

    const auto& from_vertex = split.from_vertex;
    const auto& to_vertex   = split.to_vertex;

    for(auto ev: split.equatorial_vertices){
        SubTriangleFace spoke_triangle({from_vertex,
                                       to_vertex,
                                       ev});

        //std::cout<<" - erasing spoke triangle "<<spoke_triangle<<std::endl;
        sub_triangle_map_.erase(spoke_triangle);

    }


    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" -----------------"<<std::endl;

    return true;
}


#if 0


int Uncollapser::experimental_unsplit(const EdgeSplit& split){
    std::cout<<" -------------------------------------------------- "<<std::endl;
    std::cout<<" =========================> UNDOING SPLIT "<<split<<std::endl;

    //first, gather all halfedges lying on all triangles around the edge to unsplit
    int triangle_index(0);
    std::vector<std::pair<HalfEdgeHandle,int>> halfedges_to_unsplit;
    for(auto ev: split.equatorial_vertices){
        SubTriangleFace spoke_triangle({split.from_vertex,
                                        split.to_vertex,
                                        ev});

        gather_triangle_halfedges(spoke_triangle,
                                  triangle_index,
                                  halfedges_to_unsplit);
        triangle_index++;
    }

    gather_split_center_subhalfedges(split,
                                     halfedges_to_unsplit);




    const int total_unsplits_to_do_count = halfedges_to_unsplit.size();

    int iteration_count(0);
    int unsplit_count(0);
    int total_unsplit_count(0);
    //then try to unsplit them iteratively
    do{

        unsplit_count = 0;
        std::cout<<" ------------------"<<std::endl;
        std::cout<<" - iteration "<<iteration_count<<std::endl;
        std::cout<<" - hes to unsplit: "<<halfedges_to_unsplit.size()<<std::endl;

        std::vector<std::pair<HalfEdgeHandle,int>> failed_unsplits;

        for(auto current_pair: halfedges_to_unsplit){
            const auto& current_he = current_pair.first;

            std::cout<<" -- checking he "<<mesh_.halfedge(current_he)<<std::endl;

            const auto current_e = mesh_.edge_handle(current_he);

            if(!mesh_.is_deleted(current_he)){
                if(topo_helper_.isCollapsible(current_e) &&
                        GeometryHelper::is_collapsible(mesh_, current_he) &&
                        is_unsplittable(current_pair, split)){
                    std::cout<<" ---> fully unsplittable, performing unsplit"<<std::endl;

                    mesh_.collapse_edge(current_he);
                    if(!register_unsplit_in_subtriangle_map(current_pair, split)){
                        std::cerr<<" ERROR - COULDN'T PROPERLY REGISTER UNSPLIT UN MAP."<<std::endl;
                        return -1;
                    }

                    //unsplit_hes.push_back(current_he);
                    unsplit_count++;
                    total_unsplit_count++;

                }else{
                    failed_unsplits.push_back(current_pair);
                    std::cout<<" ---> not unsplittable ";
                    if(!topo_helper_.isCollapsible(current_e)){
                        std::cout<<" (non-link)";
                    }
                    if(!GeometryHelper::is_collapsible(mesh_, current_he)){
                        std::cout<<" (non geometrically-collapsible)";
                    }
                    if(!is_unsplittable(current_pair, split)){
                        std::cout<<" (non topologically-collapsible)";
                    }

                    std::cout<<std::endl;
                }
            }else{
                std::cout<<" --> he is deleted, skipping"<<std::endl;
            }
        }
        std::cout<<" - successful unsplits: "<<unsplit_count<<"/"<<total_unsplits_to_do_count<<std::endl;
        std::cout<<" -     failed unsplits: "<<failed_unsplits.size()<<"/"<<total_unsplits_to_do_count<<std::endl;

        halfedges_to_unsplit = failed_unsplits;
        iteration_count++;

    }while(unsplit_count && !halfedges_to_unsplit.empty());



    std::cout<<" -- DONE UNSPLIT "<<split<<" AFTER "<<iteration_count<<" ITERATIONS"<<std::endl;
    std::cout<<" -------------------------------------------------- "<<std::endl;

    return unsplit_count;
}








bool Uncollapser::gather_triangle_halfedges(const SubTriangleFace& triangle,
                                            const int& triangle_index,
                                            std::vector<std::pair<HalfEdgeHandle,int>>& halfedges_on_triangle){
    std::cout<<" ----------"<<std::endl;
    std::cout<<" gathering triangle "<<triangle<<" internal halfedges"<<std::endl;

    SubFaceSet subfaces;
    sub_triangle_map_.extract_all_subfaces(triangle,
                                           subfaces);

    auto boundary_vertices_prop = mesh_.request_vertex_property<bool>();
    mark_equatorial_disc_boundary_vertices(subfaces,
                                           boundary_vertices_prop);

    auto triangle_vertices = recover_subface_set_vertices(subfaces);
    std::cout<<" - triangle vertices: "<<triangle_vertices<<std::endl;
    std::cout<<" - boundary vertices: ";
    for(auto v: triangle_vertices){
        if(boundary_vertices_prop[v]){
            std::cout<<" "<<v;
        }
    }
    std::cout<<std::endl;


    auto triangle_vertices_prop = mesh_.request_vertex_property<bool>();
    for(auto current_vertex: triangle_vertices){
        triangle_vertices_prop[current_vertex] = true;
    }


    int gathered_halfedges_count(0);
    for(auto current_vertex: triangle_vertices){

        std::cout<<" -- checking vertex "<<current_vertex<<" neighborhood"<<std::endl;
        if(!boundary_vertices_prop[current_vertex]){
            for(auto out_he_it: mesh_.outgoing_halfedges(current_vertex)){
                auto other_vertex = mesh_.to_vertex_handle(out_he_it);
                std::cout<<" --- checking other vertex "<<other_vertex<<std::endl;
                if(triangle_vertices_prop[other_vertex]){
                    std::cout<<" ----> found halfedge "<<mesh_.halfedge(out_he_it)<<std::endl;
                    halfedges_on_triangle.push_back({out_he_it, triangle_index});
                    gathered_halfedges_count++;
                }else{
                    std::cout<<" ----> not on the triangle "<<std::endl;
                }
            }
        }
    }

    std::cout<<"... done, gathered "<<gathered_halfedges_count<<" halfedges"<<std::endl;
    std::cout<<" ----------"<<std::endl;

    return true;
}




bool Uncollapser::gather_split_center_subhalfedges(const EdgeSplit& split,
                                                   std::vector<std::pair<HalfEdgeHandle,int>>& center_subhalfedges){

    std::cout<<" ----------"<<std::endl;
    std::cout<<" gathering split center sub-halfedges..."<<std::endl;


    SubTriangleFace first_spoke({split.from_vertex,
                                 split.to_vertex,
                                 split.equatorial_vertices[0]});


    SubTriangleFace second_spoke({split.from_vertex,
                                 split.to_vertex,
                                 split.equatorial_vertices[1]});

    std::set<VertexHandle> center_vertices;

    sub_triangle_map_.extract_subvertices_intersection(first_spoke,
                                                       second_spoke,
                                                       center_vertices);

    auto on_center_halfedge_prop = mesh_.request_vertex_property<bool>();

    for(auto center_v: center_vertices){
        on_center_halfedge_prop[center_v] = true;
    }
    std::cout<<" - center vertices: "<<center_vertices<<std::endl;


    auto visited = mesh_.request_vertex_property<bool>();


    std::vector<HalfEdgeHandle> half_center_subhalfedges;


    VertexHandle current_vertex = split.from_vertex;
    while(current_vertex.idx() != -1 &&
          half_center_subhalfedges.size() < center_vertices.size()){

        visited[current_vertex] = true;

        VertexHandle next_vertex(-1);

        for(auto out_he_it: mesh_.outgoing_halfedges(current_vertex)){
            auto to_vertex = mesh_.to_vertex_handle(out_he_it);
            if(!visited[to_vertex] &&
                    on_center_halfedge_prop[to_vertex]){

                next_vertex = to_vertex;
                half_center_subhalfedges.push_back(out_he_it);
                std::cout<<" -- found halfedge "<<mesh_.halfedge(out_he_it)<<std::endl;
                break;
            }
        }
        current_vertex = next_vertex;
    }

    if((half_center_subhalfedges.size() + 1) != center_vertices.size()){
        std::cerr<<" ERROR - gathered "<<half_center_subhalfedges.size()<<" halfedges but there are "<<center_vertices.size()<<" center vertices"<<std::endl;
        return false;
    }

    //add the center halfedge and their opposite halfedges
    //but leave out the he going from the from-vertex and
    //the he foind from the to-vertex
    for(auto he: half_center_subhalfedges){
        if(mesh_.from_vertex_handle(he) != split.from_vertex){
            center_subhalfedges.push_back({he,
                                           split.equatorial_vertices.size()});
        }

        if(mesh_.to_vertex_handle(he) != split.to_vertex){
            center_subhalfedges.push_back({mesh_.opposite_halfedge_handle(he),
                                           split.equatorial_vertices.size()});
        }
    }


    std::cout<<"... done, gathered "<<center_subhalfedges.size()<<" halfedges"<<std::endl;
    std::cout<<" ----------"<<std::endl;

    return true;

}



bool Uncollapser::is_unsplittable(const std::pair<HalfEdgeHandle, int>& heh,
                                  const EdgeSplit& split) const{

    std::cout<<" ----------"<<std::endl;
    const int ev_count = split.equatorial_vertices.size();
    const int triangle_index = heh.second;


    std::cout<<" checking wether halfedge "<<heh.first<<" on spoke triangle "<<triangle_index<<" is unsplittable"<<std::endl;

    if(triangle_index < 0){
        std::cerr<<" ERROR - negative triangle index edge"<<std::endl;
        return false;
    }
    if(triangle_index > ev_count){
        std::cerr<<" ERROR - triangle index "<<heh.second<<" with "<<ev_count<<" equatorial vertices"<<std::endl;
        return false;
    }

    const auto he_from_vertex = mesh_.from_vertex_handle(heh.first);
    const auto he_to_vertex = mesh_.to_vertex_handle(heh.first);

    bool unsplittable(false);

    if(triangle_index == ev_count){
        std::cout<<" - center halfedge, checking all sub-maps..."<<std::endl;

        for(auto ev: split.equatorial_vertices){
            SubTriangleFace spoke_triangle({split.from_vertex,
                                           split.to_vertex,
                                           ev});
            std::cout<<" -- checking sub-map of triangle "<<spoke_triangle<<std::endl;
            if(!sub_triangle_map_.extract_sub_map(spoke_triangle).replace_vertex_in_map_and_update(he_from_vertex,
                                                                                                   he_to_vertex)){
                std::cout<<" ---> created bad sub-map, stopping"<<std::endl;
                break;
            }
        }

        std::cout<<" --> all sub-maps are good, he "<<heh.first<<" is unsplittable!"<<std::endl;

        unsplittable = true;
    }else{
        SubTriangleFace spoke_triangle({split.from_vertex,
                                            split.to_vertex,
                                            split.equatorial_vertices[triangle_index]});

        std::cout<<" - spoke halfedge, checking sub-map of triangle "<<spoke_triangle<<std::endl;


        unsplittable = sub_triangle_map_.extract_sub_map(spoke_triangle).replace_vertex_in_map_and_update(he_from_vertex,
                                                                                                          he_to_vertex);

        if(unsplittable){
            std::cout<<" --> sub-map is good, he "<<heh.first<<" is unsplittable"<<std::endl;
        }else{
            std::cout<<" --> created bad sub-map, he "<<heh.first<<" is not unsplittable"<<std::endl;
        }
    }


    std::cout<<"... done"<<std::endl;
    std::cout<<" ----------"<<std::endl;

    return unsplittable;
}



bool Uncollapser::register_unsplit_in_subtriangle_map(const std::pair<HalfEdgeHandle, int>& heh,
                                                      const EdgeSplit& split){


    return false;

}

#endif


VertexHandle Uncollapser::uncollapse(const EdgeCollapse& collapse){

    //std::cout<<" -------------------------------------------------- "<<std::endl;
    //std::cout<<" =========================> UNDOING COLLAPSE ("<<collapse.from_vertex<<"-"<<collapse.to_vertex<<") -> "<<collapse.replacement_vertex<<std::endl;

    auto uncollapse_result_prop = mesh_.request_vertex_property<VertexHandle>("uncollapse_result");
    auto equatorial_disc_boundary_prop = mesh_.request_vertex_property<bool>("equatorial disc boundary");


    auto equatorial_disc = recover_equatorial_disc(collapse);


    /*std::cout<<" equatorial disc faces: ";
    for(auto ef: equatorial_disc){
        std::cout<<" "<<ef;
    }
    std::cout<<std::endl;*/


    if(!mark_equatorial_disc_boundary_vertices(equatorial_disc,
                                               equatorial_disc_boundary_prop)){
        std::cerr<<" AN ERROR OCCURED WHILE MARKING THE EQUATORIAL DISC BOUNDARY"<<std::endl;
        return VertexHandle(-1);
    }

    std::set<VertexHandle> north_hemispherical_vertices;


    std::queue<EdgeCollapse> to_uncollapse;
    to_uncollapse.push(collapse);


    if(!insert_spoke_secondary_uncollapses_in_queue(collapse,
                                                    to_uncollapse)){
        std::cerr<<" AN ERROR OCCURED WHILE INSERTING THE SPOKE SECONDARY UNCOLLAPSES IN THE QUEUE"<<std::endl;
        return VertexHandle(-1);
    }

    while(!to_uncollapse.empty()){

        auto current_collapse = to_uncollapse.front();
        to_uncollapse.pop();

        //if the vertex wasn't already uncollapse
        if(uncollapse_result_prop[current_collapse.replacement_vertex].idx() == -1){

            //if there's no from-vertex yet, we add a new one
            //this happens iff it's a secondary collapse
            if(current_collapse.from_vertex.idx() == -1){
                //std::cout<<" -> secondary collapse, adding new vertex"<<std::endl;
                current_collapse.from_vertex = mesh_.add_vertex({0,0,0});
            }


            auto equatorial_neighborhood = find_equatorial_neighborhood(current_collapse,
                                                                        equatorial_disc);

            auto result = uncollapse_helper(current_collapse,
                                            equatorial_neighborhood,
                                            uncollapse_result_prop);

            if(result.idx() == -1){
                return result;
            }

            //std::cout<<" CHECKING FOR SECONDARY uncollapseS..."<<std::endl;
            for(auto ef: equatorial_neighborhood){
                for(auto eq_neighbor: ef.get_vertices()){
                    //std::cout<<" - checking eq neighbor "<<eq_neighbor<<"...";
                    if(!equatorial_disc_boundary_prop[eq_neighbor]){
                        if(uncollapse_result_prop[eq_neighbor].idx() == -1){

                            to_uncollapse.push(make_secondary_collapse(eq_neighbor,
                                                                       collapse));
                            //std::cout<<" - added secondary uncollapse with to-vertex "<<eq_neighbor<<" to queue"<<std::endl;

                        }else{
                            //std::cout<<" ...which was already visited"<<std::endl;
                        }
                    }else{
                        //std::cout<<" ...which is not a replacement vertex"<<std::endl;
                    }
                }
            }
        }
    }


    if(!update_subface_map(collapse,
                           uncollapse_result_prop)){

        std::cerr<<" AN ERROR OCCURRED WHILE UPDATING THE SURFACE MAP."<<std::endl;
        return VertexHandle(-1);
    }


    /*std::cout<<" performed uncollapse,"
               " swapping back replacement_vertex "<<collapse.replacement_vertex<<
               " and to_vertex "<<collapse.to_vertex<<std::endl;*/


    //replace back the to-vertex
    sub_triangle_map_.replace_vertex_in_all_faces_containing_vertex(collapse.replacement_vertex,
                                                                    collapse.to_vertex);


    mesh_.swap_vertex_indices(collapse.to_vertex,
                              collapse.replacement_vertex);



    //std::cout<<" -- DONE UNCOLLAPSE "<<collapse.replacement_vertex<<" -> ("<<collapse.from_vertex<<"-"<<collapse.to_vertex<<")"<<std::endl;
    //std::cout<<" -------------------------------------------------- "<<std::endl;

    return collapse.from_vertex;
}


bool Uncollapser::extract_sorted_spoke_vertices(const EdgeCollapse& collapse,
                                                const int& spoke_index,
                                                std::vector<VertexHandle>& sorted_spoke_vertices) const{

    const auto& ev_count      = collapse.equatorial_vertices.size();
    const auto& center_vertex = collapse.replacement_vertex;
    const auto& previous_ev   = collapse.equatorial_vertices[(spoke_index + ev_count - 1) % ev_count];
    const auto& current_ev    = collapse.equatorial_vertices[spoke_index % ev_count];
    const auto& next_ev       = collapse.equatorial_vertices[(spoke_index + 1) % ev_count];

    /*std::cout<<" previous index : "<<((spoke_index + eq_vertices_count - 1) % eq_vertices_count)<<std::endl;
    std::cout<<" current index  : "<<(spoke_index % eq_vertices_count)<<std::endl;
    std::cout<<" next index     : "<<((spoke_index + 1) % eq_vertices_count)<<std::endl;

    std::cout<<" center vertex: "<<center_vertex<<std::endl;
    std::cout<<" previous_ev  : "<<previous_ev<<std::endl;
    std::cout<<" current_ev   : "<<current_ev<<std::endl;
    std::cout<<" next_ev      : "<<next_ev<<std::endl;
    std::cout<<" from vertex  : "<<collapse.from_vertex<<std::endl;*/

    SubTriangleFace spoke_triangle({collapse.from_vertex,
                                    center_vertex,
                                    current_ev});

    SubTriangleFace previous_triangle({center_vertex,
                                       previous_ev,
                                       current_ev});

    SubTriangleFace next_triangle({center_vertex,
                                   current_ev,
                                   next_ev});

    //std::cout<<" -> spoke triangle "<<spoke_triangle<<" between triangles "<<previous_triangle<<" and "<<next_triangle<<"..."<<std::endl;


    if(sub_triangle_map_.contains(spoke_triangle)){
        std::cout<<" ERROR - sub-face map already contains an entry for spoke triangle "<<spoke_triangle<<std::endl;
        return false;
    }

    std::set<VertexHandle> spoke_vertices;

    sub_triangle_map_.extract_subvertices_intersection(previous_triangle,
                                                       next_triangle,
                                                       spoke_vertices);

    if(!spoke_vertices.size()){
        std::cerr<<" ERROR - no common vertices between triangles "<<previous_triangle<<" and "<<next_triangle<<std::endl;
        return false;
    }

    //std::cout<<" - spoke vertices: "<<spoke_edge_vertices<<std::endl;

    if(spoke_vertices.size() > 2){

        auto visited = mesh_.request_vertex_property<bool>("visited");
        auto on_spoke_edge = mesh_.request_vertex_property<bool>("spoke_edge");
        for(auto spv: spoke_vertices){
            on_spoke_edge[spv] = true;
        }


        //std::cout<<" sorting them in the right order..."<<std::endl;
        //explore the eq edge vertices to get them in the right order
        VertexHandle current_vertex = collapse.replacement_vertex;
        while(current_vertex.idx() != -1){
            VertexHandle next_vertex(-1);
            visited[current_vertex] = true;
            sorted_spoke_vertices.push_back(current_vertex);

            for(auto out_he_it: mesh_.outgoing_halfedges(current_vertex)){
                auto to_vertex = mesh_.to_vertex_handle(out_he_it);

                if(!visited[to_vertex] &&
                        on_spoke_edge[to_vertex]){
                    next_vertex = to_vertex;
                }
            }
            current_vertex = next_vertex;
        }

        //std::cout<<" sorted spoke n°"<<spoke_index<<" vertices: "<<sorted_spoke_vertices<<std::endl;

        if(sorted_spoke_vertices.back() != current_ev){
            std::cerr<<" ERROR - last spoke vertex "<<sorted_spoke_vertices.back()<<
                       " is not the current equatorial vertex "<<current_ev<<std::endl;
            return false;
        }
    }
    return true;
}



bool Uncollapser::insert_spoke_secondary_uncollapses_in_queue(const EdgeCollapse& collapse,
                                                              std::queue<EdgeCollapse>& queue){

    //std::cout<<" ------"<<std::endl;
    //std::cout<<" inserting spoke secondary uncollapses in queue..."<<std::endl;

    std::vector<std::vector<VertexHandle>> all_sorted_spoke_vertices;

    const int& eq_vertices_count = collapse.equatorial_vertices.size();
    size_t max_spoke_vertices_count(0);

    for(int i(0); i<eq_vertices_count; i++){

        const auto& spoke_index = i;

        std::vector<VertexHandle> sorted_spoke_vertices;

        if(!extract_sorted_spoke_vertices(collapse,
                                          spoke_index,
                                          sorted_spoke_vertices)){
            std::cerr<<" AN ERROR OCCURED WHILE PRE-QUEUEING THE SPOKE VERTICES"<<std::endl;
            return false;
        }

        if(sorted_spoke_vertices.size()){
            all_sorted_spoke_vertices.push_back(sorted_spoke_vertices);

            max_spoke_vertices_count = std::max(max_spoke_vertices_count,
                                                sorted_spoke_vertices.size());
        }
    }

    if(max_spoke_vertices_count > 2){

        //now that we gathered all the spoke vertices we add them to the queue in a BFS-like manner
        //starting from 1 since the first one is the center vertex
        //and stopping before the last one, which is an equatorial vertex
        //and thus shouldn't be uncollapse
        for(size_t i(1); i<max_spoke_vertices_count-1; i++){
            //std::cout<<" - adding "<<i<<"-th spoke vertices..."<<std::endl;

            for(size_t j(0); j<all_sorted_spoke_vertices.size(); j++){
                if(all_sorted_spoke_vertices[j].size() &&
                        i<(all_sorted_spoke_vertices[j].size()-1)){

                    //std::cout<<" --> added vertex "<<sorted_spoke_vertices[j][i]<<" on spoke n°"<<j<<std::endl;
                    queue.push(make_secondary_collapse(all_sorted_spoke_vertices[j][i],
                                                       collapse));
                }
            }
        }
    }


    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------"<<std::endl;

    return true;
}




EdgeCollapse Uncollapser::make_secondary_collapse(const VertexHandle& secondary_to_vertex,
                                                  const EdgeCollapse& primary_collapse) const{

    EdgeCollapse secondary_collapse = primary_collapse;
    secondary_collapse.from_vertex = VertexHandle(-1);
    secondary_collapse.to_vertex = secondary_to_vertex;
    secondary_collapse.replacement_vertex = secondary_to_vertex;

    return secondary_collapse;
}




VertexHandle Uncollapser::uncollapse_helper(const EdgeCollapse&    collapse,
                                            const SubFaceSet&      equatorial_neighborhood,
                                            VertexPropertyT<VertexHandle>& uncollapse_result_prop){

    //std::cout<<" -------------------------------------------------- "<<std::endl;
    //std::cout<<" =========================> UNDOING COLLAPSE "<<collapse<<std::endl;



    //vertex being uncollapse
    const VertexHandle center_vertex(collapse.replacement_vertex);

    if(uncollapse_result_prop[center_vertex].idx() != -1){
        //std::cout<<" ==> ALREADY DONE, SKIPPING"<<std::endl;
        //std::cout<<" -------------------------------------------------- "<<std::endl;
        return uncollapse_result_prop[center_vertex];
    }

    uncollapse_result_prop[center_vertex] = collapse.from_vertex;
    //std::cout<<" --> registered secondary uncollapse vertex "<<center_vertex<<
    //           " as "<<uncollapse_result_prop[center_vertex]<<std::endl;



    /*std::cout<<" equatorial neighborhood: ";
    for(auto ef: equatorial_neighborhood){
        std::cout<<" "<<ef;
    }
    std::cout<<std::endl;*/



    //recover faces originally opposite to from_vertex
    auto north_hemispherical_halffaces = PEHelpers::find_hemispherical_faces(mesh_,
                                                                             center_vertex,
                                                                             equatorial_neighborhood,
                                                                             false);


    //std::cout<<" - recovered north hfs: "<<north_hemispherical_halffaces.size()<<std::endl;
    /*for(auto hf: north_hemispherical_halffaces){
        std::cout<<" -- ";
        for(auto v: hf){
            std::cout<<" "<<v;
        }
        std::cout<<std::endl;
    }*/


    //recover faces originally opposite to to_vertex
    auto south_hemispherical_halffaces = PEHelpers::find_hemispherical_faces(mesh_,
                                                                             center_vertex,
                                                                             equatorial_neighborhood,
                                                                             true);


    //std::cout<<" - recovered south hfs: "<<south_hemispherical_halffaces.size()<<std::endl;
    /*for(auto hf: south_hemispherical_halffaces){
        std::cout<<" -- ";
        for(auto v: hf){
            std::cout<<" "<<v;
        }
        std::cout<<std::endl;
    }*/


    //the is_boundary() check is there because uncollapse done on the boundary might have
    //empty hemispheres
    if(!mesh_.is_boundary(center_vertex)){
        if(!north_hemispherical_halffaces.size() ||
                !south_hemispherical_halffaces.size()){

            std::cerr<<" ERROR - couldn't recover opposite halffaces to center vertices"<<std::endl;
            std::cout<<" opposite hfs to from_vertex: "<<north_hemispherical_halffaces.size()<<std::endl;
            std::cout<<" opposite hfs to to_vertex  : "<<south_hemispherical_halffaces.size()<<std::endl;

            return VertexHandle(-1);
        }
    }


    ACG::Vec3d new_from_vertex_position;
    ACG::Vec3d new_to_vertex_position;


    if(!GeometryHelper::find_new_vertices_positions(mesh_,
                                                    center_vertex,
                                                    collapse.equatorial_vertices,
                                                    north_hemispherical_halffaces,
                                                    south_hemispherical_halffaces,
                                                    new_from_vertex_position,
                                                    new_to_vertex_position)){
        std::cerr<<" ERROR - couldn't find a position for from-vertex "<<collapse.from_vertex<<std::endl;
        return VertexHandle(-1);
    }



    //delete existing vertex to be uncollapse to create a void
    mesh_.delete_vertex(center_vertex);

    //replace the uncollapse vertex with a new one
    auto new_to_vertex = mesh_.add_vertex(new_to_vertex_position);

    if(center_vertex.idx() != new_to_vertex.idx()){
        uncollapse_result_prop[new_to_vertex] = uncollapse_result_prop[center_vertex];

        mesh_.swap_vertex_indices(new_to_vertex, center_vertex);
        new_to_vertex = center_vertex;

    }

    //new vertex
    auto new_from_vertex = mesh_.add_vertex(new_from_vertex_position);
    if(collapse.from_vertex.idx() != new_from_vertex.idx()){
        mesh_.swap_vertex_indices(new_from_vertex, collapse.from_vertex);
        new_from_vertex = collapse.from_vertex;
    }




    //std::cout<<" - adding north hemisphere cells..."<<std::endl;
    //add all cells originally connected to v_new
    if(!add_cells_from_vertex_and_opposite_halffaces(new_from_vertex,
                                                     north_hemispherical_halffaces)){
        return VertexHandle(-1);
    }

    //std::cout<<" - adding south hemisphere cells..."<<std::endl;
    //add all cells originally connected to v_uncollapse
    if(!add_cells_from_vertex_and_opposite_halffaces(new_to_vertex,
                                                     south_hemispherical_halffaces)){
        return VertexHandle(-1);
    }


    std::set<VertexHandle> equatorial_neighbor_vertices;

    auto equatorial_neighbor_prop = mesh_.request_vertex_property<bool>("equatorial neighbor");
    //add the equatorial cells and gather the equatorial neighbor vertices
    for(auto ef: equatorial_neighborhood){
        auto vertices = ef.opposite_face().get_vertices();
        //we also mark the neighbors as such for the next operation
        for(auto en: vertices){
            equatorial_neighbor_prop[en] = true;
        }

        vertices.push_back(new_from_vertex);

        auto ch_new = mesh_.add_cell(vertices, true);
        if(ch_new.idx() == -1){
            std::cerr<<" ERROR - couldn't create new equatorial cell "<<ef<<" + "<<new_from_vertex<<std::endl;
            return VertexHandle(-1);
        }
        //std::cout<<" added new equatorial cell "<<ef<<" + "<<v_new<<std::endl;
    }

    //std::cout<<" -- DONE UNCOLLAPSE "<<collapse.replacement_vertex<<" -> ("<<collapse.from_vertex<<"-"<<collapse.to_vertex<<")"<<std::endl;
    //std::cout<<" -------------------------------------------------- "<<std::endl;


    return new_from_vertex;
}





bool Uncollapser::mark_equatorial_disc_boundary_vertices(const EquatorialDisc& equatorial_disc,
                                                         VertexPropertyT<bool>& equatorial_disc_boundary_prop) const{

    //std::cout<<" ------------"<<std::endl;
    //std::cout<<" marking equatorial disc boundary vertices..."<<std::endl;


    TriMesh equatorial_disc_mesh;

    //add the vertices to the tri-mesh
    auto ovm_to_om_vertex_prop = mesh_.request_vertex_property<OpenMesh::SmartVertexHandle>("ovm to om");
    OpenMesh::VPropHandleT<VertexHandle> om_to_ovm_prop;
    equatorial_disc_mesh.add_property(om_to_ovm_prop);
    auto visited = mesh_.request_vertex_property<bool>("visited");
    for(auto subface: equatorial_disc){
        //std::cout<<" converting tet-face "<<subface<<std::endl;

        std::vector<OpenMesh::SmartVertexHandle> tri_mesh_face_vertices;
        for(auto v: subface.get_vertices()){
            if(v.idx() == -1){
                std::cerr<<" ERROR - found invalid vertex in subface "<<subface<<std::endl;
                return false;
            }
            if(!visited[v]){
                visited[v] = true;
                ovm_to_om_vertex_prop[v] = equatorial_disc_mesh.add_vertex(mesh_.vertex(v));
                equatorial_disc_mesh.property(om_to_ovm_prop, ovm_to_om_vertex_prop[v]) = v;
            }
            tri_mesh_face_vertices.push_back(ovm_to_om_vertex_prop[v]);
        }


        //add the face to the surface mesh
        auto tri_face = equatorial_disc_mesh.add_face(tri_mesh_face_vertices);

        if(tri_face.idx() == -1){
            std::cerr<<" ERROR - could not add subface "<<subface<<" to ED tri-mesh"<<std::endl;
            return false;
        }
    }


    //std::cout<<" --> boundary vertices: ";

    //then mark all the boundary vertices in the TriMesh
    for(auto v: equatorial_disc_mesh.vertices()){
        if(equatorial_disc_mesh.is_boundary(v)){
            const auto& tet_mesh_v = equatorial_disc_mesh.property(om_to_ovm_prop, v);
            //std::cout<<" "<<tet_mesh_v;
            equatorial_disc_boundary_prop[tet_mesh_v] = true;
        }
    }
    //std::cout<<std::endl;



    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------------"<<std::endl;

    return true;
}



bool Uncollapser::update_subface_map(const EdgeCollapse& collapse,
                                     const VertexPropertyT<VertexHandle>& uncollapse_result_prop){

    //std::cout<<" =================================================================================== "<<std::endl;
    //std::cout<<" UPDATING SUBFACE MAP'S NORTHERN EQUATORIAL DISC, SPOKES AND NORTHERN HEMISPHERE..."<<std::endl;
    /*std::cout<<" equatorial vertices: ";
    for(auto v: collapse.equatorial_vertices){
        std::cout<<" "<<v;
    }
    std::cout<<std::endl;*/


    //std::cout<<" ======================== EQUATORIAL TRIANGLES AND SPOKES"<<std::endl;

    if(!update_equatorial_triangles_and_spokes(collapse, uncollapse_result_prop)){
        std::cerr<<" AN ERROR OCCURED WHILE UPDATING THE EQUATORIAL TRIANGLES AND SPOKES, STOPPING PROCESS"<<std::endl;
        return false;
    }


    //std::cout<<" ======================== NORTH HEMISPHERE TRIANGLES"<<std::endl;

    if(!update_north_hemisphere_triangles(collapse,
                                          uncollapse_result_prop)){
        std::cerr<<" AN ERROR OCCURED WHILE uncollapseTING A NORTHERN TRIANGLE, STOPPING PROCESS"<<std::endl;
        return false;
    }


    //std::cout<<" =====  north triangles connected to two (non-consecutive) eq vertices..."<<std::endl;

    if(!update_bypass_triangles(collapse, uncollapse_result_prop)){
        std::cerr<<" AN ERROR OCCURED WHILE UPDATING THE BYPASS TRIANGLES, STOPPING PROCESS"<<std::endl;
        return false;
    }

    //std::cout<<" ... DONE"<<std::endl;
    //std::cout<<" =================================================================================== "<<std::endl;

    return true;
}




bool Uncollapser::update_equatorial_triangles_and_spokes(const EdgeCollapse& collapse,
                                                         const VertexPropertyT<VertexHandle>& uncollapse_result_prop){


    //std::cout<<" ------"<<std::endl;
    //std::cout<<" updating equatorial triangles and spokes..."<<std::endl;


    const auto& from_vertex = collapse.from_vertex;
    const auto& replacement_vertex = collapse.replacement_vertex;

    const auto& ev_count = collapse.equatorial_vertices.size();

    //go through each equatorial triangle
    for(size_t i(0); i<collapse.equatorial_vertices.size(); i++){
        auto current_ev = collapse.equatorial_vertices[i];
        auto next_ev = collapse.equatorial_vertices[(i+1) % ev_count];


        auto southern_eq_triangle = SubTriangleFace({replacement_vertex,
                                                     current_ev,
                                                     next_ev});

        //std::cout<<" ======= HANDLING NEW TRIANGLE "<<southern_eq_triangle<<std::endl;


        if(sub_triangle_map_.contains(southern_eq_triangle)){

            //uncollapse all the "northern" side triangles of the equator as a mirror of the southern side
            auto northern_eq_triangle = SubTriangleFace({from_vertex,
                                                         current_ev,
                                                         next_ev});

            if(!copy_subfaces_and_replace_vertices_with_uncollapse_result(southern_eq_triangle,
                                                                          northern_eq_triangle,
                                                                          uncollapse_result_prop)){
                std::cerr<<" ERROR while updating surface map"<<std::endl;
                return false;
            }


            //std::cout<<" ===== SPOKE"<<std::endl;

            //uncollapse the "spoke" triangles and update the
            //north hemisphere triangles containing the equatorial edge to contain the
            //north equatorial edge
            if(!partition_spoke_triangle(collapse, i,
                                     uncollapse_result_prop)){
                std::cout<<" AN ERROR OCCURED WHILE uncollapseTING A SPOKE TRIANGLE, STOPPING PROCESS"<<std::endl;
                return false;
            }


        }else{
            //std::cout<<" this triangle wasn't uncollapse, skipping"<<std::endl;
        }

        //std::cout<<" ======= DONE WITH TRIANGLE"<<std::endl;

    }


    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------"<<std::endl;

    return true;
}




bool Uncollapser::update_bypass_triangles(const EdgeCollapse& collapse,
                                          const VertexPropertyT<VertexHandle>& uncollapse_result_prop){

    //std::cout<<" ------"<<std::endl;
    //std::cout<<" updating bypass triangles..."<<std::endl;


    std::vector<SubTriangleFace> bypass_triangles;

    const auto& from_vertex = collapse.from_vertex;
    const auto& replacement_vertex = collapse.replacement_vertex;

    auto eq_vertex_index_prop = mesh_.request_vertex_property<int>();

    int i(0);
    for(auto eq_v: collapse.equatorial_vertices){
        //i+1 so zero means it's not an EV
        eq_vertex_index_prop[eq_v] = i+1;
        i++;
    }

    const int& ev_count = collapse.equatorial_vertices.size();
    //std::cout<<" EV count: "<<ev_count<<std::endl;

    for(auto north_hf: collapse.north_hemispherical_hfs){
        //std::cout<<"----------------------------"<<std::endl;
        //std::cout<<" - checking north hf "<<north_hf<<std::endl;

        std::vector<std::pair<VertexHandle, VertexHandle>> bypass_edges;

        bool found_bypass_edge(false);
        for(int i(0); i<3; i++){
            //std::cout<<"------"<<std::endl;
            int current_index = i;
            int next_index    = (i+1) % 3;
            //std::cout<<" -- current index: "<<current_index<<", next index: "<<next_index<<std::endl;

            /*auto current_vertex = north_hf[current_index];
            auto next_vertex = north_hf[next_index];
            std::cout<<" -- current vertex: "<<current_vertex<<", next vertex: "<<next_vertex<<std::endl;*/

            int current_ev_index = eq_vertex_index_prop[north_hf[current_index]];
            int next_ev_index = eq_vertex_index_prop[north_hf[next_index]];
            //std::cout<<" -- current ev index: "<<current_ev_index<<", next ev index: "<<next_ev_index<<std::endl;

            //if the difference between the two indices is greater than one,
            //AND not equal to the number of EVs minus 1
            //(in the case where current = last EV and next = first EV)
            //then it means the halfface contains a bypass edge
            if(std::abs(current_ev_index - next_ev_index) > 1 &&
                    std::abs(current_ev_index - next_ev_index) != (ev_count-1)){
                found_bypass_edge = true;

                bypass_edges.push_back({north_hf[current_index], north_hf[next_index]});

                //std::cout<<" ---> found bypass edge ("<<bypass_edges.back().first<<"-"<<bypass_edges.back().second<<")"<<std::endl;
            }

        }


        for(auto bypass_edge: bypass_edges){

            SubTriangleFace north_bypass_triangle({bypass_edge.first,
                                                   bypass_edge.second,
                                                   replacement_vertex});
            //std::cout<<" -- checking bypass triangle "<<north_bypass_triangle<<std::endl;

            if(sub_triangle_map_.contains(north_bypass_triangle)){

                //std::cout<<" ---> bypass triangle was uncollapse, updating"<<std::endl;
                SubFaceSet subfaces;
                sub_triangle_map_.extract_all_subfaces(north_bypass_triangle,
                                                       subfaces);


                bypass_triangles.push_back(north_bypass_triangle);

                SubTriangleFace updated_north_bypass_triangle({bypass_edge.first,
                                                               bypass_edge.second,
                                                               from_vertex});

                if(!copy_subfaces_and_replace_vertices_with_uncollapse_result(north_bypass_triangle,
                                                                         updated_north_bypass_triangle,
                                                                         uncollapse_result_prop)){
                    std::cerr<<" ERROR while handling 'bypass' north hemisphere triangles "<<std::endl;
                    return false;
                }

                //erase the old entry now that it's updated
                //sub_triangle_map_.erase(north_bypass_triangle);

            }else{
                //std::cout<<" ---> bypass triangle wasn't uncollapse, skipping"<<std::endl;
            }
        }
    }


    //erase the now-obsolete bypass triangles
    for(auto bypass_triangle: bypass_triangles){
        //std::cout<<bypass_triangle<<" ";
        sub_triangle_map_.erase(bypass_triangle);
    }
    //std::cout<<std::endl;


    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------"<<std::endl;

    return true;
}




bool Uncollapser::update_north_hemisphere_triangles(const EdgeCollapse& collapse,
                                                    const VertexPropertyT<VertexHandle>& uncollapse_result_prop){


    //std::cout<<" ------"<<std::endl;
    //std::cout<<" updating north hemisphere triangles..."<<std::endl;

    const auto& center_vertex = collapse.replacement_vertex;

    //mark the equatorial vertices and northern vertices as such
    auto equatorial_vertex_prop = mesh_.request_vertex_property<bool>();
    auto north_hemispherical_vertex_prop = mesh_.request_vertex_property<bool>();

    auto visited = mesh_.request_vertex_property<bool>();


    for(auto eq_v: collapse.equatorial_vertices){
        equatorial_vertex_prop[eq_v] = true;
    }

    std::set<VertexHandle> north_hemispherical_vertices;
    for(auto north_hf: collapse.north_hemispherical_hfs){
        //std::cout<<" - gathering north hf "<<north_hf<<" vertices "<<std::endl;
        for(auto north_v: north_hf){
            if(!equatorial_vertex_prop[north_v]){
                north_hemispherical_vertex_prop[north_v] = true;
                north_hemispherical_vertices.insert(north_v);
            }
        }
    }

    //std::cout<<" north hemispherical vertices: "<<north_hemispherical_vertices<<std::endl;

    //handle the north vertices
    for(auto north_v: north_hemispherical_vertices){
        //std::cout<<" ---------"<<std::endl;
        //std::cout<<" checking north hemispherical vertex "<<north_v<<std::endl;

        visited[north_v] = true;

        //check for triangles connecting an equatorial vertex,
        //the center vertex and a north vertex
        for(auto eq_v: collapse.equatorial_vertices){
            SubTriangleFace north_triangle({north_v,
                                            center_vertex,
                                            eq_v});

            //std::cout<<" -- checking north hemisphere triangle "<<north_triangle<<std::endl;
            if(sub_triangle_map_.contains(north_triangle)){
                //std::cout<<" ---> in subface map, uncollapseting..."<<std::endl;


                SubTriangleFace north_triangle_upper_half({north_v,
                                                           collapse.from_vertex,
                                                           eq_v});

                //replace vertices on the equatorial edge with the new ones
                if(!copy_subfaces_and_replace_vertices_with_uncollapse_result(north_triangle,
                                                                         north_triangle_upper_half,
                                                                         uncollapse_result_prop)){
                    std::cerr<<" ERROR while handling north hemisphere triangles"<<std::endl;
                    return false;
                }

                //std::cout<<" --- north triangle upper half: "<<north_triangle_upper_half<<std::endl;


                //erase the current entry for the north triangle
                sub_triangle_map_.erase(north_triangle);




            }else{
                //std::cout<<" --> not in subface map, skipping"<<std::endl;
            }
        }

        //and then check for triangles connecting the center vertex and
        //two neighboring north vertices
        for(auto other_north_v: north_hemispherical_vertices){
            if(!visited[other_north_v]){
                //std::cout<<" - checking other north vertex "<<other_north_v<<std::endl;
                SubTriangleFace north_triangle({north_v,
                                                center_vertex,
                                                other_north_v});

                //std::cout<<" -- checking north hemisphere triangle "<<north_triangle<<std::endl;
                if(sub_triangle_map_.contains(north_triangle)){
                    //std::cout<<" ---> in subface map, uncollapseting..."<<std::endl;

                    SubTriangleFace updated_north_triangle({north_v,
                                                            collapse.from_vertex,
                                                            other_north_v});

                    //replace vertices on the equatorial edge with the new ones
                    if(!copy_subfaces_and_replace_vertices_with_uncollapse_result(north_triangle,
                                                                             updated_north_triangle,
                                                                             uncollapse_result_prop)){
                        std::cerr<<" ERROR while handling north hemisphere triangles"<<std::endl;
                        return false;
                    }

                    //std::cout<<" --- updated north triangle: "<<updated_north_triangle<<std::endl;


                    //erase the current entry for the north triangle
                    sub_triangle_map_.erase(north_triangle);

                }else{
                    //std::cout<<" --> not in subface map, skipping"<<std::endl;
                }
            }
        }
    }

    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------"<<std::endl;

    return true;
}




bool Uncollapser::partition_spoke_triangle(const EdgeCollapse& collapse,
                                       int spoke_index,
                                       const VertexPropertyT<VertexHandle>& uncollapse_result_prop){


    //std::cout<<" ------"<<std::endl;
    //std::cout<<" uncollapseting spoke n°"<<spoke_index<<"..."<<std::endl;

    const auto& ev_count      = collapse.equatorial_vertices.size();
    const auto& center_vertex = collapse.replacement_vertex;
    const auto& current_ev    = collapse.equatorial_vertices[spoke_index % ev_count];

    /*std::cout<<" previous index : "<<((spoke_index + eq_vertices_count - 1) % eq_vertices_count)<<std::endl;
    std::cout<<" current index  : "<<(spoke_index % eq_vertices_count)<<std::endl;
    std::cout<<" next index     : "<<((spoke_index + 1) % eq_vertices_count)<<std::endl;

    std::cout<<" center vertex: "<<center_vertex<<std::endl;
    std::cout<<" previous_ev  : "<<previous_ev<<std::endl;
    std::cout<<" current_ev   : "<<current_ev<<std::endl;
    std::cout<<" next_ev      : "<<next_ev<<std::endl;
    std::cout<<" from vertex  : "<<collapse.from_vertex<<std::endl;*/


    std::vector<VertexHandle> sorted_spoke_edge_vertices;

    if(!extract_sorted_spoke_vertices(collapse,
                                      spoke_index,
                                      sorted_spoke_edge_vertices)){
        std::cerr<<" ERROR - could not extract the sorted vertices of spoke n°"<<spoke_index<<std::endl;
        return false;
    }

    if(!sorted_spoke_edge_vertices.size()){
        //std::cout<<" spoke wasn't uncollapse"<<std::endl;
        return true;
    }

    //std::cout<<" sorted spoke vertices: "<<sorted_spoke_edge_vertices<<std::endl;

    /*decompose the equatorial spoke face as follows:
     *  0'
     *  |``-1'..2'
     *  |  /|  /|`.
     *  |/  |/  |  `.
     *  0'''1'''2''''3
     * where
     *  assuming spoke_triangle = {0', 0, 3}
     *  {0,1,2,3} = sorted_spoke_edge_vertices
     *  0  = center_vertex
     *  0' = from_vertex (assuming collapse of edge (from_vertex-to_vertex)
     *  1',2' = secondary uncollapse from_vertices
     *  3  = current equatorial vertex (current_ev)
     * this is done by taking the successive "squares" {i', i, (i+1), (i+1)'}
     * and decomposing them as the triangles {i', i, (i+1)'} and {i, (i+1), (i+1)'}
     * that for each i from 0 to k-2, k being the number of vertices on the
     * equatorial edge 0-3 (including both) */
    std::vector<SubTriangleFace> subfaces;
    for(size_t i(0); i<sorted_spoke_edge_vertices.size() - 2; i++){
        const auto& eev_i  = sorted_spoke_edge_vertices[i];
        const auto& eev_i1 = sorted_spoke_edge_vertices[i+1];

        /*std::cout<<" --- i="<<i<<std::endl;
        std::cout<<" eev_i'   = "<<uncollapse_result_prop[eev_i]<<std::endl;
        std::cout<<" eev_i    = "<<eev_i<<std::endl;
        std::cout<<" eev_i+1  = "<<eev_i1<<std::endl;
        std::cout<<" eev_i+1' = "<<uncollapse_result_prop[eev_i1]<<std::endl;*/

        SubTriangleFace first_triangle({uncollapse_result_prop[eev_i],
                                        eev_i,
                                        uncollapse_result_prop[eev_i1]});

        SubTriangleFace second_triangle({eev_i,
                                         eev_i1,
                                         uncollapse_result_prop[eev_i1]});

        //std::cout<<" adding triangles "<<first_triangle<<" and "<<second_triangle<<std::endl;

        subfaces.push_back(first_triangle);
        subfaces.push_back(second_triangle);
    }

    //add the "tip" triangle (2', 2, 3)
    auto last_vertex_before_eq_v = sorted_spoke_edge_vertices[sorted_spoke_edge_vertices.size()-2];
    auto last_vertex_before_eq_v_prime = uncollapse_result_prop[last_vertex_before_eq_v];

    SubTriangleFace tip_triangle({last_vertex_before_eq_v_prime,
                                  last_vertex_before_eq_v,
                                  current_ev});

    //std::cout<<" tip triangle: "<<tip_triangle<<std::endl;
    subfaces.push_back(tip_triangle);


    SubTriangleFace spoke_triangle({collapse.from_vertex,
                                    center_vertex,
                                    current_ev});

    /*std::cout<<" spoke subfaces: ";
    for(auto subface: subfaces){
        std::cout<<subface<<" ";
    }
    std::cout<<std::endl;*/

    sub_triangle_map_.insert(spoke_triangle, subfaces);



    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------"<<std::endl;

    return true;
}







bool Uncollapser::copy_subfaces_and_replace_vertices_with_uncollapse_result(const SubTriangleFace& source_face,
                                                                       const SubTriangleFace& target_face,
                                                                       const VertexPropertyT<VertexHandle>& uncollapse_result_prop){

    /*std::cout<<" ------"<<std::endl;
    std::cout<<" copying all subfaces of triangle "<<source_face<<
               " to triangle "<<target_face<<
               "..."<<std::endl;*/


    SubFaceSet source_subfaces;
    sub_triangle_map_.extract_all_subfaces(source_face,
                                           source_subfaces);

    sub_triangle_map_.erase(target_face);

    /*std::cout<<" source subfaces: ";
    for(auto subface: source_subfaces){
        std::cout<<subface<<" ";
    }
    std::cout<<std::endl;*/

    for(auto source_subface: source_subfaces){
        //std::cout<<" - updating source subface "<<subface<<"... "<<std::endl;;
        auto source_subface_vertices = source_subface.get_vertices();
        for(auto subvertex: source_subface_vertices){
            //std::cout<<" -- subvertex: "<<subvertex<<", secondary uncollapse vertex: "<<uncollapse_result_prop[subvertex]<<std::endl;

            //if the vertex of the subface was part of a secondary uncollapse
            //we replace it with the vertex it was uncollapse to
            if(uncollapse_result_prop[subvertex].idx() != -1
                    && !source_subface.contains_vertex(uncollapse_result_prop[subvertex])){

                source_subface.replace_vertex(subvertex,
                                              uncollapse_result_prop[subvertex]);
            }
        }

        sub_triangle_map_.insert(target_face, source_subface);

    }

    /*std::cout<<" ----"<<std::endl;
    std::cout<<" source subfaces : ";
    for(auto source_subface: source_subfaces){
        std::cout<< source_subface<<" ";
    }
    std::cout<<std::endl;

    std::cout<<" target subfaces : ";
    for(auto target_subface: sub_triangle_map_.find(target_face)->second){
        std::cout<<target_subface<<" ";
    }
    std::cout<<std::endl;*/


    /*std::cout<<" ...done"<<std::endl;
    std::cout<<" ------"<<std::endl;*/

    return true;

}




SubFaceSet Uncollapser::find_equatorial_neighborhood(const EdgeCollapse& collapse,
                                                                    const EquatorialDisc& equatorial_disc) const{
    //std::cout<<" ----------------------------"<<std::endl;
    //std::cout<<" recovering equatorial neighborhood..."<<std::endl;

    SubFaceSet equatorial_neighborhood;

    const auto& center_vertex = collapse.replacement_vertex;

    for(auto vhf_it = mesh_.vhf_iter(center_vertex); vhf_it.valid(); vhf_it++){
        std::vector<VertexHandle> hf_vertices;
        for(auto hfv_it = mesh_.hfv_iter(*vhf_it); hfv_it.valid(); hfv_it++){
            hf_vertices.push_back(*hfv_it);
        }

        SubTriangleFace face(hf_vertices);

        //std::cout<<" - checking face "<<face<<std::endl;

        if(equatorial_disc.find(face) != equatorial_disc.end()){
            //std::cout<<" --> found eq face "<<face<<std::endl;
            equatorial_neighborhood.insert(face);
        }
    }

    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ----------------------------"<<std::endl;
    return equatorial_neighborhood;
}




std::set<VertexHandle> Uncollapser::recover_subface_set_vertices(const SubFaceSet& subface_set){

    std::set<VertexHandle> subface_set_vertices;

    for(auto ef: subface_set){
        for(auto v: ef.get_vertices()){
            subface_set_vertices.insert(v);
        }
    }

    return subface_set_vertices;
}


SubFaceSet Uncollapser::recover_equatorial_disc(const EdgeCollapse& collapse){

    /*std::cout<<" ----------------------------"<<std::endl;
    std::cout<<" recovering equatorial disc subfaces..."<<std::endl;
    std::cout<<" equatorial vertices: ";
    for(auto v: collapse.equatorial_vertices){
        std::cout<<" "<<v;
    }
    std::cout<<std::endl;*/

    SubFaceSet equatorial_disc_faces;

    const auto& center_vertex = collapse.replacement_vertex;

    for(size_t i(0); i<collapse.equatorial_vertices.size(); i++){
        auto ev_i = collapse.equatorial_vertices[i];
        auto ev_i1 = collapse.equatorial_vertices[(i+1) % collapse.equatorial_vertices.size()];

        SubTriangleFace equatorial_triangle({center_vertex, ev_i, ev_i1});
        //std::cout<<" - gathering eq triangle "<<equatorial_triangle<<" subfaces"<<std::endl;

        sub_triangle_map_.extract_all_subfaces(equatorial_triangle,
                                               equatorial_disc_faces);
    }

    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ----------------------------"<<std::endl;

    return equatorial_disc_faces;
}





bool Uncollapser::add_cells_from_vertex_and_opposite_halffaces(const OpenVolumeMesh::VertexHandle& center_vertex,
                                                               const std::vector<std::vector<OpenVolumeMesh::VertexHandle>>& opposite_halffaces){

    for(auto hf: opposite_halffaces){

        //and add the same but connected to v_new
        hf.push_back(center_vertex);

        /*std::cout<<" adding cell: ";
        for(auto v: hf){
            std::cout<<" "<<v;
        }
        std::cout<<std::endl;*/

        auto ch_new = mesh_.add_cell(hf, true);
        if(ch_new.idx() == -1){
            std::cout<<" ERROR - couldn't create new cell "<<hf<<std::endl;
            return false;
        }

        /*std::cout<<" added cell "<<ch_new<<": ";
        for(auto cv_it = mesh_.cv_iter(ch_new); cv_it.valid(); cv_it++){
            std::cout<<*cv_it<<" ";
        }
        std::cout<<std::endl;*/
    }

    return true;
}



void Uncollapser::print_subfaces(std::vector<VertexHandle> face){
    SubFaceSet subfaces;
    sub_triangle_map_.extract_all_subfaces(SubTriangleFace(face),
                                           subfaces);
    std::cout<<" subfaces of face ("<<face<<") : "<<std::endl;
    for(auto subface: subfaces){
        std::cout<<subface<<std::endl;
    }
    std::cout<<std::endl;
}



}
