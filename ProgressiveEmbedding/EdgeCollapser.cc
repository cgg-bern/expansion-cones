#include "EdgeCollapser.hh"


using namespace std;


using Point = TetrahedralMesh::PointT;

namespace OpenVolumeMesh{

EdgeCollapser::EdgeCollapser(TetrahedralMesh& mesh) :
    mesh_(mesh),
    link_condition_prop_(mesh_.request_edge_property<bool>("link_condition")),
    topo_helper_(mesh){

    initial_mesh_ = mesh_;
}


#define DEBUG_EDGE_COLLAPSER true

#define DEBUG_PRINT(cout_arg) if(DEBUG_EDGE_COLLAPSER){std::cout<<cout_arg;}

#define MAX_EDGE_VALENCE_FOR_HISTOGRAM 10


int EdgeCollapser::twoStepCollapse(TetrahedralMesh& mesh){

    EdgeCollapser edge_collapser(mesh);

    return edge_collapser.twoStepCollapse(false);
}



int EdgeCollapser::threeStepCollapse(TetrahedralMesh& mesh){

    EdgeCollapser edge_collapser(mesh);

    return edge_collapser.threeStepCollapse();
}




bool EdgeCollapser::split(const EdgeSplit& split){
    auto halfedge_to_split = topo_helper_.halfedge_exists(split.from_vertex,
                                                          split.to_vertex);
    if(halfedge_to_split.idx() == -1){
        std::cerr<<" ERROR - requesting split "<<split<<
                   " but halfedge ("<<split.from_vertex<<"-"<<split.to_vertex<<")"<<
                   " doesn't exist"<<std::endl;
        return false;
    }

    auto resplit = splitEdge(halfedge_to_split);

    if(resplit.middle_vertex != split.middle_vertex){
        std::cerr<<" ERROR - requested split is "<<split<<
                   " but new split is "<<resplit<<std::endl;
        return false;
    }

    return true;
}


bool EdgeCollapser::collapse(const EdgeCollapse& collapse){

    auto halfedge_to_collapse = topo_helper_.halfedge_exists(collapse.from_vertex,
                                                             collapse.to_vertex);
    if(halfedge_to_collapse.idx() == -1){
        std::cerr<<" ERROR - requesting collapse "<<collapse<<
                   " but halfedge ("<<collapse.from_vertex<<"-"<<collapse.to_vertex<<")"<<
                   " doesn't exist"<<std::endl;
        return false;
    }

    auto recollapse = collapseEdge(halfedge_to_collapse);

    if(recollapse.replacement_vertex != collapse.replacement_vertex){
        std::cerr<<" ERROR - requested collapse is "<<collapse<<
                   " but new collapse is "<<recollapse<<std::endl;
        return false;
    }

    return true;
}



COLLAPSE_RESULT_STATUS EdgeCollapser::twoStepCollapse(bool random_order){

    std::cout<<" ============================================================== "<<std::endl;

    std::cout<<" RUNNING FULL TWO-STEP EDGE COLLAPSE.... "<<std::endl;


    auto bad_tets = BadTetFinder::findBadTets(mesh_);
    std::cout<<" - initial bad tets : "<<std::endl;
    std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
    std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS"<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;

    COLLAPSE_RESULT_STATUS overall_result(COLLAPSE_FAILURE);

    COLLAPSE_RESULT_STATUS full_collapse(COLLAPSE_SUCCESS);
    COLLAPSE_RESULT_STATUS experimental_collapse(COLLAPSE_SUCCESS);

    int iteration_count(0);

    auto start_time = std::chrono::high_resolution_clock::now();

    while(overall_result != COLLAPSE_ERROR &&
          (full_collapse == COLLAPSE_SUCCESS || experimental_collapse == COLLAPSE_SUCCESS) &&
          iteration_count < 50){

        std::cout<<" ==================================== ITERATION "<<iteration_count<<std::endl;


        int current_deg_count(0);


        if(random_order){
            full_collapse = collapseAllCollapsibleEdgesInRandomOrder();

        }else{
            current_deg_count = 0;
            full_collapse = collapseAllCollapsibleEdges();
        }


        if(full_collapse == COLLAPSE_SUCCESS){
            std::cout<<" --> COLLAPSIBLE COLLAPSES SUCCESFUL"<<std::endl;
        }else if(full_collapse == COLLAPSE_ERROR){
            std::cout<<" ERROR - STOPPING TWO-STEP COLLAPSE"<<std::endl;
            overall_result = COLLAPSE_ERROR;
            break;
        }


        experimental_collapse = breakTriTets(random_order);
        //experimental_collapse = splitAndCollapseCycles(MIN_VALENCE_SELECTION);

        if(experimental_collapse == COLLAPSE_SUCCESS){
            std::cout<<" --> UNCOLLAPSIBLE CONFIGS BREAK SUCCESFUL"<<std::endl;
        }else if(experimental_collapse == COLLAPSE_ERROR){
            std::cout<<" ERROR - STOPPING TWO-STEP COLLAPSE"<<std::endl;
            overall_result = COLLAPSE_ERROR;
            break;
        }

        iteration_count++;

    }


    auto end_time = std::chrono::high_resolution_clock::now();

    bool only_one_resulting_interior_vertex(false);

    int remaining_interior_vertices_count = topo_helper_.interiorVerticesCount();

    if(overall_result == COLLAPSE_ERROR){
        std::cout<<" ERROR - SOMETHING BAD HAPPENED DURING COLLAPSE AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" ===========> REMAINING INTERIOR VERTICES: "<<remaining_interior_vertices_count<<std::endl;

    }else if(remaining_interior_vertices_count == 1){
        std::cout<<" SUCCESFULLY COLLAPSED MESH TO A SINGLE INTERIOR POINT AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" iteration count: "<<iteration_count<<std::endl;
        only_one_resulting_interior_vertex = true;

    }else{
        std::cout<<" STOPPED COLLAPSING WITHOUT GETTING A SINGLE INTERIOR VERTEX AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" ===========> REMAINING INTERIOR VERTICES: "<<remaining_interior_vertices_count<<std::endl;
    }

    std::cout<<" ------------------------------ other info ----------"<<std::endl;
    std::cout<<" =====>  FACE-VALENCE != CELL-VALENCE EDGES: "<<findFaceValenceNotEqualCellValenceEdges().size()<<std::endl;
    std::cout<<" =====>            VALENCE-3 INTERIOR EDGES: "<<topo_helper_.valence3InteriorEdgesCount()<<std::endl;
    std::cout<<" =====>         VALENCE-4 INTERIOR VERTICES: "<<topo_helper_.valence4InteriorVerticesCount()<<std::endl;
    std::cout<<" ----------------------------------------------------"<<std::endl;

    std::cout<<" ============================================================== "<<std::endl;

    setLinkConditionPropStatus(false);

    //mesh_.collect_garbage();

    return only_one_resulting_interior_vertex ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;
}



COLLAPSE_RESULT_STATUS EdgeCollapser::collapseUntilPositiveVolume(){
    std::cout<<" --------------------------------"<<std::endl;
    std::cout<<" COLLAPSING ALL COLLAPSIBLE EDGES.... "<<std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    bool collapsed_at_least_one_edge(false);

    auto bad_tets = BadTetFinder::findBadTets(mesh_);
    std::cout<<" -- initial bad tets: "<<std::endl;
    std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
    std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS: "<<std::endl;
    for(auto flipped_tet: bad_tets.second){
        std::cout<<flipped_tet<<" ";
    }
    std::cout<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;

    int previous_deg_tets_count = bad_tets.first.size();
    int previous_flipped_tets_count = bad_tets.second.size();

    auto edge_to_collapse = findNonBoundaryCollapsibleEdgeOnBadTet(bad_tets);

    int collapse_count(0);


    while(edge_to_collapse.idx() != -1){

        collapsed_at_least_one_edge = true;
        //safe_collapseEdge(edge_to_collapse);
        auto collapse = collapseEdge(edge_to_collapse);

        if(collapse.result_status == COLLAPSE_ERROR){
            std::cout<<" ERROR WHILE COLLAPSING ALL EDGES. STOPPING PROCESS"<<std::endl;
            return COLLAPSE_ERROR;
        }

        EdgeSplit empty_split;
        split_and_collapse_sequence_.push_back({empty_split, collapse});
        collapse_sequence_.push_back(collapse);

        collapse_count++;

        bad_tets = BadTetFinder::findBadTets(mesh_);
        edge_to_collapse = findNonBoundaryCollapsibleEdgeOnBadTet(bad_tets);


        int current_deg_tets_count = bad_tets.first.size();
        int current_flipped_tets_count = bad_tets.second.size();



        /*if(!(collapse_count % 100)){
            std::cout<<" ... collapsed "<<collapse_count<<" edges "<<std::endl;
            std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
            std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS"<<std::endl;
            std::cout<<" --------------------------------"<<std::endl;
        }*/

        std::cout<<" ... collapsed "<<collapse_count<<" edges "<<std::endl;
        std::cout<<" ----- "<<current_deg_tets_count<<" DEGENERATE TETS"<<std::endl;
        std::cout<<" ----- "<<current_flipped_tets_count<<" FLIPPED TETS: "<<std::endl;
        for(auto flipped_tet: bad_tets.second){
            std::cout<<flipped_tet<<" ";
        }
        std::cout<<std::endl;
        std::cout<<" --------------------------------"<<std::endl;

        if(current_deg_tets_count > previous_deg_tets_count){
            std::cout<<" --> collapse created more degenerate tets. stopping"<<std::endl;
            //break;
        }

        if(current_flipped_tets_count > previous_flipped_tets_count){
            std::cout<<" --> collapse created more degenerate tets. stopping"<<std::endl;
            //break;
        }

        previous_deg_tets_count = current_deg_tets_count;
        previous_flipped_tets_count = current_flipped_tets_count;

    }

    //bad_cycles_count_.push_back({0, findNonLinkInteriorEdgeCycles().size()});


    total_collapse_count_ += collapse_count;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout<<" COLLAPSED ALL COLLAPSIBLE EDGES AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
    std::cout<<"   ===> REMAINING INTERIOR VERTICES: "<<topo_helper_.interiorVerticesCount()<<std::endl;

    bad_tets = BadTetFinder::findBadTets(mesh_);
    std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
    std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS"<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;


    setLinkConditionPropStatus(false);

    return (bad_tets.first.size() || bad_tets.second.size()) ? COLLAPSE_FAILURE : COLLAPSE_SUCCESS;


}



COLLAPSE_RESULT_STATUS EdgeCollapser::collapseAllCollapsibleEdges(){

    std::cout<<" --------------------------------"<<std::endl;
    std::cout<<" COLLAPSING ALL COLLAPSIBLE EDGES.... "<<std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    bool collapsed_at_least_one_edge(false);

    auto edge_to_collapse = findAnyNonBoundaryCollapsibleEdge();

    int collapse_count(0);

    while(edge_to_collapse.idx() != -1){


        collapsed_at_least_one_edge = true;
        //safe_collapseEdge(edge_to_collapse);
        auto collapse = collapseEdge(edge_to_collapse);

        if(collapse.result_status == COLLAPSE_ERROR){
            std::cout<<" ERROR WHILE COLLAPSING ALL EDGES. STOPPING PROCESS"<<std::endl;
            return COLLAPSE_ERROR;
        }

        EdgeSplit empty_split;
        split_and_collapse_sequence_.push_back({empty_split, collapse});
        collapse_sequence_.push_back(collapse);

        collapse_count++;

        edge_to_collapse = findAnyNonBoundaryCollapsibleEdge();

        if(!(collapse_count % 100)){
            std::cout<<" ... collapsed "<<collapse_count<<" edges "<<std::endl;
        }
    }

    //bad_cycles_count_.push_back({0, findNonLinkInteriorEdgeCycles().size()});


    total_collapse_count_ += collapse_count;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout<<" COLLAPSED ALL COLLAPSIBLE EDGES AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
    std::cout<<"   ===> REMAINING INTERIOR VERTICES: "<<topo_helper_.interiorVerticesCount()<<std::endl;

#if 0
    auto bad_tets = bad_tet_finder_.findBadTets();
    bad_tets = bad_tet_finder_.findBadTets();
    std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
    std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS"<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;
    degenerates_count = bad_tets.first.size();
#endif


    setLinkConditionPropStatus(false);

    return collapsed_at_least_one_edge ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;

}








COLLAPSE_RESULT_STATUS EdgeCollapser::collapseAllCollapsibleBoundaryEdges(){

    std::cout<<" --------------------------------"<<std::endl;
    std::cout<<" COLLAPSING ALL COLLAPSIBLE BOUNDARY EDGES.... "<<std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    bool collapsed_at_least_one_edge(false);

    auto edge_to_collapse = findAnyBoundaryCollapsibleEdge();

    int collapse_count(0);

    while(edge_to_collapse.idx() != -1){


        collapsed_at_least_one_edge = true;
        //safe_collapseEdge(edge_to_collapse);
        auto collapse = collapseEdge(edge_to_collapse);

        if(collapse.result_status == COLLAPSE_ERROR){
            std::cout<<" ERROR WHILE COLLAPSING ALL BOUNDARY EDGES. STOPPING PROCESS"<<std::endl;
            return COLLAPSE_ERROR;
        }

        EdgeSplit empty_split;
        split_and_collapse_sequence_.push_back({empty_split, collapse});
        collapse_sequence_.push_back(collapse);


        collapse_count++;

        edge_to_collapse = findAnyBoundaryCollapsibleEdge();

        if(!(collapse_count % 100)){
            std::cout<<" ... collapsed "<<collapse_count<<" boundary edges "<<std::endl;
        }
    }

    //bad_cycles_count_.push_back({0, findNonLinkInteriorEdgeCycles().size()});


    total_collapse_count_ += collapse_count;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout<<" COLLAPSED ALL COLLAPSIBLE BOUNDARY EDGES AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
    std::cout<<"   ===> REMAINING INTERIOR VERTICES: "<<topo_helper_.interiorVerticesCount()<<std::endl;
    std::cout<<"   ===> REMAINING BOUNDARY VERTICES: "<<topo_helper_.boundaryVerticesCount()<<std::endl;

#if 0
    auto bad_tets = bad_tet_finder_.findBadTets();
    bad_tets = bad_tet_finder_.findBadTets();
    std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
    std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS"<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;
    degenerates_count = bad_tets.first.size();
#endif


    setLinkConditionPropStatus(false);

    return collapsed_at_least_one_edge ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;

}




COLLAPSE_RESULT_STATUS EdgeCollapser::collapseAllCollapsibleEdgesInRandomOrder(){

    std::cout<<" --------------------------------"<<std::endl;
    std::cout<<" COLLAPSING ALL COLLAPSIBLE EDGES IN RANDOM ORDER.... "<<std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    bool collapsed_at_least_one_edge(false);
    bool collapsed_at_least_one_edge_during_iteration(false);

    int collapse_count(0);
    int iteration_count(0);

    do{
        collapsed_at_least_one_edge_during_iteration = false;
        std::vector<EdgeHandle> to_collapse;

        for(auto e: mesh_.edges()){
            if(topo_helper_.isCollapsible(e)){
                to_collapse.push_back(e);
            }
        }
        std::cout<<" -- edges to collapse at iteration "<<iteration_count<<": "<<to_collapse.size()<<std::endl;

        shuffleEdgeVector(to_collapse);

        for(auto e: to_collapse){
            if(!mesh_.is_deleted(e) && topo_helper_.isCollapsible(e)){
                //collapseEdge(mesh_.halfedge_handle(e, 0));
                collapseEdge(mesh_.halfedge_handle(e, 0));
                collapse_count++;
                collapsed_at_least_one_edge_during_iteration = true;
                collapsed_at_least_one_edge = true;

                if(!(collapse_count % 100)){
                    std::cout<<" ...collapsed "<<collapse_count<<" edges"<<std::endl;
                }
            }
        }

        iteration_count++;

    }while(collapsed_at_least_one_edge_during_iteration && iteration_count < 50);

    total_collapse_count_ += collapse_count;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout<<" COLLAPSED ALL COLLAPSIBLE EDGES AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;

    auto bad_tets = BadTetFinder::findBadTets(mesh_);
    std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
    std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS"<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;


    setLinkConditionPropStatus(false);

    return collapsed_at_least_one_edge ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;

}



COLLAPSE_RESULT_STATUS EdgeCollapser::breakTriTets(bool random_order){



    std::cout<<" --------------------------------"<<std::endl;
    std::cout<<" BREAKING ALL TRI-TETS"<<(random_order ? " IN RANDOM ORDER" : "")<< ".... "<<std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<EdgeHandle> to_split;

    bool collapsed_at_least_one_edge(false);

    int collapse_count(0);

    int non_link_to_split_count(0);

    for(auto e: mesh_.edges()){
        if(!mesh_.is_boundary(e) /*&& e.idx() <= original_edges_last_idx_*/){
            int valence(0);
            for(auto c_it = mesh_.ec_iter(e); c_it.valid(); c_it++){
                valence++;
            }

            if(valence == 3){
                to_split.push_back(e);
            }
        }
    }

    std::cout<<" ----> EDGES TO SPLIT: "<<to_split.size()<<std::endl;
    std::cout<<"   among which "<<non_link_to_split_count<<" were non-link"<<std::endl;

    int non_valence_3_split_attempts(0);

    if(random_order){
        shuffleEdgeVector(to_split);
    }

    int failed_count(0);
    int failed_splits_because_non_link_count(0);
    int failed_splits_because_non_interior_count(0);
    int failed_splits_because_non_link_and_non_interior_count(0);
    int semi_boundary_failed_edge_count(0);
    int semi_boundary_succeeded_edge_count(0);

    /*int link_failed_edge_count(0);
    int link_succeeded_edge_count(0);
    int non_link_failed_edge_count(0);
    int non_link_succeeded_edge_count(0);*/



    for(auto e: to_split){


        /*if(collapse_sequence_.size() >= 87){
            std::cout<<" WARNING - HARD-CODED BREAK"<<std::endl;
            //return COLLAPSE_ERROR;
        }*/

        if(mesh_.is_deleted(e)){
            std::cout<<" ABOUT TO SPLIT A DELETED EDGE !"<<std::endl;
            return COLLAPSE_ERROR;
        }

        if(topo_helper_.edgeValence(e) != 3){
            non_valence_3_split_attempts++;
            continue;
        }

        bool semi_boundary = mesh_.is_boundary(mesh_.edge(e).from_vertex()) || mesh_.is_boundary(mesh_.edge(e).to_vertex());
        //bool link_edge = link_condition(mesh_, e);

        COLLAPSE_RESULT_STATUS collapse_result = splitAndCollapseNewNeighbor(e);

        collapsed_at_least_one_edge |= (collapse_result == COLLAPSE_SUCCESS);

        switch(collapse_result){
        case COLLAPSE_ERROR:{
            std::cout<<" ERROR WHILE PERFORMING EDGE COLLAPSE. STOPPING PROCESS"<<std::endl;
            return COLLAPSE_ERROR;
        }
        case COLLAPSE_SUCCESS: {
            //std::cout<<" -- collapsing edge "<<edge_to_collapse<<": "<<mesh_.halfedge(edge_to_collapse)<<std::endl;
            collapse_count++;

            if(semi_boundary){
                semi_boundary_succeeded_edge_count++;
            }

            /*if(link_edge){
                link_succeeded_edge_count++;
            }else{
                non_link_succeeded_edge_count++;
            }*/

            if(!(collapse_count % 100)){
                std::cout<<" ... collapsed "<<collapse_count<<" edges "<<std::endl;
            }
            break;
        }
        case NON_INTERIOR_EDGE_FAILURE:{
            failed_splits_because_non_interior_count++;
            break;
        }
        case NON_LINK_EDGE_FAILURE:{
            failed_splits_because_non_link_count++;
            break;
        }
        case NON_INTERIOR_AND_NON_LINK_EDGE_FAILURE:{
            failed_splits_because_non_link_and_non_interior_count++;
            break;
        }
        default:{}
        }

        if(collapse_result != COLLAPSE_SUCCESS){

            failed_count++;

            if(!semi_boundary){
                semi_boundary_failed_edge_count++;
            }

            /*if(link_edge){
                link_failed_edge_count++;
            }else{
                non_link_failed_edge_count++;
            }*/

        }
    }

    //bad_cycles_count_.push_back({1, bad_cycles_count});

    auto end_time = std::chrono::high_resolution_clock::now();

    total_collapse_count_ += collapse_count;

    std::cout<<" FULLY BROKE TRI-TETS AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
    /*std::cout<<" -----> non valence-3 edge splits: "<<non_valence_3_split_attempts<<std::endl;
    std::cout<<" -----> failed edge splits: "<<failed_count<<std::endl;
    std::cout<<" -----> failed edge splits because non-link: "<<failed_splits_because_non_link_count<<std::endl;
    std::cout<<" -----> failed edge splits because non-interior: "<<failed_splits_because_non_interior_count<<std::endl;
    std::cout<<" -----> failed edge splits because non-link AND non-interior: "<<failed_splits_because_non_link_and_non_interior_count<<std::endl;
    std::cout<<" ------------------------"<<std::endl;
    std::cout<<" -----> semi-boundary edges failed: "<<semi_boundary_failed_edge_count<<std::endl;
    std::cout<<" -----> semi-boundary edges succeded: "<<semi_boundary_succeeded_edge_count<<std::endl;
    std::cout<<" ------------------------"<<std::endl;*/
    /*std::cout<<" -----> failed    non-link edges: "<<non_link_failed_edge_count<<std::endl;
    std::cout<<" -----> failed    link     edges: "<<link_failed_edge_count<<std::endl;
    std::cout<<" -----> succeeded non-link edges: "<<non_link_succeeded_edge_count<<std::endl;
    std::cout<<" -----> succeeded link     edges: "<<link_succeeded_edge_count<<std::endl;
*/

    setLinkConditionPropStatus(false);

    return collapsed_at_least_one_edge ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;
}







COLLAPSE_RESULT_STATUS EdgeCollapser::splitAndCollapseNewNeighbor(OpenVolumeMesh::EdgeHandle eh,
                                                                                 COLLAPSE_SELECTION_HEURISTIC collapse_heuristic){

    setLinkConditionPropStatus(false);

    auto saved_mesh = mesh_;

    //std::set<EdgeHandle> reference_edge_link;

    //std::cout<<" -- splitting edge "<<eh<<": "<<mesh_.edge(eh)<<std::endl;
    //printOutHalfedgeInfo(mesh_.halfedge_handle(eh, 0));

    auto original_from_vertex = mesh_.edge(eh).from_vertex();
    auto original_to_vertex   = mesh_.edge(eh).to_vertex();

    auto split = splitEdge(eh);

    //std::cout<<" -- new vertex: "<<new_vertex<<std::endl;

    std::vector<HalfEdgeHandle> collapsible_new_edges;
    HalfEdgeHandle edge_to_collapse_if_failed(-1);


    COLLAPSE_RESULT_STATUS collapse_result(COLLAPSE_FAILURE);

    int new_halfedge_count(0);
    int partially_boundary_new_halfedge_count(0);
    int non_link_new_halfedge_count(0);
    int back_up_new_halfedge_count(0);


    auto out_he_it = mesh_.outgoing_halfedges(split.middle_vertex).first;
    while(out_he_it.valid()){

        if(collapse_heuristic == RANDOM_SELECTION &&
                collapsible_new_edges.size()){
            break;
        }

        auto to_vertex = mesh_.to_vertex_handle(*out_he_it);

        //ignoring new halfedge touching the boundary
        if(!mesh_.is_boundary(to_vertex)){

            //if it's connected to either original vertices, then it's a back-up halfedge
            if(to_vertex == original_from_vertex ||
                    to_vertex == original_to_vertex){

                if(!link_condition(mesh_, mesh_.edge_handle(*out_he_it))){
                    std::cout<<" ERROR - Newly created edge is not link-condition"<<std::endl;
                    return COLLAPSE_ERROR;
                }else{
                    edge_to_collapse_if_failed = *out_he_it;
                    back_up_new_halfedge_count++;
                    //std::cout<<" ------> found back-up halfedge "<<*out_he_it<<std::endl;
                }
            }else{//otherwise it's a potentially collapsible edge
                if(topo_helper_.isCollapsible(mesh_.edge_handle(*out_he_it))){
                    //std::cout<<" ------> found new collapsible halfedge "<<*out_he_it<<std::endl;
                    collapsible_new_edges.push_back(*out_he_it);
                }else{
                    non_link_new_halfedge_count++;
                    //std::cout<<" ------> new edge is not collapsible "<<std::endl;
                }
            }
        }else{
            //std::cout<<" -----> uncollapsible because attached to the boundary"<<std::endl;
            partially_boundary_new_halfedge_count++;
        }

        new_halfedge_count++;
        out_he_it++;
    }

    HalfEdgeHandle edge_to_collapse(-1);

    //edge to collapse selection
    if(collapsible_new_edges.size()){
        //std::cout<<"--------------------------------------------"<<std::endl;
        //std::cout<<" collapsible new edges: "<<collapsible_new_edges.size()<<std::endl;
        switch(collapse_heuristic){
        case RANDOM_SELECTION:{
            //std::cout<<"  -> using random heuristic"<<std::endl;

            edge_to_collapse = collapsible_new_edges[0];
            break;
        }
#warning TODO: test heuristic on to-vertex valence
        case MIN_VALENCE_SELECTION:{
            std::cout<<"  -> using min valence heuristic"<<std::endl;

            auto min_valence(std::numeric_limits<size_t>::max());

            for(auto he: collapsible_new_edges){
                auto to_vertex = mesh_.to_vertex_handle(he);
                std::cout<<" vertex "<<to_vertex<<" valence "<<mesh_.valence(to_vertex)<<std::endl;
                if(mesh_.valence(to_vertex) < min_valence){
                    std::cout<<"    --> updated edge to collapse to "<<he<<": "<<mesh_.halfedge(he)<<std::endl;
                    min_valence = mesh_.valence(to_vertex);
                    edge_to_collapse = he;
                }
            }
            break;

        }
        default:{
            //TODO
        }

        }
    }


    if(edge_to_collapse.idx() != -1){
        if(mesh_.is_deleted(edge_to_collapse)){
            std::cout<<" ABOUT TO COLLAPSE A DELETED EDGE !"<<std::endl;
            return COLLAPSE_ERROR;
        }

        auto collapse = collapseEdge(edge_to_collapse);

        collapse.TEMP_after_split = true;

        if(collapse.result_status == COLLAPSE_SUCCESS){
            collapse_result = COLLAPSE_SUCCESS;
            split_and_collapse_sequence_.push_back({split, collapse});
            collapse_sequence_.push_back(collapse);

        }else{
            std::cout<<" ERROR - COLLAPSE FAILED"<<std::endl;
            return COLLAPSE_ERROR;
        }
    }else if(edge_to_collapse_if_failed.idx() != -1){

        if(mesh_.is_deleted(edge_to_collapse_if_failed)){
            std::cout<<" ABOUT TO COLLAPSE A BACKUP DELETED EDGE !"<<std::endl;
            return COLLAPSE_ERROR;
        }

        auto collapse = collapseEdge(edge_to_collapse_if_failed);
        split_and_collapse_sequence_.push_back({split, collapse});
        collapse_sequence_.push_back(collapse);


        //std::cout<<" --------------- > CANCELLED SAC "<<std::endl;

        //#warning  brute SAC cancel
        //mesh_ = saved_mesh;


        if(partially_boundary_new_halfedge_count && non_link_new_halfedge_count){
            collapse_result = NON_INTERIOR_AND_NON_LINK_EDGE_FAILURE;
        }else if(partially_boundary_new_halfedge_count){
            collapse_result = NON_INTERIOR_EDGE_FAILURE;
        }else if(non_link_new_halfedge_count){
            collapse_result = NON_LINK_EDGE_FAILURE;
        }else{
            std::cout<<" ERROR - FAILED TO COLLAPSE BUT NEITHER BECAUSE OF NON-LINK OR NON-INTERIOR"<<std::endl;
            std::cout<<" partially boundary new halfedges: "<<partially_boundary_new_halfedge_count<<std::endl;
            std::cout<<"           non link new halfedges: "<<non_link_new_halfedge_count<<std::endl;
            std::cout<<" -----------------------------------------------"<<std::endl;

            return COLLAPSE_ERROR;
        }

        //std::cout<<" ----- canceled collapse"<<std::endl;

    }else{
        std::cout<<" COULDN'T FIND A BACK-UP EDGE, WHICH IS WEIRD"<<std::endl;
        return COLLAPSE_ERROR;
    }


    /*auto bad_tets = bad_tet_finder.findBadTets();
    if(bad_tets.first.size()){
        std::cout<<" split-and-collapse created a bad tet: "<<std::endl;
        std::cout<<" ------: "<<bad_tets.first.size()<<" degenerates and "<<std::endl;
        std::cout<<" ------: "<<bad_tets.second.size()<<" flipped"<<std::endl;
        std::cout<<" --------------------------------------------------------"<<std::endl;
        return COLLAPSE_ERROR;
    }*/


    return collapse_result;

}


COLLAPSE_RESULT_STATUS EdgeCollapser::breakCycle(const EdgeCycle& cycle){

    auto links_intersection = link(mesh_, cycle[0]).intersection(link(mesh_, cycle[1])).intersection(link(mesh_, cycle[2]));

    auto potential_fan_tips = links_intersection.vertices();

    COLLAPSE_RESULT_STATUS result(COLLAPSE_FAILURE);

    for(auto fan_tip_vertex: potential_fan_tips){
        if(result == COLLAPSE_SUCCESS){
            break;
        }

        //now try to collapse one of the edges connecting the fan tip
        auto cycle_vertices = cycle.toVertexSet(mesh_);

        for(auto out_he_it = mesh_.outgoing_halfedges(fan_tip_vertex).first; out_he_it.valid(); out_he_it++){
            if(result == COLLAPSE_SUCCESS){
                break;
            }
            if(cycle_vertices.find(mesh_.to_vertex_handle(*out_he_it)) != cycle_vertices.end() &&
                    topo_helper_.isCollapsible(mesh_.edge_handle(*out_he_it))){

                collapseEdge(*out_he_it);
                result = COLLAPSE_SUCCESS;
            }
        }
    }


    return potential_fan_tips.size() ? result : COLLAPSE_ERROR;
}



void EdgeCollapser::splitAllBoundaryEdges(){

    std::vector<HalfEdgeHandle> halfedges;
    for(auto e: mesh_.edges()){
        if(mesh_.is_boundary(e)){
            halfedges.push_back(mesh_.halfedge_handle(e, 0));
        }
    }

    for(auto he: halfedges){
        mesh_.split_edge(he);
    }
}











COLLAPSE_RESULT_STATUS EdgeCollapser::splitAndCollapseCycles(COLLAPSE_SELECTION_HEURISTIC collapse_heuristic){

    std::cout<<" -------------------------------------------------"<<std::endl;
    std::cout<<" BREAKING ALL CYCLES WITH THE 'SPLIT-AND-COLLAPSE' METHOD USING HEURISTIC "<<collapse_heuristic<<".... "<<std::endl;


    auto start_time = std::chrono::high_resolution_clock::now();

    bool collapsed_at_least_one_edge(false);

    std::vector<int> histogram = std::vector<int>(MAX_EDGE_VALENCE_FOR_HISTOGRAM, 0);
    split_edge_valence_histogram.push_back(histogram);

    //each pair contains <priority_value, edge_handle>
    //NOTE: max-priority
    std::priority_queue<std::pair<int, EdgeHandle>> to_split;


    for(auto e: mesh_.edges()){

        if(!mesh_.is_boundary(e) &&
                !link_condition(mesh_, e) /*&&
                        e.idx() <= original_edges_last_idx_*/){

            int priority(0);
            switch(collapse_heuristic){
            case MIN_VALENCE_SELECTION:{
                //using negative valence so lower valence are treated first
                priority = -mesh_.valence(e);
                break;
            }
            case MAX_VALENCE_SELECTION:{
                priority = mesh_.valence(e);
                break;
            }
            default:{}//leave priority constant -> no particular heuristic
            }

            to_split.push({priority, e});
        }
    }
    std::cout<<" --> EDGES TO SPLIT: "<<to_split.size()<<" edges"<<std::endl;



    int attempted_count(0);
    int succeeded_count(0);
    int skipped_count(0);
    int failed_count(0);
    int semi_boundary_succeeded_edge_count(0);
    int semi_boundary_failed_edge_count(0);
    int failed_splits_because_non_interior_count(0);
    int failed_splits_because_non_link_count(0);
    int failed_splits_because_non_link_and_non_interior_count(0);

    while(!to_split.empty()){
        attempted_count++;

        auto top = to_split.top();
        auto e = top.second;
        to_split.pop();

        //std::cout<<" trying to split-and-collapse edge "<<mesh_.edge(e)<<", valence "<<top.first<<std::endl;

        if(mesh_.is_deleted(e)){
            std::cout<<" ABOUT TO SPLIT A DELETED EDGE !"<<std::endl;
            return COLLAPSE_ERROR;
        }

        if(link_condition(mesh_, e)){
            skipped_count++;
            //std::cout<<"  --> skipped"<<std::endl;
            continue;
        }

        //analysis stuff
        bool semi_boundary = mesh_.is_boundary(mesh_.edge(e).from_vertex()) || mesh_.is_boundary(mesh_.edge(e).to_vertex());

        // histogram stuff
        auto valence = mesh_.valence(e);
        auto histogram_index = std::min(valence,
                                        size_t(MAX_EDGE_VALENCE_FOR_HISTOGRAM-1));

        //actual collapse
        COLLAPSE_RESULT_STATUS collapse_result = splitAndCollapseNewNeighbor(e);

        collapsed_at_least_one_edge |= (collapse_result == COLLAPSE_SUCCESS);

        switch(collapse_result){
        case COLLAPSE_ERROR:{
            std::cout<<" ERROR WHILE SPLIT-AND-COLLAPSING CYCLES. STOPPING PROCESS"<<std::endl;
            return COLLAPSE_ERROR;
        }
        case COLLAPSE_SUCCESS: {
            //std::cout<<"  --> succeeded"<<std::endl;
            //std::cout<<" -- collapsing edge "<<edge_to_collapse<<": "<<mesh_.halfedge(edge_to_collapse)<<std::endl;
            succeeded_count++;

            split_edge_valence_histogram.back()[histogram_index]++;

            if(semi_boundary){
                semi_boundary_succeeded_edge_count++;
            }

            if(!(succeeded_count % 100)){
                std::cout<<" ... collapsed "<<succeeded_count<<" edges "<<std::endl;
            }
            break;
        }
        case NON_INTERIOR_EDGE_FAILURE:{
            failed_splits_because_non_interior_count++;
            break;
        }
        case NON_LINK_EDGE_FAILURE:{
            failed_splits_because_non_link_count++;
            break;
        }
        case NON_INTERIOR_AND_NON_LINK_EDGE_FAILURE:{
            failed_splits_because_non_link_and_non_interior_count++;
            break;
        }
        default:{}
        }

        if(collapse_result != COLLAPSE_SUCCESS){

            //std::cout<<" --> failed"<<std::endl;
            failed_count++;

            if(semi_boundary){
                semi_boundary_failed_edge_count++;
            }
        }
    }


    auto end_time = std::chrono::high_resolution_clock::now();

    total_collapse_count_ += succeeded_count;

    std::cout<<" FULLY BROKE BAD CYCLES AFTER "<<attempted_count<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
    std::cout<<" -----> succeeded: "<<succeeded_count<<std::endl;
    std::cout<<" -----> skipped  : "<<skipped_count<<std::endl;
    std::cout<<" -----> failed   : "<<failed_count<<std::endl;
    std::cout<<" -----> failed edge splits because non-link: "<<failed_splits_because_non_link_count<<std::endl;
    std::cout<<" -----> failed edge splits because non-interior: "<<failed_splits_because_non_interior_count<<std::endl;
    std::cout<<" -----> failed edge splits because non-link AND non-interior: "<<failed_splits_because_non_link_and_non_interior_count<<std::endl;
    std::cout<<" ------------------------"<<std::endl;
    std::cout<<" -----> semi-boundary edges failed: "<<semi_boundary_failed_edge_count<<std::endl;
    std::cout<<" -----> semi-boundary edges succeded: "<<semi_boundary_succeeded_edge_count<<std::endl;


    setLinkConditionPropStatus(false);

    std::cout<<" last histogram: "<<std::endl;
    for(size_t i(0); i<split_edge_valence_histogram[0].size(); i++){
        std::cout<<i<<": ";
        for(size_t j(0); j<split_edge_valence_histogram.size(); j++){
            std::cout<<std::setw(4)<<split_edge_valence_histogram[j][i]<<" ";
        }
        std::cout<<std::endl;
    }

    std::cout<<" -------------------------------------------------"<<std::endl;

    return succeeded_count ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;

}




COLLAPSE_RESULT_STATUS EdgeCollapser::threeStepCollapse(){

    std::cout<<" ============================================================== "<<std::endl;

    std::cout<<" RUNNING THREE-STEP EDGE COLLAPSE.... "<<std::endl;



    COLLAPSE_RESULT_STATUS overall_result(COLLAPSE_FAILURE);

    int remaining_interior_vertices_count(0);
    COLLAPSE_RESULT_STATUS cycle_breaking_result(COLLAPSE_SUCCESS);
    bool collapsed_at_least_one_edge(true);
    bool collapsed_at_least_one_edge_this_iteration(true);

    int iteration_count(0);

    auto start_time = std::chrono::high_resolution_clock::now();


    while(overall_result != COLLAPSE_ERROR &&
          remaining_interior_vertices_count != 1 &&
          collapsed_at_least_one_edge_this_iteration &&
          iteration_count < 50){

        std::cout<<" ==================================== ITERATION "<<iteration_count<<std::endl;

        collapsed_at_least_one_edge_this_iteration = false;


        //---> two-step collapse
        remaining_interior_vertices_count = twoStepCollapse();
        if(remaining_interior_vertices_count == -1){
            std::cout<<" ERROR - STOPPING THREE-STEP COLLAPSE"<<std::endl;
            overall_result = COLLAPSE_ERROR;
            break;
        }

        if(remaining_interior_vertices_count == 1){
            break;
        }

        //cycles split-and-collapse
        cycle_breaking_result = splitAndCollapseCycles(MIN_VALENCE_SELECTION);
        if(cycle_breaking_result == COLLAPSE_SUCCESS){
            collapsed_at_least_one_edge_this_iteration = true;
            collapsed_at_least_one_edge = true;
            std::cout<<" --> CYCLES BREAKING SUCCESFUL"<<std::endl;
        }else if(cycle_breaking_result == COLLAPSE_ERROR){
            std::cout<<" ERROR - STOPPING THREE-STEP COLLAPSE"<<std::endl;
            overall_result = COLLAPSE_ERROR;
            break;
        }

        iteration_count++;
    }


    auto end_time = std::chrono::high_resolution_clock::now();

    bool only_one_resulting_interior_vertex(false);

    if(overall_result == COLLAPSE_ERROR){
        std::cout<<" ERROR - SOMETHING BAD HAPPENED DURING COLLAPSE AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" =====> REMAINING INTERIOR VERTICES: "<<topo_helper_.interiorVerticesCount()<<std::endl;

    }else if(topo_helper_.interiorVerticesCount() == 1){
        std::cout<<" SUCCESFULLY COLLAPSED MESH TO A SINGLE INTERIOR POINT AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" iteration count: "<<iteration_count<<std::endl;
        only_one_resulting_interior_vertex = true;

    }else{
        std::cout<<" STOPPED COLLAPSING WITHOUT GETTING A SINGLE INTERIOR VERTEX AFTER "<<total_collapse_count_<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
        std::cout<<" =====> REMAINING INTERIOR VERTICES: "<<topo_helper_.interiorVerticesCount()<<std::endl;
    }

    std::cout<<" ============================================================== "<<std::endl;


    setLinkConditionPropStatus(false);


    //mesh_.collect_garbage();

    return only_one_resulting_interior_vertex ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;
}





std::vector<EdgeCollapse> EdgeCollapser::getCollapseSequence() const{
    return collapse_sequence_;
}

std::vector<SplitAndCollapse> EdgeCollapser::getSplitAndCollapseSequence() const{
    return split_and_collapse_sequence_;
}


const SubTriangleMap& EdgeCollapser::getSubTriangleMap() const{
    return sub_triangle_map_;
}




COLLAPSE_RESULT_STATUS EdgeCollapser::breakFanTets(){


    std::cout<<" --------------------------------"<<std::endl;
    std::cout<<" BREAKING ALL FAN-TETS.... "<<std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<EdgeHandle> to_split;

    bool collapsed_at_least_one_edge(false);


    std::vector<EdgeHandle> to_collapse;

    auto candidate_cycles = findNonLinkInteriorEdgeCycles();

    int breakable_fan_tets_count(0);
    int unbreakable_fan_tets_count(0);
    int tri_tets_count(0);
    int other_configs_count(0);
    int unbreakable_because_non_interior(0);
    int unbreakable_because_non_link(0);

    for(auto cycle: candidate_cycles){
        std::cout<<" -- trying to break cycle: ";
        cycle.print();

        auto links_intersection = link(mesh_, cycle[0]).intersection(link(mesh_, cycle[1])).intersection(link(mesh_, cycle[2]));

        auto potential_fan_tips = links_intersection.vertices();

        if(potential_fan_tips.size()){
            bool found_edge_to_collapse(false);

            if(links_intersection.edges().size()){
                tri_tets_count++;
            }


            for(auto fan_tip_vertex: potential_fan_tips){
                if(found_edge_to_collapse){
                    break;
                }

                //now try to collapse one of the edges connecting the fan tip

                std::set<VertexHandle> cycle_vertices;
                for(auto e: cycle){
                    cycle_vertices.insert(mesh_.edge(e).from_vertex());
                    cycle_vertices.insert(mesh_.edge(e).to_vertex());
                }


                for(auto v: cycle_vertices){
                    if(found_edge_to_collapse){
                        break;
                    }
                    for(auto out_he_it = mesh_.outgoing_halfedges(v).first; out_he_it.valid(); out_he_it++){
                        if(found_edge_to_collapse){
                            break;
                        }
                        if(mesh_.to_vertex_handle(*out_he_it) == fan_tip_vertex){
                            //std::cout<<" ------ checking edge "<<*out_he_it<<": "<<mesh_.halfedge(*out_he_it)<<std::endl;
                            if(!link_condition(mesh_,mesh_.edge_handle(*out_he_it))){
                                unbreakable_because_non_link++;
                            }

                            if(!topo_helper_.isInterior(mesh_.edge_handle(*out_he_it))){
                                unbreakable_because_non_interior++;
                            }

                            if(topo_helper_.isCollapsible(mesh_.edge_handle(*out_he_it))){
                                //std::cout<<" ----------> collapsible"<<std::endl;
                                to_collapse.push_back(mesh_.edge_handle(*out_he_it));
                                found_edge_to_collapse = true;
                            }
                        }
                    }
                }
            }

            if(found_edge_to_collapse){
                breakable_fan_tets_count++;
            }else{
                unbreakable_fan_tets_count++;
            }
        }else{
            other_configs_count++;
        }
    }


    std::cout<<" ---------------- "<<std::endl;
    std::cout<<" OVER "<<candidate_cycles.size()<<" BAD CYCLES, WE HAVE: "<<std::endl;
    std::cout<<"   BREAKABLE FAN-TETS : "<<breakable_fan_tets_count<<std::endl;
    std::cout<<" UNBREAKABLE FAN-TETS : "<<unbreakable_fan_tets_count<<std::endl;
    std::cout<<" OTHER CONFIGURATIONS : "<<other_configs_count<<std::endl;
    std::cout<<"              ------"<<std::endl;
    std::cout<<" UNCOLLAPSIBLE EDGES BECAUSE NON-LINK    : "<<unbreakable_because_non_link<<std::endl;
    std::cout<<" UNCOLLAPSIBLE EDGES BECAUSE NON-INTERIOR: "<<unbreakable_because_non_interior<<std::endl;
    std::cout<<"              ------"<<std::endl;
    std::cout<<"           (TRI-TETS) : "<<tri_tets_count<<std::endl;
    std::cout<<" ---------------- "<<std::endl;



    int collapse_count(0);
    int now_deleted_edges_count(0);
    int no_longer_collapsible_edges_count(0);

    for(auto e: to_collapse){
        if(!mesh_.is_deleted(e) && topo_helper_.isCollapsible(e)){
            collapseEdge(mesh_.halfedge_handle(e, 0));
            collapse_count++;
        }else{
            if(mesh_.is_deleted(e)){
                now_deleted_edges_count++;
            }
            if(!topo_helper_.isCollapsible(e)){
                no_longer_collapsible_edges_count++;
            }
        }
    }


    auto end_time = std::chrono::high_resolution_clock::now();

    total_collapse_count_ += collapse_count;

    std::cout<<" FULLY BROKE FAN-TETS AFTER "<<collapse_count<<" COLLAPSES IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
    std::cout<<" ("<<now_deleted_edges_count<<") edges were deleted in the meantime and ("<<no_longer_collapsible_edges_count<<") were no longer collapsible"<<std::endl;

    auto bad_tets = BadTetFinder::findBadTets(mesh_);
    std::cout<<" ----- "<<bad_tets.first.size()<<" DEGENERATE TETS"<<std::endl;
    std::cout<<" ----- "<<bad_tets.second.size()<<" FLIPPED TETS"<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;




    int bad_cycles_count = findNonLinkInteriorEdgeCycles().size();
    bad_cycles_count_.push_back({2,bad_cycles_count});


    std::cout<<" --> bad cycles after breaking fan-tets: "<<bad_cycles_count<<std::endl;
    std::cout<<" --------------------------------"<<std::endl;

    setLinkConditionPropStatus(false);

    return collapsed_at_least_one_edge ? COLLAPSE_SUCCESS : COLLAPSE_FAILURE;

}









EdgeSplit EdgeCollapser::splitEdge(const OpenVolumeMesh::EdgeHandle& e){
    return splitEdge(mesh_.halfedge_handle(e, 0));
}



EdgeSplit EdgeCollapser::splitEdge(const OpenVolumeMesh::HalfEdgeHandle& he){

    auto from_vertex = mesh_.from_vertex_handle(he);
    auto to_vertex = mesh_.to_vertex_handle(he);

    //std::cout<<" ------ PERFORMING EDGE SPLIT ("<<from_vertex<<"-"<<to_vertex<<")"<<std::endl;

    std::vector<VertexHandle> equatorial_vertices;
    for(auto hehf_it = mesh_.hehf_iter(he); hehf_it.valid(); hehf_it++){
        for(auto hfv_it = mesh_.hfv_iter(*hehf_it); hfv_it.valid(); hfv_it++){
            if(*hfv_it != from_vertex &&
                    *hfv_it != to_vertex){
                equatorial_vertices.push_back(*hfv_it);
            }
        }
    }

    auto middle_vertex = mesh_.split_edge(he);
    //std::cout<<" middle vertex: "<<middle_vertex<<std::endl;

    EdgeSplit result = {from_vertex,
                        to_vertex,
                        middle_vertex,
                        equatorial_vertices};


    register_split_in_subtriangle_map(result);

    //std::cout<<" ------ DONE"<<std::endl;
    return result;
}



void EdgeCollapser::register_split_in_subtriangle_map(const EdgeSplit& split){

    /*std::cout<<" -------- "<<std::endl;
    std::cout<<" register split ("<<from_vertex<<"-"<<to_vertex<<") -> "<<middle_vertex<<" in subtriangle map"<<std::endl;
    std::cout<<" equatorial vertices: ";
    for(auto ev: equatorial_vertices){
        std::cout<<" "<<ev;
    }
    std::cout<<std::endl;*/

    for(auto ev: split.equatorial_vertices){
        SubTriangleFace key({split.from_vertex, split.to_vertex, ev});
        SubTriangleFace value_a({split.from_vertex, split.middle_vertex, ev});
        SubTriangleFace value_b({split.middle_vertex, split.to_vertex, ev});

        sub_triangle_map_.insert(key, value_a);
        sub_triangle_map_.insert(key, value_b);


        /*std::cout<<" subfaces for    key "<<key<<": ";
        for(auto sub_tri: sub_triangle_map_[key]){
            std::cout<<" "<<sub_tri;
        }
        std::cout<<std::endl;

        std::cout<<" subfaces for op key "<<key.opposite_face()<<": ";
        for(auto sub_tri: sub_triangle_map_[key.opposite_face()]){
            std::cout<<" "<<sub_tri;
        }
        std::cout<<std::endl;
        std::cout<<" ---"<<std::endl;*/

    }


    /*std::cout<<" ... done"<<std::endl;
    std::cout<<" -------- "<<std::endl;*/

}


EdgeCollapse EdgeCollapser::collapseEdge(const HalfEdgeHandle& _heh){

    auto from_vertex = mesh_.from_vertex_handle(_heh);
    auto to_vertex = mesh_.to_vertex_handle(_heh);
    auto replacement_vertex = mesh_.add_vertex(mesh_.vertex(to_vertex));
    mesh_.delete_vertex(replacement_vertex);


    //std::cout<<" ------ PERFORMING COLLAPSE ("<<from_vertex<<"-"<<to_vertex<<") -> "<<replacement_vertex<<std::endl;
    /*if(mesh_.valence(from_vertex) == 4){
        std::cout<<" WARNING - collapse from vertex of valence 4"<<std::endl;
    }*/


    /*if(from_vertex.idx() == 933 ){
        std::cout<<" FOUND BAD VERTEX AT COLLAPSE N° "<<collapse_sequence_.size()<<std::endl;

        exit(EXIT_FAILURE);
    }*/

    std::vector<VertexHandle> equatorial_vertices;
    for(auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); hehf_it++){
        for(auto hfv_it = mesh_.hfv_iter(*hehf_it); hfv_it.valid(); hfv_it++){
            if(*hfv_it != from_vertex &&
                    *hfv_it != to_vertex){
                equatorial_vertices.push_back(*hfv_it);
            }
        }
    }

    EquatorialDisc equatorial_disc;
    for(size_t i(0); i<equatorial_vertices.size(); i++){
        int j = (i+1) % equatorial_vertices.size();
        equatorial_disc.insert({to_vertex,
                                equatorial_vertices[i],
                                equatorial_vertices[j]});
    }


    auto result_vertex = mesh_.collapse_edge(_heh);


    EdgeCollapse result =  { from_vertex, to_vertex,
                             replacement_vertex,
                             equatorial_vertices,
                             {},
                             (result_vertex == to_vertex) ? COLLAPSE_SUCCESS : COLLAPSE_ERROR};


    //recover the original hemispherical halffaces
    result.north_hemispherical_hfs = PEHelpers::find_hemispherical_faces(mesh_,
                                                                         to_vertex,
                                                                         equatorial_disc,
                                                                         false);

    //replace the to_vertex with a new one to differentiate between before and after the collapse
    mesh_.swap_vertex_indices(to_vertex,
                              replacement_vertex);

    return result;
}




bool EdgeCollapser::checkThatCyclesAreIndeedNotFaces(){

    auto cycles = findNonLinkInteriorEdgeCycles();

    bool no_cycle_is_a_face(true);

    for(auto cycle: cycles){

        std::set<VertexHandle> cycle_vertices;
        for(auto e: cycle){
            cycle_vertices.insert(mesh_.edge(e).from_vertex());
            cycle_vertices.insert(mesh_.edge(e).to_vertex());
        }

        for(auto f: mesh_.faces()){

            int found_vertices(0);
            for(auto fv_it = mesh_.fv_iter(f); fv_it.valid(); fv_it++){
                if(cycle_vertices.find(*fv_it) != cycle_vertices.end()){
                    found_vertices++;
                }
            }

            if(found_vertices == 3){
                no_cycle_is_a_face = false;
                std::cout<<" cycle ";cycle.print(mesh_); std::cout<<" is actually a face..."<<std::endl;
                auto face_handle = mesh_.halfface(std::vector<VertexHandle>(cycle_vertices.begin(), cycle_vertices.end()));
                if(face_handle.idx() == -1){
                    std::cout<<"  but not according to OVM "<<std::endl;
                }else{
                    std::cout<<"  confirmed by OVM"<<std::endl;
                }
            }
        }
    }

    if(no_cycle_is_a_face){
        std::cout<<" ======> no cycle is a face of the mesh, all good "<<std::endl;
    }

    return no_cycle_is_a_face;
}







int EdgeCollapser::twoStepCollapseInRandomOrder(int nb_runs){

    std::cout<<" ============================================================== "<<std::endl;

    std::cout<<" RUNNING EXPERIMENTAL RANDOM ORDER TWO-STEP COLLAPSE FOR "<<nb_runs<<" RUNS.... "<<std::endl;

    TetrahedralMesh initial_mesh = mesh_;


    int iteration_count(0);
    bool single_interior_vertex(false);
    int min_interior_vertices(mesh_.n_vertices());

    int success_count(0);

    while(iteration_count < nb_runs /*&& !single_interior_vertex*/){
        std::cout<<std::endl;
        std::cout<<"---------------------------------------------"<<std::endl;
        std::cout<<"---------------------------------------------"<<std::endl;
        std::cout<<".... TRYING RANDOM ORDER "<<iteration_count<<std::endl;
        std::cout<<"---------------------------------------------"<<std::endl;
        std::cout<<"---------------------------------------------"<<std::endl;
        std::cout<<std::endl;

        iteration_count++;
        single_interior_vertex = twoStepCollapse(true);
        min_interior_vertices = std::min(min_interior_vertices, topo_helper_.interiorVerticesCount());

        if(single_interior_vertex){
            success_count++;
        }
        mesh_ = initial_mesh;

        /*if(!single_interior_vertex){
            mesh_ = initial_mesh;
        }*/

        /*if(single_interior_vertex){
            std::cout<<"  =======> SUCCESS AFTER "<<iteration_count<<" RANDOM ORDER TRIALS"<<std::endl;
        }else{
            std::cout<<"  =======> MINIMAL NUMBER OF REMAINING INTERIOR VERTICES: "<<min_interior_vertices<<std::endl;
        }*/
    }

    std::cout<<" SUCCESS COUNT: "<<success_count<<"/"<<nb_runs<<std::endl;



    std::cout<<" ============================================================== "<<std::endl;

    mesh_.collect_garbage();

    return success_count;
}

std::set<OpenVolumeMesh::EdgeHandle> EdgeCollapser::findLoneNonLinkEdges(){

    storeLinkConditionAsProperty(false);

    std::set<EdgeHandle> lone_non_link_edges;

    for(auto e: mesh_.edges()){
        if(topo_helper_.isInterior(e) && !link_condition_prop_[e]){
            auto from_vertex = mesh_.edge(e).from_vertex();
            auto to_vertex = mesh_.edge(e).to_vertex();

            auto visited = mesh_.request_vertex_property<bool>("visited");

            bool found_at_least_one_non_link_outgoing_edge(false);
            //mark all vertices neighbor to from_vertex as visited if the edge connecting them is non-link
            for(auto out_he: mesh_.outgoing_halfedges(from_vertex)){
                if(!link_condition_prop_[mesh_.edge_handle(out_he)]){
                    visited[mesh_.to_vertex_handle(out_he)] = true;
                    found_at_least_one_non_link_outgoing_edge = true;
                }
            }

            if(!found_at_least_one_non_link_outgoing_edge){
                lone_non_link_edges.insert(e);
                continue;
            }

            bool found_at_least_one_cycle(false);
            for(auto out_he: mesh_.outgoing_halfedges(to_vertex)){
                if(!link_condition_prop_[mesh_.edge_handle(out_he)] &&
                        visited[mesh_.to_vertex_handle(out_he)]){
                    found_at_least_one_cycle = true;
                }
            }

            if(!found_at_least_one_cycle){
                lone_non_link_edges.insert(e);
            }
        }
    }

    return lone_non_link_edges;
}



std::set<OpenVolumeMesh::EdgeHandle> EdgeCollapser::findFaceValenceNotEqualCellValenceEdges(){
    std::set<EdgeHandle> face_valence_not_equal_cell_valence_edges;

    for(auto e: mesh_.edges()){
        if(topo_helper_.isInterior(e)){
            if((topo_helper_.edgeValence(e) == 3) != (mesh_.valence(e) == 3)){
                face_valence_not_equal_cell_valence_edges.insert(e);
            }
        }
    }

    return face_valence_not_equal_cell_valence_edges;
}


std::set<EdgeCycle> EdgeCollapser::findNonLinkInteriorEdgeCycles(bool include_face_cycles){

    storeLinkConditionAsProperty(false);

    /*std::cout<<" -------------------------------------------- "<<std::endl;
    std::cout<<" LOOKING NON-LINK EDGE CYCLES.... "<<std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();*/

    std::set<EdgeCycle> cycles;

    for(auto e: mesh_.edges()){
        if(topo_helper_.isInterior(e) && !link_condition_prop_[e]){
            //std::cout<<" -- checking non-link edge "<<e<<": "<<mesh_.edge(e)<<std::endl;

            auto from_vertex = mesh_.edge(e).from_vertex();
            auto to_vertex   = mesh_.edge(e).to_vertex();

            std::vector<EdgeHandle> neighbor_non_link_edges;

            std::vector<std::vector<EdgeHandle>> non_link_edge_cycles;

            for(auto out_he_it = mesh_.outgoing_halfedges(from_vertex).first; out_he_it.valid(); out_he_it++){
                auto neighbor_edge = mesh_.edge_handle(*out_he_it);

                if(topo_helper_.isInterior(neighbor_edge) && !link_condition_prop_[neighbor_edge]){
                    auto tip_vertex = mesh_.to_vertex_handle(*out_he_it);

                    for(auto out_out_he_it = mesh_.outgoing_halfedges(tip_vertex).first; out_out_he_it.valid(); out_out_he_it++){
                        auto neighbor_neighbor_edge = mesh_.edge_handle(*out_out_he_it);

                        if(topo_helper_.isInterior(neighbor_neighbor_edge) &&
                                !link_condition_prop_[neighbor_neighbor_edge] &&
                                mesh_.to_vertex_handle(*out_out_he_it) == to_vertex){

                            EdgeCycle cycle({e, neighbor_edge, neighbor_neighbor_edge});
                            auto cycle_vertices = cycle.toVertexSet(mesh_);

                            if(include_face_cycles ||
                                    mesh_.halfface(std::vector<VertexHandle>(cycle_vertices.begin(), cycle_vertices.end())).idx() == -1){
                                cycles.insert(cycle);
                            }
                        }
                    }
                }
            }
        }
    }


    /*auto end_time = std::chrono::high_resolution_clock::now();
    std::cout<<" ANALYSIS PERFORMED IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS. RESULTS: "<<std::endl;
    std::cout<<" ---- non-link edge cycles count: "<<cycles.size()<<std::endl;
    std::cout<<" -------------------------------------------- "<<std::endl;*/

    return cycles;

}






bool EdgeCollapser::containsNonLinkInteriorCycle(){

    storeLinkConditionAsProperty(false);

    for(auto e: mesh_.edges()){
        if(topo_helper_.isInterior(e) && !link_condition_prop_[e]){
            //std::cout<<" -- checking non-link edge "<<e<<": "<<mesh_.edge(e)<<std::endl;

            auto from_vertex = mesh_.edge(e).from_vertex();
            auto to_vertex   = mesh_.edge(e).to_vertex();

            std::vector<EdgeHandle> neighbor_non_link_edges;

            std::vector<std::vector<EdgeHandle>> non_link_edge_cycles;

            for(auto out_he_it = mesh_.outgoing_halfedges(from_vertex).first; out_he_it.valid(); out_he_it++){
                auto neighbor_edge = mesh_.edge_handle(*out_he_it);

                if(topo_helper_.isInterior(neighbor_edge) && !link_condition_prop_[neighbor_edge]){
                    auto tip_vertex = mesh_.to_vertex_handle(*out_he_it);

                    for(auto out_out_he_it = mesh_.outgoing_halfedges(tip_vertex).first; out_out_he_it.valid(); out_out_he_it++){
                        auto neighbor_neighbor_edge = mesh_.edge_handle(*out_out_he_it);

                        if(topo_helper_.isInterior(neighbor_neighbor_edge) &&
                                !link_condition_prop_[neighbor_neighbor_edge] &&
                                mesh_.to_vertex_handle(*out_out_he_it) == to_vertex){

                            bool actually_a_face(false);
                            //check that this is not actually a face
                            for(auto ef_it = mesh_.ef_iter(e); ef_it.valid(); ef_it++){
                                int found_edges(0);

                                for(auto fe_it = mesh_.fe_iter(*ef_it); fe_it.valid(); fe_it++){
                                    if(*fe_it == neighbor_edge ||
                                            *fe_it == neighbor_neighbor_edge){
                                        found_edges++;
                                    }
                                }

                                // actually_a_face |= (found_edges == 2);
                            }
                            if(!actually_a_face){
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    return false;

}





std::set<EdgeCycle> EdgeCollapser::findNonLinkEdgeCyclesWithAtLeastOneInteriorEdge(){

    storeLinkConditionAsProperty(false);

    //std::cout<<" -------------------------------------------- "<<std::endl;
    //std::cout<<" LOOKING NON-LINK EDGE CYCLES.... "<<std::endl;
    //auto start_time = std::chrono::high_resolution_clock::now();

    std::set<EdgeCycle> cycles;

    for(auto e: mesh_.edges()){
        if(!link_condition_prop_[e]){
            //std::cout<<" -- checking non-link edge "<<e<<": "<<mesh_.edge(e)<<std::endl;

            auto from_vertex = mesh_.edge(e).from_vertex();
            auto to_vertex   = mesh_.edge(e).to_vertex();

            std::vector<EdgeHandle> neighbor_non_link_edges;

            std::vector<std::vector<EdgeHandle>> non_link_edge_cycles;

            for(auto out_he_it = mesh_.outgoing_halfedges(from_vertex).first; out_he_it.valid(); out_he_it++){
                auto neighbor_edge = mesh_.edge_handle(*out_he_it);

                if(!link_condition_prop_[neighbor_edge]){
                    auto tip_vertex = mesh_.to_vertex_handle(*out_he_it);

                    for(auto out_out_he_it = mesh_.outgoing_halfedges(tip_vertex).first; out_out_he_it.valid(); out_out_he_it++){
                        auto neighbor_neighbor_edge = mesh_.edge_handle(*out_out_he_it);

                        if(!link_condition_prop_[neighbor_neighbor_edge] &&
                                mesh_.to_vertex_handle(*out_out_he_it) == to_vertex){

                            cycles.insert({e, neighbor_edge, neighbor_neighbor_edge});

#if 0
                            if(isInterior(e) || isInterior(neighbor_edge) || isInterior(neighbor_neighbor_edge))){
                                cycles.insert({e, neighbor_edge, neighbor_neighbor_edge});
                            }
#endif

                        }
                    }
                }
            }
        }
    }


    //auto end_time = std::chrono::high_resolution_clock::now();
    //std::cout<<" ANALYSIS PERFORMED IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS. RESULTS: "<<std::endl;
    //std::cout<<" ---- non-link edge cycles count: "<<cycles.size()<<std::endl;
    //std::cout<<" -------------------------------------------- "<<std::endl;

    return cycles;

}






std::set<EdgeCycle> EdgeCollapser::findNonLinkInteriorEdgeCyclesAroundVertex(const OpenVolumeMesh::VertexHandle& vertex,
                                                                             bool count_special_config_cycles){


    //std::cout<<" -------------------------------------------- "<<std::endl;
    //std::cout<<" LOOKING NON-LINK EDGE CYCLES.... "<<std::endl;


    //auto start_time = std::chrono::high_resolution_clock::now();

    std::set<EdgeCycle> cycles;

    std::set<EdgeHandle> edges_to_check;

    auto link_vertices = link(mesh_, vertex).vertices();

    for(auto v: link_vertices){
        for(auto out_he_it = mesh_.outgoing_halfedges(v).first; out_he_it.valid(); out_he_it++){
            edges_to_check.insert(mesh_.edge_handle(*out_he_it));
        }
    }

    auto visited = mesh_.request_edge_property<bool>("visited");

    for(auto e: edges_to_check){
        visited[e] = true;
        if(topo_helper_.isInterior(e) && !link_condition(mesh_, e)){
            //std::cout<<" -- checking non-link edge "<<e<<": "<<mesh_.edge(e)<<std::endl;

            auto from_vertex = mesh_.edge(e).from_vertex();
            auto to_vertex   = mesh_.edge(e).to_vertex();

            auto vertex_visited = mesh_.request_vertex_property<EdgeHandle>("visited", EdgeHandle(-1));
            std::vector<EdgeHandle> neighbor_non_link_edges;

            std::vector<std::vector<EdgeHandle>> non_link_edge_cycles;

            for(auto out_he_it = mesh_.outgoing_halfedges(from_vertex).first; out_he_it.valid(); out_he_it++){
                auto neighbor_edge = mesh_.edge_handle(*out_he_it);

                if(!visited[neighbor_edge] &&
                        topo_helper_.isInterior(neighbor_edge) &&
                        !link_condition(mesh_, neighbor_edge)){

                    auto tip_vertex = mesh_.to_vertex_handle(*out_he_it);

                    for(auto out_out_he_it = mesh_.outgoing_halfedges(tip_vertex).first; out_out_he_it.valid(); out_out_he_it++){
                        auto neighbor_neighbor_edge = mesh_.edge_handle(*out_out_he_it);

                        if(!visited[neighbor_neighbor_edge] &&
                                topo_helper_.isInterior(neighbor_neighbor_edge) &&
                                !link_condition(mesh_, neighbor_neighbor_edge) &&
                                mesh_.to_vertex_handle(*out_out_he_it) == to_vertex){

                            bool actually_a_face(false);
                            //check that this is not actually a face
                            for(auto ef_it = mesh_.ef_iter(e); ef_it.valid(); ef_it++){
                                int found_edges(0);

                                for(auto fe_it = mesh_.fe_iter(*ef_it); fe_it.valid(); fe_it++){
                                    if(*fe_it == neighbor_edge ||
                                            *fe_it == neighbor_neighbor_edge){
                                        found_edges++;
                                    }
                                }

                                actually_a_face |= (found_edges == 2);
                            }
                            if(!actually_a_face){
                                bool insert_cycle(true);
                                if(!count_special_config_cycles){
                                    auto link_intersection = link(mesh_, e).intersection(link(mesh_, neighbor_edge).intersection(link(mesh_, neighbor_neighbor_edge)));
                                    insert_cycle = !(link_intersection.edges().size() && topo_helper_.edgeValence(*link_intersection.edges().begin()) == 3);
                                }

                                if(insert_cycle){
                                    cycles.insert({e, neighbor_edge, neighbor_neighbor_edge});
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    //auto end_time = std::chrono::high_resolution_clock::now();

    /*std::cout<<" ANALYSIS PERFORMED IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS. RESULTS: "<<std::endl;
    std::cout<<" ---- non-link edge cycles count: "<<cycles.size()<<std::endl;

    std::cout<<" -------------------------------------------- "<<std::endl;*/

    return cycles;

}

HalfEdgeHandle EdgeCollapser::findNonBoundaryCollapsibleEdgeOnBadTet(const BadTetFinder::BadTetList& bad_tets){



    HalfEdgeHandle edge_to_collapse(-1);

    if(bad_tets.first.size()){
        std::cout<<" - searching for a collapsible edge on a degenerate tet..."<<std::endl;
        for(auto deg_tet: bad_tets.first){
            if(edge_to_collapse.idx() != -1){
                break;
            }

            for(auto ce_it = mesh_.ce_iter(deg_tet); ce_it.valid(); ce_it++){
                if(topo_helper_.isCollapsible(*ce_it)){
                    edge_to_collapse = mesh_.halfedge_handle(*ce_it, 0);
                }
            }
        }
    }else if(bad_tets.second.size()){
        std::cout<<" - searching for a collapsible edge on a flipped tet..."<<std::endl;
        for(auto flipped_tet: bad_tets.second){
            if(edge_to_collapse.idx() != -1){
                std::cout<<" - stopping search"<<std::endl;
                break;
            }

            for(auto ce_it = mesh_.ce_iter(flipped_tet); ce_it.valid(); ce_it++){
                std::cout<<" -- checking edge "<<mesh_.edge(*ce_it)<<std::endl;
                if(topo_helper_.isCollapsible(*ce_it)){
                    edge_to_collapse = mesh_.halfedge_handle(*ce_it, 0);
                    std::cout<<" ---> found edge "<<mesh_.halfedge(edge_to_collapse)<<std::endl;
                }else{
                    std::cout<<" ---> uncollapsible:"<<std::endl;
                    std::cout<<"    link condition : "<<link_condition(mesh_, mesh_.halfedge_handle(*ce_it, 0))<<std::endl;
                    std::cout<<"          boundary : "<<mesh_.is_boundary(*ce_it)<<std::endl;
                    std::cout<<"           valence : "<<mesh_.valence(*ce_it)<<std::endl;
                    std::cout<<" from-vertex boundary : "<<mesh_.is_boundary(mesh_.edge(*ce_it).from_vertex())<<std::endl;
                    std::cout<<"   to-vertex boundary : "<<mesh_.is_boundary(mesh_.edge(*ce_it).to_vertex())<<std::endl;
                }
                std::cout<<" ---------------"<<std::endl;

            }
            std::cout<<" =========================== "<<std::endl;
        }
    }else{
        std::cout<<" - there are no bad tets anymore"<<std::endl;
    }

    return edge_to_collapse;
}



OpenVolumeMesh::HalfEdgeHandle EdgeCollapser::findAnyNonBoundaryCollapsibleEdge(){


    bool found_link_condition(false);
    HalfEdgeHandle edge_to_collapse(-1);
    auto e_it = mesh_.edges().first;

    while(!found_link_condition && e_it != mesh_.edges().second){
        if(topo_helper_.isCollapsible(*e_it)){
            found_link_condition = true;
            edge_to_collapse = mesh_.halfedge_handle(*e_it, 0);
        }
        e_it++;
    }

    return edge_to_collapse;

}




OpenVolumeMesh::HalfEdgeHandle EdgeCollapser::findAnyBoundaryCollapsibleEdge(){


    bool found_link_condition(false);
    HalfEdgeHandle edge_to_collapse(-1);
    auto e_it = mesh_.edges().first;

    while(!found_link_condition && e_it != mesh_.edges().second){
        if(mesh_.is_boundary(*e_it) && link_condition(mesh_, *e_it)){
            found_link_condition = true;
            edge_to_collapse = mesh_.halfedge_handle(*e_it, 0);
        }
        e_it++;
    }

    return edge_to_collapse;

}






void EdgeCollapser::storeLinkConditionAsProperty(bool print_out_status){

    const int total_edge_count = mesh_.n_logical_edges();
    int edge_count(0);
    int percentage_increase(1);
    int next_percentage(percentage_increase);

    auto start_time = std::chrono::high_resolution_clock::now();

    for(auto e: mesh_.edges()){
        link_condition_prop_[e] = link_condition(mesh_, e);

        /*if(!link_condition_prop_[e]){
            non_link_condition_edges_.push_back(e);
        }*/

        edge_count++;

        int percentage_done = (edge_count * 100) / total_edge_count;

        if( percentage_done > next_percentage){
            auto current_time = std::chrono::high_resolution_clock::now();

            auto current_duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
            auto estimated_remaining_time = (current_duration_ms * 100) / percentage_done - current_duration_ms;

            if(print_out_status){
                std::cout<<" checked link condition of "<<percentage_done<<"% of all edges after "<<current_duration_ms <<" milliseconds. ERT: "<<estimated_remaining_time<<std::endl;
            }
            next_percentage = percentage_done + percentage_increase;
        }
    }

    setLinkConditionPropStatus(true);
}


void EdgeCollapser::shuffleEdgeVector(std::vector<EdgeHandle>& edge_vector){

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(edge_vector.begin(), edge_vector.end(), g);
}

bool EdgeCollapser::linkConditionUpToDate() const{
    return link_condition_is_up_to_date_;
}




void EdgeCollapser::setLinkConditionPropStatus(bool up_to_date){
    link_condition_is_up_to_date_ = up_to_date;

}



std::ostream& operator<<(std::ostream& os,
                         const EdgeCollapse& collapse){
    os<<"C(("<<collapse.from_vertex<<"-"<<collapse.to_vertex<<") -> "<<collapse.replacement_vertex<<")";
    return os;
}

std::ostream& operator<<(std::ostream& os,
                         const EdgeSplit& split){


    os<<"S(("<<split.from_vertex<<"-"<<split.to_vertex<<") -> "<<
        "("<<split.from_vertex<<"-"<<split.middle_vertex<<"-"<<split.to_vertex<<"))";

    return os;
}

std::ostream& operator<<(std::ostream& os,
                         const SplitAndCollapse& sac){



    os<<"SAC[ ";
    if(sac.split.from_vertex.idx() == -1 ||
            sac.split.middle_vertex.idx() == -1 ||
            sac.split.to_vertex.idx() == -1){
        os<<"S(-)";
    }else{
        os<<sac.split;
    }
    os<<" | "<<sac.collapse<<" ]";

    return os;
}

}
