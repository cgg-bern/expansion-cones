#include "StarShapifyableExpansionCone.hh"
#include "ProgEmbeddingHelpers.hh"


#define ENABLE_SS_DEBUG 0
#if ENABLE_SS_DEBUG
#define SS_DEBUG(msg) SS_DEBUG(msg<<std::endl);
#define IF_SS_DEBUG(statement) statement;
#else
#define SS_DEBUG(msg)
#define IF_SS_DEBUG(statement)
#endif

#define ENABLE_PRECISION_REDUCTION 1


#define ENABLE_BAD_TETS_CHECKS_OUT_OF_LOOPS 0
#define ENABLE_BAD_TETS_CHECKS_IN_LOOPS 0
#define ENABLE_TIMINGS 0

namespace OpenVolumeMesh{

/*int StarShapifyer::star_shapify(ExpansionCone& expansion_cone,
                                TetrahedralMesh& mesh){
    StarShapifyer star_shapifyer(expansion_cone, mesh);

    return star_shapifyer.star_shapify();
}*/

StarShapifyableExpansionCone::StarShapifyableExpansionCone(const ExpansionCone& expansion_cone,
                                                           const int max_allocated_time_s) :
    ExpansionCone(expansion_cone),
    collapsed_prop_(request_vertex_property<bool>()),
    secondary_vertex_prop_(request_vertex_property<VertexHandle>()),
    witness_vertex_1_ring_cell_prop_(request_cell_property<bool>()),
    witness_vertex_1_ring_base_vertex_prop_(request_vertex_property<bool>()),
    original_edge_prop_(request_edge_property<bool>()),
    mid_vertices_positions_prop_(request_vertex_property<std::vector<VertexPosition>>()),
    mid_vertices_prop_(request_vertex_property<std::vector<VertexHandle>>()),
    original_base_vertex_prop_(request_vertex_property<VertexHandle>("", VertexHandle(-1))),
    //attached_vertices_f_prop_(request_face_property<std::pair<VertexHandle, int>>()),
    //attached_vertices_e_prop_(request_edge_property<std::pair<VertexHandle, int>>()),
    base_vertex_to_spoke_edge_prop_(request_vertex_property<HalfEdgeHandle>()),
    remaining_trimmed_copy_vertices_at_last_RS_iteration_(std::numeric_limits<int>::max()),
    contracted_prop_(request_vertex_property<bool>()),
    max_allocated_time_s_(max_allocated_time_s),
    start_time_(std::chrono::high_resolution_clock::now())
  /*, base_edge_to_spoke_face_prop_(request_edge_property<FaceHandle>())*/{

    primary_reverse_shelling_sequence_.clear();
    full_reverse_shelling_sequence_.clear();
    collapse_list_.clear();
    split_list_.clear();
    vertices_to_contract_.clear();
}



void StarShapifyableExpansionCone::operator=(const StarShapifyableExpansionCone& other){

    ExpansionCone::operator=(other);

    for(auto v: vertices()){
        collapsed_prop_[v] = false;
        secondary_vertex_prop_[v] = VertexHandle(-1);
        witness_vertex_1_ring_base_vertex_prop_[v] = false;
        base_vertex_to_spoke_edge_prop_[v] = HalfEdgeHandle(-1);
    }

    remaining_trimmed_copy_vertices_at_last_RS_iteration_ = std::numeric_limits<int>::max();
    full_reverse_shelling_sequence_ = other.full_reverse_shelling_sequence_;
    collapse_list_ = other.collapse_list_;
    split_list_ = other.split_list_;
    vertices_to_contract_ = other.vertices_to_contract_;
    //max_allocated_time_s_ = other.max_allocated_time_s_;
}



STAR_SHAPIFICATION_RESULT StarShapifyableExpansionCone::is_star_shapifyable(){

    //function deprecated but kept for now to maintain Expander implementation as it is
    return SS_SUCCESS;
}



void StarShapifyableExpansionCone::set_vertex(const VertexHandle& vh,
                                              const VertexPosition& pos,
                                              bool update_max_size){

    if(update_max_size){
        int size = byte_size(pos);
        new_position_max_precision_ = std::max(new_position_max_precision_, size);
        new_position_precision_record_.push_back(size);
    }

    ExpansionCone::set_vertex(vh, pos);
}


STAR_SHAPIFICATION_RESULT StarShapifyableExpansionCone::star_shapify(bool use_faster_but_unsafe_RS,
                                                                     bool reduce_precision){

    SS_DEBUG(" =============================================================="<<std::endl);
    SS_DEBUG(" =============================================================="<<std::endl);
    SS_DEBUG(" =============== STAR-SHAPIFYING CONE "<<(*this)<<"..."<<std::endl);

#if !ENABLE_PRECISION_REDUCTION
    SS_DEBUG(" WARNING - disabled precision reduction"<<std::endl);
    reduce_precision = false;
#endif


    if(!full_reverse_shelling_sequence_.empty()){
        SS_DEBUG(" ERROR - reverse shelling sequence is not empty at 'star_shapify' start"<<std::endl);
        return SS_ERROR;
    }


    if(!split_list_.empty()){
        SS_DEBUG(" ERROR - split list is not empty at 'star_shapify' start"<<std::endl);
        return SS_ERROR;
    }


    if(!collapse_list_.empty()){
        SS_DEBUG(" ERROR - collapse list is not empty at 'star_shapify' start"<<std::endl);
        return SS_ERROR;
    }

    /*if(use_faster_but_unsafe_RS){
        SS_DEBUG(" ---------- original cone details: "; print_details();
    }*/


    //TBD if this parameter makes sense.
    //star_shapifying_cluster = false;

    IF_SS_DEBUG(auto original_cone = *this);


#if ENABLE_BAD_TETS_CHECKS_OUT_OF_LOOPS
    //degenerate tets are allowed when star-shapifying a cluster because
    //tets incident to multiple tip vertices are bound to be degenerate
    if(/*!star_shapifying_cluster &&
            */ExactBadTetFinder::meshContainsDegenerateTets(*this, vertex_position_prop_)){
        SS_DEBUG(" =================================================="<<std::endl);
        SS_DEBUG(" ERROR - cone initially contains degenerate tets: "<<std::endl);
        auto bad_tets = ExactBadTetFinder::findBadTets(*this, vertex_position_prop_);
        for(auto t: bad_tets.first){
            SS_DEBUG(" - "<<t<<" -> "<<cone_to_mesh_handle(t)<<": "<<get_cell_vertices(t)<<std::endl);
        }
        print_details();
        SS_DEBUG(" =================================================="<<std::endl);
        return SS_ERROR;
    }
#endif

    const int initial_n_vertices(n_vertices());

    //to ensure that cells are deleted when deleting vertices
    this->enable_fast_deletion(false);

    //1. identify best witness vertex
    //TODO: put this in a function
    int max_c_valence(0);
    witness_vertex_ = VertexHandle(-1);

    for(auto v: vertices()){
        if(!is_cone_tip(v)){
            auto c_val = cell_valence(v);
            if(c_val > max_c_valence){
                max_c_valence = c_val;
                witness_vertex_ = v;
            }
        }
    }


    SS_DEBUG(" - witness vertex: "<<witness_vertex_<<" with valence = "<<valence(witness_vertex_));
    if(witness_vertex_.idx() == -1){
        SS_DEBUG(" ERROR - couldn't find the witness vertex"<<std::endl);
        return SS_ERROR;
    }



    //2. split the base to remove any tri-tet
    int base_split_result = split_base_non_link_edges_not_connected_to_witness();
    if(base_split_result){
        SS_DEBUG(" ERROR - couldn't split base non-link edges"<<std::endl);
        return SS_ERROR;
    }

    SS_DEBUG(" ---------- cone after base clean-up details: ");
    IF_SS_DEBUG(print_details(2));


    //mark the witness vertex' 1-ring as such
    for(auto vc_it = vc_iter(witness_vertex_); vc_it.valid(); vc_it++){
        witness_vertex_1_ring_cell_prop_[*vc_it] = true;
        SS_DEBUG(" -- set cell "<<get_cell_vertices(*vc_it)<<" as part of witness 1-ring");
    }
    for(auto vv_it = vv_iter(witness_vertex_); vv_it.valid(); vv_it++){
        witness_vertex_1_ring_base_vertex_prop_[*vv_it] = !tip_vertices_prop_[*vv_it];
    }

    //2nd pass
    //4. unlock the 2-ring neighborhood vertices connected to >1 1-ring vertices
    /*base_split_result = split_base_1_2_ring_edges();
    if(base_split_result){
        SS_DEBUG(" ERROR - couldn't unlock the 2-ring vertices"<<std::endl);
        return SS_ERROR;
    }*/


    //3. find the initial candidates to remove
    std::vector<VertexHandle> vertices_to_collapse;
    gather_vertices_outside_of_1_ring_neighborhood_of_witness_vertex(vertices_to_collapse);
    SS_DEBUG(" - vertices to remove ("<<vertices_to_collapse.size()<<"): "<<vertices_to_collapse);
    if(!vertices_to_collapse.size()){
        SS_DEBUG(" ERROR - couldn't find any vertices to collapse"<<std::endl);
        return SS_NO_VERTICES_TO_COLLAPSE;
    }

    for(auto e: edges()){
        original_edge_prop_[e] = true;
    }

    //5. run the reverse shelling sequence to split any necessary base edges
    const int max_iteration_count(n_edges() - valence(witness_vertex_));


    auto reverse_shelling_result = run_reverse_shelling_sequence(max_iteration_count,
                                                                 use_faster_but_unsafe_RS,
                                                                 vertices_to_collapse);



    if(reverse_shelling_result){
        SS_DEBUG(" ERROR - couldn't run the reverse shelling sequence"<<std::endl);
        SS_DEBUG(" - original cone details:");
        //IF_SS_DEBUG(original_cone.print_details());
        return SS_ERROR;
    }



    //SS_DEBUG(" WARNING - MANUAL STOP HERE"<<std::endl);
    //return SS_GENERAL_FAILURE;


    //6. reset the trimmed copy
    trimmed_copy_ = *this;
    trimmed_copy_.enable_fast_deletion(false);


    //7. set-up the spoke tip/base vertex -> spoke property
    for(auto tip_v: cone_tip_vertices()){
        for(auto out_he: outgoing_halfedges(tip_v)){
            base_vertex_to_spoke_edge_prop_[to_vertex_handle(out_he)] = out_he;
            original_base_vertex_prop_[to_vertex_handle(out_he)] = to_vertex_handle(out_he);
            /*SS_DEBUG(" - set edge "<<cone_to_mesh_handle(from_vertex_handle(out_he))<<
                       "-"<<cone_to_mesh_handle(to_vertex_handle(out_he))<<
                       " as spoke of base vertex "<<expansion_cone.cone_to_mesh_handle(to_vertex_handle(out_he))<<std::endl);
            */
        }
    }

    //8. update the list of candidates to collapse (necessary because the indices of the new vertices
    //   might not match between the trimmed copy and the original cone)
    gather_vertices_outside_of_1_ring_neighborhood_of_witness_vertex(vertices_to_collapse);
    SS_DEBUG(" - vertices to collapse ("<<vertices_to_collapse.size()<<"): "<<vertices_to_collapse);
    if(!vertices_to_collapse.size()){
        SS_DEBUG(" ERROR - couldn't find any vertices to collapse after running reverse-shelling sequence"<<std::endl);
        return SS_NO_VERTICES_TO_COLLAPSE;
    }


#if 0
    //temp check
    ExpansionCone updated_cone;
    set_up_1_ring_neighborhood_as_expansion_cone(*this,
                                                 vertex_position_prop_,
                                                 *cone_tip_vertices().begin(),
                                                 updated_cone);
    VertexPosition pos;
    auto geo_exp_result = updated_cone.is_geo_expandable(pos);

    if(geo_exp_result){
        SS_DEBUG(" ERROR - temp check for geo-exp failed before the collapse sequence"<<std::endl);
        SS_DEBUG(" updated cone: "<<std::endl);
        updated_cone.print_details();
        SS_DEBUG(" cone: "<<std::endl);
        print_details();
        return SS_ERROR;
    }
#endif



    //9. collapse all necessary vertices in reverse-shelling order
    auto collapse_result = run_collapse_sequence(vertices_to_collapse);
    if(collapse_result){
        SS_DEBUG(" - failed to run the collapse sequence"<<std::endl);
        return SS_ERROR;
    }

    //10. and run the contraction sequence to have all positive-volume tets
    auto contraction_result = run_contraction_sequence(reduce_precision);
    if(contraction_result){
        SS_DEBUG(" - failed to run the contraction sequence"<<std::endl);
        return SS_ERROR;
    }

    SS_DEBUG(" .... DONE!"<<std::endl);
//#if ENABLE_BAD_TETS_CHECKS_OUT_OF_LOOPS
    auto bad_tets = ExactBadTetFinder::findBadTets(*this, vertex_position_prop_);
    if(bad_tets.first.size()){
        SS_DEBUG(" ERROR - cone contains degenerate tets after contraction:"<<std::endl);
        for(auto bad_tet: bad_tets.first){
            SS_DEBUG(" - "<<bad_tet<<" :" <<get_cell_vertices(bad_tet)<<std::endl);
        }
        return SS_ERROR;
    }
    if(bad_tets.second.size()){
        SS_DEBUG(" ERROR - cone contains flipped tets after contraction:"<<std::endl);
        for(auto bad_tet: bad_tets.second){
            SS_DEBUG(" - "<<bad_tet<<" :" <<get_cell_vertices(bad_tet)<<std::endl);
        }
        return SS_ERROR;
    }
//#endif

    if(found_deg_tets_){
        SS_DEBUG(" ERROR - deg tets created during contraction but now they're gone..."<<std::endl);
        return SS_ERROR;
    }


    int original_edges_count(0), new_edges_count(0);
    count_base_edges(original_edges_count, new_edges_count);


    SS_DEBUG(" DONE WITH STAR-SHAPIFICATION"<<std::endl);
    SS_DEBUG(" ---------------------------------------"<<std::endl);
    SS_DEBUG(" initial vertex count = "<<initial_n_vertices<<std::endl);
    SS_DEBUG("   final vertex count = "<<n_logical_vertices()<<std::endl);
    SS_DEBUG("         growth ratio = "<<((double)n_logical_vertices() / (double)initial_n_vertices)<<std::endl);
    SS_DEBUG(" ---------------------------------------"<<std::endl);
    SS_DEBUG(" original edges (pieces) count = "<<original_edges_count<<std::endl);
    SS_DEBUG("      new edges (pieces) count = "<<new_edges_count<<std::endl);
    SS_DEBUG("            edges growth ratio = "<<((double)new_edges_count / (double)original_edges_count)<<std::endl);
    SS_DEBUG(" ---------------------------------------"<<std::endl);
    SS_DEBUG("   max position byte precision = "<<new_position_max_precision_<<std::endl);
    SS_DEBUG("      new position record size = "<<new_position_precision_record_.size()<<std::endl);
    SS_DEBUG("             total saved bytes = "<<saved_position_bytes_<<std::endl);
    SS_DEBUG("  10-points approximation of precision evolution: ");
    double sum(0);
    int mod = std::max((int)new_position_precision_record_.size() / 10, 1);
    for(int i(0); i<(int)new_position_precision_record_.size(); i++){
        sum += new_position_precision_record_[i];
        if(i && !(i % mod)){
            SS_DEBUG(" "<<(sum/mod));
            sum = 0;
        }
    }
    SS_DEBUG(std::endl);
    SS_DEBUG(" =============================================================="<<std::endl);

    return SS_SUCCESS;
}


int StarShapifyableExpansionCone::find_split_path_along_RS_sequence(const EdgeHandle& start,
                                                                    const EdgePropertyT<EdgeHandle>& collapsed_to_prop,
                                                                    std::deque<EdgeHandle>& path_to_split){


    //SS_DEBUG("------------------------------------------------------------------");
    //SS_DEBUG("----- looking for split path starting from edge "<<edge(start));

    path_to_split.clear();

    auto current_edge = start;

    int i(0);
    while(current_edge.is_valid() && i<(int)trimmed_copy_.n_edges()){

        //SS_DEBUG(" ----------------- split-path back-tracking iteration "<<(i+1));
        //SS_DEBUG(" - current edge "<<edge(current_edge));

        path_to_split.push_back(current_edge);

        current_edge = collapsed_to_prop[current_edge];
        i++;
    }



    //SS_DEBUG("------------------------------------------------------------------");

    return 0;
}


int StarShapifyableExpansionCone::find_shortest_dual_path_to_boundary(const EdgeHandle& start,
                                                                       std::deque<EdgeHandle>& path_to_split){


    //SS_DEBUG("------------------------------------------------------------------");
    //SS_DEBUG("----- looking for shortest dual path from edge "<<edge(start));


    //find the 2-ring blocking vertex
    //auto blocking_2_ring_vertex = witness_vertex_1_ring_prop[from_vertex_handle(start)] ? to_vertex_handle(start).to_vertex() : edge(start).from_vertex();
    //auto blocking_2_ring_vertex = to_vertex_handle(start);
    //SS_DEBUG(" - identified 2-ring blocking vertex as "<<blocking_2_ring_vertex);



    //auto from_vertex = edge(start).from_vertex();
    //auto to_vertex = edge(start).to_vertex();

    path_to_split.clear();

    if(is_base_boundary_edge(start)){
        SS_DEBUG(" --> edge "<<edge(start)<<" is an original base boundary tet -> done");
        path_to_split.push_back(start);
        return 0;
    }


#warning TODO: use prop instead
    auto distance_prop = request_cell_property<int>("", std::numeric_limits<int>::max());
    std::vector<CellHandle> to_visit;
    for(auto vc_it = vc_iter(witness_vertex_); vc_it.valid(); vc_it++){
        distance_prop[*vc_it] = 0;
        auto c_vertices = get_cell_vertices(*vc_it);

        //SS_DEBUG(" - cell "<<get_cell_vertices(*vc_it)<<" is incident to the witness (distance 0)");
    }

    //setting distance of all remaining tets of the trimmed copy to 0 so we don't go through those
    for(auto c: trimmed_copy_.cells()){
        //SS_DEBUG(" - removed tet "<<get_cell_vertices(c)<<" from path candidates"<<std::endl);
        distance_prop[c] = 0;
    }

    //pick whichever incident cell is not part of the trimmed copy
    CellHandle start_tet(-1);
    for(auto ec_it = ec_iter(start); ec_it.valid(); ec_it++){
        SS_DEBUG(" - checking incident tet "<<get_cell_vertices(*ec_it)<<std::endl);
        if(distance_prop[*ec_it]){
            start_tet = *ec_it;
        }
    }
    if(!start_tet.is_valid()){
        SS_DEBUG(" ERROR - couldn't find a starting tet for the shortest dual path-finding algorithm"<<std::endl);
        return -1;
    }

    distance_prop[start_tet] = 1;
    to_visit.push_back(start_tet);

    SS_DEBUG(" - starting tet: "<<get_cell_vertices(start_tet));


    CellHandle closest_boundary_tet(-1);
    int current_distance(0);
    auto previous_cell_on_path = request_cell_property<CellHandle>("", CellHandle(-1));
    while(closest_boundary_tet.idx() == -1 &&
          !to_visit.empty() &&
          current_distance < (int)n_cells()){//hard limit to avoid infinite loops

        //SS_DEBUG(" ------");
        //SS_DEBUG(" - visiting  ring "<<(current_distance+1)<<" with "<<to_visit.size()<<" cells to visit");


        std::vector<CellHandle> next_ring;
        for(auto c: to_visit){
            int neighbors_count(0);
            for(auto cc_it = cc_iter(c); cc_it.valid(); cc_it++){
                if(distance_prop[*cc_it] > current_distance){
                    distance_prop[*cc_it] = current_distance + 1;
                    next_ring.push_back(*cc_it);
                    previous_cell_on_path[*cc_it] = c;
                   // SS_DEBUG(" -- cell "<<get_cell_vertices(*cc_it)<<" is on next ring, set previous cell as "<<get_cell_vertices(previous_cell_on_path[*cc_it]));
                }
                neighbors_count++;
            }
            if(neighbors_count < 3){
                //SS_DEBUG(" ---> found boundary tet "<<get_cell_vertices(c)<<", reached end of path");
                closest_boundary_tet = c;
                break;
            }
        }
        to_visit = next_ring;
        current_distance++;
    }

    if(closest_boundary_tet.idx() == -1){
        SS_DEBUG(" ERROR - couldn't find closest boundary tet"<<std::endl);
        return -1;
    }

    //first, get the boundary edge
    for(auto ce_it = ce_iter(closest_boundary_tet); ce_it.valid(); ce_it++){
        //SS_DEBUG(" - checking edge "<<edge(*ce_it)<<" with valence "<<valence(*ce_it));
        if(valence(*ce_it) == 2 &&
                !tip_vertices_prop_[edge(*ce_it).from_vertex()] &&
                !tip_vertices_prop_[edge(*ce_it).to_vertex()]){
            //SS_DEBUG(" --> found boundary edge "<<edge(*ce_it));
            path_to_split.push_back( *ce_it);
            break;
        }

    }

    //Then get edges to split from going back from the boundary to the starting tet
    //SS_DEBUG(" - backtracking from tet "<<get_cell_vertices(closest_boundary_tet)<<" to get full path");
    auto current_tet = closest_boundary_tet;
    int i(0);
    while(current_tet != start_tet &&
          i < current_distance){
        //SS_DEBUG(" -------- ");
        //SS_DEBUG(" - current tet: "<<get_cell_vertices(current_tet));

        auto previous_tet = previous_cell_on_path[current_tet];


        auto common_edge = common_base_edge(current_tet, previous_tet);
        if(common_edge.idx() == -1){
            SS_DEBUG(" ERROR - couldn't find common edge between tets "<<get_cell_vertices(current_tet)<<
                       " and "<<get_cell_vertices(previous_tet)<<std::endl);
            return -1;
        }

        path_to_split.push_front(common_edge);

        current_tet = previous_tet;

        i++;
    }

    //path_to_split.push_back(start);
    if(path_to_split.front() != start){
        SS_DEBUG("WARNING - added start edge 'manually'");
        path_to_split.push_front(start);
    }

    //SS_DEBUG("----- done. Path length = "<<path_to_split.size());
    //SS_DEBUG("------------------------------------------------------------------");

    return 0;
}



EdgeHandle StarShapifyableExpansionCone::common_base_edge(const CellHandle& c1,
                                                          const CellHandle& c2) const{

    EdgeHandle common_edge(-1);

    HalfFaceHandle common_face(-1);
    //find common face between two tets
    for(auto chf_it = chf_iter(c1); chf_it.valid(); chf_it++){
        auto op_hf = opposite_halfface_handle(*chf_it);
        if(!is_boundary(op_hf)){
            auto neighbor_tet = incident_cell(op_hf);
            if(neighbor_tet == c2){
                //SS_DEBUG(" -- found common face with previous tet: "<<get_halfface_vertices(op_hf));
                common_face = op_hf;
                break;
            }
        }
    }
    if(common_face.idx() == -1){
        SS_DEBUG(" ERROR - couldn't find common face between tets "<<get_cell_vertices(c1)<<
                   " and "<<get_cell_vertices(c2)<<std::endl);
        return common_edge;
    }

    //find common edge between two tets
    for(auto hfe_it = hfe_iter(common_face); hfe_it.valid(); hfe_it++){
        if(!tip_vertices_prop_[edge(*hfe_it).from_vertex()] &&
                !tip_vertices_prop_[edge(*hfe_it).to_vertex()]){
            common_edge = *hfe_it;
            //SS_DEBUG(" -- found common edge "<<edge(common_edge));
            break;
        }
    }

    return common_edge;
}

int StarShapifyableExpansionCone::split_base_non_link_edges_not_connected_to_witness(){
    SS_DEBUG(" ------------------------------------");
    SS_DEBUG(" splitting non-link base edge to remove tri-tets...");

    ExpansionCone cone_base = *this;
    for(auto tip_v: cone_tip_vertices()){
        cone_base.delete_vertex(tip_v);
    }

    int i(0);
    EdgeHandle to_split(-1);
    do{
        to_split = EdgeHandle(-1);

        for(auto e: cone_base.edges()){
            auto from_v = cone_base.edge(e).from_vertex();
            auto to_v   = cone_base.edge(e).to_vertex();

            if(from_v != witness_vertex_ &&
                    to_v != witness_vertex_ &&
                    !link_condition(cone_base, e)){
                to_split = e;
                break;
            }

        }

        if(to_split.idx() != -1){

            auto cone_mid_vertex = cone_base.split_edge(to_split);
            auto mid_vertex = split_base_edge_and_register_split(to_split);
            cone_base.set_vertex(cone_mid_vertex, vertex(mid_vertex));

            SS_DEBUG(" --> split "<<edge(to_split)<<" -> "<<mid_vertex<<" at "<<vec2vec(vertex(mid_vertex)));
        }


        i++;
    }while(to_split.idx() != -1 && i<(int)n_edges());

    if(i-1){
        SS_DEBUG(" fully split base after "<<(i-1)<<" splits");
    }


    SS_DEBUG(" ------------------------------------");
    return 0;
}


int StarShapifyableExpansionCone::split_base_1_2_ring_edges(){

    SS_DEBUG(" ------------------------------------");
    SS_DEBUG(" splitting 1-2-ring edges from base");

    std::vector<VertexHandle> one_ring_vertices;
    for(auto v: vertices()){
        if(witness_vertex_1_ring_base_vertex_prop_[v]){
            one_ring_vertices.push_back(v);
        }
    }
    SS_DEBUG(" - 1-ring vertices: "<<one_ring_vertices);


    std::vector<EdgeHandle> one_two_ring_edges;
    auto edge_to_sub_edge_prop = request_edge_property<EdgeHandle>("",EdgeHandle(-1));
    for(auto one_ring_v: one_ring_vertices){
        for(auto out_he: /*trimmed_copy_.*/outgoing_halfedges(one_ring_v)){
            auto neighbor = to_vertex_handle(out_he);
            if(!witness_vertex_1_ring_base_vertex_prop_[neighbor] &&
                    !tip_vertices_prop_[neighbor] &&
                    neighbor != witness_vertex_){
                one_two_ring_edges.push_back(edge_handle(out_he));
                edge_to_sub_edge_prop[edge_handle(out_he)] = edge_handle(out_he);
                SS_DEBUG(" - found 1-2-ring edge "<<halfedge(out_he))
            }
        }
    }

    //TODO: check if possible to remove (any) two of them

    for(auto one_two_ring_e: one_two_ring_edges){

        SS_DEBUG("------------------------------------");
        SS_DEBUG(" -- splitting edge "<<edge(one_two_ring_e));

        //then split this 1-2-ring edge
        auto start_edge = one_two_ring_e;
        auto next_edge = edge_to_sub_edge_prop[one_two_ring_e];
        int i(0);
        while(start_edge != next_edge && i<(int)n_edges()){
            start_edge = next_edge;
            next_edge = edge_to_sub_edge_prop[next_edge];
            SS_DEBUG(" -- next edge: "<<next_edge);
            i++;
        }
        if(i == (int)n_edges()){
            SS_DEBUG(" ERROR - this shouldn't happen"<<std::endl);
            return -1;
        }
        SS_DEBUG(" -- using starting edge "<<edge(start_edge));

        if(start_edge.idx() == -1){
            SS_DEBUG(" ERROR - no sub-edge for 1-2-ring edge "<<edge(one_two_ring_e)<<std::endl);
            return -1;
        }
        if(is_deleted(start_edge)){
            SS_DEBUG(" ERROR - sub-edge "<<edge(start_edge)<<" to original edge "<<edge(one_two_ring_e)<<" is deleted"<<std::endl);
            return -1;
        }

        std::deque<EdgeHandle> path_to_split;
        auto sp_result = find_shortest_dual_path_to_boundary(start_edge,
                                                             path_to_split);

        if(sp_result){
            SS_DEBUG(" ERROR - couldn't find path to split starting from edge "<<edge(start_edge)<<std::endl);
            return SS_ERROR;
        }

        SS_DEBUG(" -- applying "<<path_to_split.size()<<" split path");
        for(auto to_split: path_to_split){

            split_base_edge_and_register_split(to_split);
            auto mid_vertex = split_list_.back().cone_new_vertex;
            SS_DEBUG(" --- split edge "<<edge(to_split));

            //replace the sub-edge if necessary
            if(edge_to_sub_edge_prop[to_split].idx() != -1){
                SS_DEBUG(" --> edge "<<edge(to_split)<<" was split, replacing its sub-edge");
                //NOTE. using random sub-edge for now
                auto sub_edge = edge_handle(find_halfedge(mid_vertex, edge(to_split).from_vertex()));
                if(sub_edge.idx() == -1){
                    SS_DEBUG(" ERROR - couldn't find sub-edge "<<mid_vertex<<"-"<<edge(to_split).from_vertex()<<" to edge "<<edge(to_split)<<std::endl);
                    return -1;
                }
                SS_DEBUG(" --> found sub-edge "<<edge(sub_edge));
                edge_to_sub_edge_prop[to_split] = sub_edge;
                edge_to_sub_edge_prop[sub_edge] = sub_edge;
            }
        }
    }


    SS_DEBUG(" ------------------------------------");
    return 0;
}


int StarShapifyableExpansionCone::split_base_non_extended_link_edges_not_connected_to_witness(){

    SS_DEBUG(" ------------------------------------");
    SS_DEBUG(" splitting non extended link 1-2-ring edges from base");

    auto one_ring_valence_prop = request_vertex_property<int>();
    std::vector<VertexHandle> one_ring_vertices;
    for(auto v: vertices()){
        if(witness_vertex_1_ring_base_vertex_prop_[v]){
            one_ring_vertices.push_back(v);
        }else if(!witness_vertex_1_ring_base_vertex_prop_[v] &&
                 v != witness_vertex_ &&
                 !tip_vertices_prop_[v]){
            SS_DEBUG(" - found non 1-ring, non witness, non tip vertex "<<v)
            for(auto vv_it = vv_iter(v); vv_it.valid(); vv_it++){
                one_ring_valence_prop[v] += witness_vertex_1_ring_base_vertex_prop_[*vv_it];
            }
            SS_DEBUG(" --> 1-ring valence = "<<one_ring_valence_prop[v]);
        }
    }

    SS_DEBUG(" - 1-ring vertices: "<<one_ring_vertices);

    int i(0);
    int max_iteration_count(n_edges() -(int)valence(*cone_tip_vertices_.begin()));

    SS_DEBUG(" max iteration count: "<<max_iteration_count<<std::endl);
    bool found_edge_to_split(true);
    while(found_edge_to_split && i<max_iteration_count){
        found_edge_to_split = false;

        for(auto one_ring_v: one_ring_vertices){
            SS_DEBUG(" - checking neighbors of 1-ring vertex "<<one_ring_v);
            for(auto out_he: outgoing_halfedges(one_ring_v)){
                const auto two_ring_v = to_vertex_handle(out_he);
                SS_DEBUG(" -- 1-ring vertex: "<<one_ring_v<<": checking 2-ring vertex "<<two_ring_v<<" with 1-ring valence "<<one_ring_valence_prop[two_ring_v]);
                if(one_ring_valence_prop[two_ring_v] > 1){
                    SS_DEBUG(" --- 2-ring vertex "<<(two_ring_v)<<" has 1-ring valence "<<one_ring_valence_prop[two_ring_v]);
                    //check if it's not 2 because there's a triangle (a,b,*vv_it) where a and b are part of the 1-ring

                    if(one_ring_valence_prop[two_ring_v] == 2){
                        SS_DEBUG(" ---- checking if not simply on a 1-2-ring face");
                        auto one_ring_neighborhood_prop = request_vertex_property<bool>();
                        std::vector<VertexHandle> one_ring_neighbors;
                        for(auto vv_it2 = vv_iter(two_ring_v); vv_it2.valid(); vv_it2++){
                            SS_DEBUG(" ----- checking 2-ring vertex "<<two_ring_v<<" neighbor "<<(*vv_it2));
                            if(witness_vertex_1_ring_base_vertex_prop_[*vv_it2]){
                                one_ring_neighborhood_prop[*vv_it2] = true;
                                one_ring_neighbors.push_back(*vv_it2);
                            }
                        }
                        if(one_ring_neighbors.size() != 2){
                            SS_DEBUG(" ERROR - 1-ring neighbors to 2-ring vertex "<<two_ring_v<<" are "<<one_ring_neighbors<<" but prop says there should only be 2"<<std::endl);
                            return -1;
                        }

                        for(auto vv_it2 = vv_iter(one_ring_neighbors[0]); vv_it2.valid(); vv_it2++){
                            if(one_ring_neighborhood_prop[*vv_it2]){
                                SS_DEBUG(" -----> triangle ("<<two_ring_v<<","<<(*vv_it2)<<","<<one_ring_neighbors[0]<<") exists, skipping 2-ring vertex "<<two_ring_v);
                                break;
                            }
                        }
                    }


                    SS_DEBUG(" ----> couldn't find 1-2-ring face incident to 2-ring vertex "<<two_ring_v<<", splitting one of its incident 1-2-ring edges");

                    //then split this 1-2-ring edge
                    auto start_edge = edge_handle(out_he);
                    SS_DEBUG(" -- using starting edge "<<edge(start_edge));

                    std::deque<EdgeHandle> path_to_split;
                    auto sp_result = find_shortest_dual_path_to_boundary(start_edge,
                                                                         path_to_split);

                    if(sp_result){
                        SS_DEBUG(" ERROR - couldn't find path to split starting from edge "<<edge(start_edge)<<std::endl);
                        return SS_ERROR;
                    }

                    SS_DEBUG(" -- applying "<<path_to_split.size()<<" split path");
                    for(auto to_split: path_to_split){

                        //decreasing the 1-ring valence of the two ring vertex
                        auto from_v = edge(to_split).from_vertex();
                        auto to_v   = edge(to_split).to_vertex();
                        if(witness_vertex_1_ring_base_vertex_prop_[from_v] && one_ring_valence_prop[to_v] > 0){
                            one_ring_valence_prop[to_v]--;
                            SS_DEBUG(" --> decreased 1-ring valence for 2-ring vertex "<<to_v<<" to "<<one_ring_valence_prop[to_v]);
                        }else if(witness_vertex_1_ring_base_vertex_prop_[to_v] && one_ring_valence_prop[from_v] > 0){
                            one_ring_valence_prop[from_v]--;
                            SS_DEBUG(" --> decreased 1-ring valence for 2-ring vertex "<<from_v<<" to "<<one_ring_valence_prop[from_v]);
                        }

                        split_base_edge_and_register_split(to_split);
                        auto mid_vertex = split_list_.back().cone_new_vertex;
                    }
                    found_edge_to_split = true;
                    break;
                }
            }
        }
        i++;
    }

    if(i == max_iteration_count){
        SS_DEBUG(" ERROR - this shouldn't happen"<<std::endl);
        return -1;
    }

    SS_DEBUG(" ------------------------------------");
    SS_DEBUG(" ------------------------------------");
    return 0;
}



int StarShapifyableExpansionCone::run_reverse_shelling_sequence(const int max_iteration_count,
                                                                bool backtrack_splits_along_shortest_path,
                                                                std::vector<VertexHandle>& vertices_to_collapse){


    SS_DEBUG(" -----------------------------------------------"<<std::endl);
    SS_DEBUG(" ------ RUNNING REVERSE-SHELLING SEQUENCE ------"<<std::endl);
    SS_DEBUG(" -----------------------------------------------"<<std::endl);

    int reverse_shelling_result(0);
    int iteration_count(0);


    max_vertex_index_at_last_RS_iterations_ = -1;

    do{
        SS_DEBUG(" - running iteration "<<iteration_count<<"/"<<max_iteration_count<<" of trying to reverse-shell")
        SS_DEBUG(" - updated vertices count = "<<n_vertices());
        SS_DEBUG(" - updated vertices to collapse count = "<<vertices_to_collapse.size());

        trimmed_copy_ = *this;
        trimmed_copy_.enable_fast_deletion(false);
        reverse_shelling_result = run_reverse_shelling_sequence_internal(backtrack_splits_along_shortest_path,
                                                                         vertices_to_collapse);
        if(reverse_shelling_result == -1){
            SS_DEBUG(" - an error occured while running the reverse-shelling sequence"<<std::endl);
            return SS_ERROR;
        }


        iteration_count++;

        SS_DEBUG(" ------------------- stats after "<<iteration_count<<" RS iterations: ");
        SS_DEBUG(" (witness vertex is "<<witness_vertex_<<" with valence "<<valence(witness_vertex_)<<")");
        SS_DEBUG(" ------------------------------------------------------------");
        SS_DEBUG("                    vertices to collapse: "<<vertices_to_collapse.size());
        SS_DEBUG("                     trimmed cone #verts: "<<trimmed_copy_.n_logical_vertices());
        SS_DEBUG("     remaining base vertices to collapse: "<<(trimmed_copy_.n_logical_vertices() - valence(witness_vertex_) -1));
        SS_DEBUG("                                  status: "<<reverse_shelling_result<<std::endl);
        SS_DEBUG(" ------------------------------------------------------------");
    }while(reverse_shelling_result &&
           iteration_count < max_iteration_count);


    if(!reverse_shelling_result){
        SS_DEBUG(" - SUCCESFULLY RAN REVERSE-SHELLING SEQUENCE AFTER "<<iteration_count<<" ATTEMPTS USING "<<(backtrack_splits_along_shortest_path ? "SHORTEST-PATH" : "SAFE")<<" METHOD "<<std::endl);
        SS_DEBUG(" - original cone after reverse-shelling sequence:");
        IF_SS_DEBUG(print_details());
        SS_DEBUG(" ===============================================================");
    }

    return reverse_shelling_result;
}


int StarShapifyableExpansionCone::run_reverse_shelling_sequence_internal(bool backtrack_splits_along_shortest_path,
                                                                         std::vector<VertexHandle>& vertices_to_collapse){

    SS_DEBUG(" -----------------------------------------------------------------------------------");
    SS_DEBUG(" ----------------------- RUNNING REVERSE-SHELLING SEQUENCE -------------------------");
    SS_DEBUG(" -----------------------------------------------------------------------------------");
    SS_DEBUG(" --> back-tracking splits along shortest path: "<<backtrack_splits_along_shortest_path);

    //used as hard limit for valence reduction
    //const int initial_n_edges(n_edges());

    splits_count_ = 0;
    backtracking_splits_count_ = 0;

    auto cell_valence_at_deletion_time_prop = trimmed_copy_.request_vertex_property<int>("c_val_at_deletion");
    trimmed_copy_.set_persistent(cell_valence_at_deletion_time_prop);

    SS_DEBUG(" - current reverse-shelling sequence: "<<full_reverse_shelling_sequence_);


    //maps an edge to another edge that was collapsed onto it.
    auto collapsed_edge_prop = request_edge_property<EdgeHandle>("", EdgeHandle(-1));
    auto collapse_depth_prop = request_edge_property<int>();

    primary_reverse_shelling_sequence_ = full_reverse_shelling_sequence_;

    //removing the reverse-shelling sequence to prioritize the new boundary vertices
    if(!backtrack_splits_along_shortest_path){
        for(auto v: primary_reverse_shelling_sequence_){
            collapsed_prop_[v] = true;
        }
    }

    int current_sequence_index(0);

    full_reverse_shelling_sequence_.clear();

    int base_collapse_count(0);
    int removed_vertices_count(0);
    while(removed_vertices_count < (int)vertices_to_collapse.size()){

        check_for_timeout();

        SS_DEBUG(" ================= removing "<<(removed_vertices_count+1)<<"-th vertex...");

        if(removed_vertices_count && !(removed_vertices_count % 1000)){
            SS_DEBUG(" - removed "<<removed_vertices_count<<
                       "/"<<vertices_to_collapse.size()<<
                       " vertices, diff = "<<(vertices_to_collapse.size() - removed_vertices_count)<<std::endl);
        }

        VertexHandle vertex_to_collapse(-1);
        int c_val(0);

        if(!backtrack_splits_along_shortest_path){
            int min_index(-1);
            if(current_sequence_index < (int)primary_reverse_shelling_sequence_.size()){
                SS_DEBUG(" - checking if next vertex in sequence "<<primary_reverse_shelling_sequence_[current_sequence_index]<<" can be removed");
                c_val = trimmed_copy_.cell_valence(primary_reverse_shelling_sequence_[current_sequence_index]);
                if(c_val == 1 || c_val == 2){
                    vertex_to_collapse = primary_reverse_shelling_sequence_[current_sequence_index];
                    //c_val = trimmed_copy_.cell_valence(vertex_to_collapse);
                    //SS_DEBUG(" --> next vertex in sequence "<<vertex_to_collapse<<" has valence "<<c_val<<", using this one and updating the sequence index to "<<(current_sequence_index+1));
                    current_sequence_index++;
                }
                SS_DEBUG(" - sequence not re-done, using max vertex index before split-backpropagation "<<
                         max_vertex_index_at_last_RS_iterations_<<
                         " as min index for new vertex to remove"<<std::endl);

                min_index = max_vertex_index_at_last_RS_iterations_;
            }else{
                SS_DEBUG(" - sequence done, lifting min index for new vertex to remove"<<std::endl);
            }


            if(!vertex_to_collapse.is_valid()){
                SS_DEBUG(" --> next vertex in sequence is not removable, looking for a new one");
                vertex_to_collapse = trimmed_copy_.find_removable_vertex_with_cell_valence_lower_than_3(vertices_to_collapse,
                                                                                                        collapsed_prop_,
                                                                                                        min_index);
                if(vertex_to_collapse.idx() == -1){
                    SS_DEBUG(" ERROR - couldn't find next vertex to collapse among "<<vertices_to_collapse<<std::endl);
                    SS_DEBUG(" trimmed copy details: "); trimmed_copy_.print_details();
                    return -1;
                }

                IF_SS_DEBUG(c_val = trimmed_copy_.cell_valence(vertex_to_collapse);)
                SS_DEBUG(" - min-valence vertex is "<<vertex_to_collapse<<" with cell-valence "<<c_val);
            }

            //vertex_to_collapse = trimmed_copy_.find_removable_vertex_with_cell_valence_lower_than_3(vertices_to_collapse,
            //                                                                                        collapsed_prop_);

            if(!vertex_to_collapse.is_valid()){
                SS_DEBUG(" ERROR - couldn't find next vertex to collapse among "<<vertices_to_collapse<<std::endl);
                SS_DEBUG(" trimmed copy details: "); trimmed_copy_.print_details();
                return -1;
            }
        }else{

            vertex_to_collapse = trimmed_copy_.find_removable_vertex_with_cell_valence_lower_than_3(vertices_to_collapse,
                                                                                                    collapsed_prop_);

        }
        c_val = trimmed_copy_.cell_valence(vertex_to_collapse);



        if(c_val != 1 && c_val != 2){
            SS_DEBUG(" - couldn't find a vertex of valence < 3 to remove, unlocking situation and restarting RS sequence "<<full_reverse_shelling_sequence_<<std::endl<<" after "<<removed_vertices_count<<" deletions");
            max_vertex_index_at_last_RS_iterations_ = n_vertices();
            //return SS_RESET_REQUESTED;

            VertexHandle mid_vertex;
            //NOTE: removing this for now, because I'm not sure if this would make the sequence re-doing more complicated
            //auto split_result = split_non_blocking_base_boundary_edge();
            int split_result = 1;
            if(split_result){

                //SS_DEBUG(" -- couldn't split an edge to unlock the situation after performing RS sequence "<<full_reverse_shelling_sequence_<<". Resetting the reverse-shelling sequence after "<<removed_vertices_count<<" deletions ");

                //SS_DEBUG(" -- looking for any base boundary edge to split and back-propagate");


                auto to_split = find_trim_base_boundary_edge_to_start_split_path();

                if(to_split.idx() ==-1){
                    SS_DEBUG(" ERROR - couldn't find base boundary edge to split"<<std::endl);
                    return SS_ERROR;
                }

                SS_DEBUG(" -- found edge to split "<<edge(to_split)<<", backtracking to get split path");
                //IF_SS_DEBUG(trimmed_copy_.print_details());


                int remaining_vertices_in_trimmed_copy = trimmed_copy_.n_logical_vertices();

                if(remaining_vertices_in_trimmed_copy > remaining_trimmed_copy_vertices_at_last_RS_iteration_){
                    SS_DEBUG(" ERROR - there were "<<remaining_trimmed_copy_vertices_at_last_RS_iteration_<<" remaining vertices in the trimmed copy but now there are "<<remaining_vertices_in_trimmed_copy<<std::endl);
                    return SS_ERROR;
                }
                remaining_trimmed_copy_vertices_at_last_RS_iteration_ = remaining_vertices_in_trimmed_copy;

                std::deque<EdgeHandle> path_to_split;
                auto sp_result(-1);

                if(backtrack_splits_along_shortest_path){
                    sp_result = find_shortest_dual_path_to_boundary(to_split,
                                                                    path_to_split);

                }else{
                    sp_result = find_split_path_along_RS_sequence(to_split,
                                                                   collapsed_edge_prop,
                                                                   path_to_split);
                }


                if(sp_result){
                    SS_DEBUG(" ERROR - couldn't find path to split starting from edge "<<edge(to_split)<<
                               "using "<<(backtrack_splits_along_shortest_path ? "shortest-path " : "safe ")<<"method"<<std::endl);
                    return SS_ERROR;
                }

                SS_DEBUG(" -- applying "<<path_to_split.size()<<" splits path");
                for(auto to_split: path_to_split){
                    split_base_edge_and_register_split(to_split);
                    SS_DEBUG(" -- performed split "<<split_list_.back());
                    SS_DEBUG(" -- added vertex "<<vertices_to_collapse.back()<<" to vertices to collapse");

                    if(original_edge_prop_[to_split]){
                        vertices_to_collapse.push_back(split_list_.back().cone_new_vertex);
                    }else{

                        HalfEdgeHandle to_collapse(-1);
                        auto mid_vertex = split_list_.back().cone_new_vertex;
                        SS_DEBUG(" -- split edge is not an original one, collapsing the mid-vertex");
                        //SS_DEBUG(" -- split edge "<<edge(to_split)<<" is not an original one, collapsing the mid-vertex"<<std::endl);;

                        base_collapse_count++;
                        for(auto out_he: outgoing_halfedges(mid_vertex)){
                            auto neighbor = to_vertex_handle(out_he);
                            if(neighbor.idx() + 1 == mid_vertex.idx()){
                                SS_DEBUG(" ---> neighbor "<<neighbor<<" is the previous mid-vertex, collapsing edge "<<halfedge(out_he));
                                //SS_DEBUG(" ---> neighbor "<<neighbor<<" is the previous mid-vertex, collapsing edge "<<halfedge(out_he)<<std::endl);

                                to_collapse = out_he;
                            }
                        }

                        if(!to_collapse.is_valid()){
                            SS_DEBUG(" ERROR - couldn't find an outoing halfedge from "<<mid_vertex<<" to collapse"<<std::endl);
                            return SS_ERROR;
                        }

                        auto collapsed_to_vertex = collapse_edge(to_collapse);
                        split_list_.back().cone_collapsed_to_vertex = collapsed_to_vertex;

                    }

                }

                SS_DEBUG(" -- last split was "<<split_list_.back());
                return SS_RESET_REQUESTED;
                //SS_DEBUG(" ERROR - couldn't split an edge to unlock the situation"<<std::endl);
                //return SS_EDGE_SPLIT_ERROR;
            }else{
                vertices_to_collapse.push_back(split_list_.back().cone_new_vertex);
                //SS_DEBUG(" - cell valence of new vertex: "<<trimmed_copy_.cell_valence(split_list_.back().cone_new_vertex));
            }

            //SS_DEBUG("--> updated vertices to collapse to "<<vertices_to_collapse);
        }

        if(c_val != 1 && c_val != 2){
            SS_DEBUG(" ERROR - couldn't find a valence-1 or valence-2 vertex to collapse"<<std::endl);

            SS_DEBUG(" trimmed cone details:"); trimmed_copy_.print_details();
            return SS_ERROR;
        }

        cell_valence_at_deletion_time_prop[vertex_to_collapse] = c_val;
        //SS_DEBUG(" - set deletion cell-valence of "<<vertex_to_collapse<<" at "<<c_val);


        //SS_DEBUG(" - removed vertex "<<vertex_to_collapse<<" from trimmed cone after "<<removed_vertices_count<<" deletions");
        full_reverse_shelling_sequence_.push_back(vertex_to_collapse);

        if(!backtrack_splits_along_shortest_path){
            SS_DEBUG(" - updating collapse edge prop with vertex to collapse "<<vertex_to_collapse<<std::endl);
            auto update_result = update_collapsed_edge_prop(vertex_to_collapse,
                                                            c_val,
                                                            collapsed_edge_prop,
                                                            collapse_depth_prop);

            if(update_result){
                SS_DEBUG(" ERROR - couldn't update target edge prop"<<std::endl);
                return SS_ERROR;
            }
        }

        //SS_DEBUG("---------------");

        //then delete the old base vertex from the trimmed copy
        trimmed_copy_.delete_vertex(vertex_to_collapse);


        removed_vertices_count++;

    }


    if(base_collapse_count){
        SS_DEBUG(" - edges collapsed during RSS: "<<base_collapse_count<<std::endl);
    }
    SS_DEBUG(" --------------------------------------------------");
    SS_DEBUG(" ------ DONE WITH REVERSE SHELLING SEQUENCE -------");
    SS_DEBUG(" --------------------------------------------------");

    return 0;
}


#if 0
int StarShapifyableExpansionCone::run_direct_reverse_shelling_sequence(std::vector<VertexHandle>& vertices_to_collapse){


    SS_DEBUG(" -----------------------------------------------------------------------------------");
    SS_DEBUG(" ----------------------- RUNNING DIRECT REVERSE-SHELLING SEQUENCE -------------------------");
    SS_DEBUG(" -----------------------------------------------------------------------------------");

    //used as hard limit for valence reduction
    //const int initial_n_edges(n_edges());

    auto cell_valence_at_deletion_time_prop = trimmed_copy_.request_vertex_property<int>("c_val_at_deletion");
    trimmed_copy_.set_persistent(cell_valence_at_deletion_time_prop);


    //maps an edge to another edge that was collapsed onto it.
    auto collapsed_edge_prop = request_edge_property<EdgeHandle>("", EdgeHandle(-1));
    auto collapse_depth_prop = request_edge_property<int>();

    int removed_vertices_count(0);
    while(removed_vertices_count < (int)vertices_to_collapse.size()){

        //SS_DEBUG(" ================= removing "<<(removed_vertices_count+1)<<"-th vertex...");

        VertexHandle vertex_to_collapse(-1);
        int c_val(0);

        vertex_to_collapse = trimmed_copy_.find_removable_vertex_with_cell_valence_lower_than_3(vertices_to_collapse,
                                                                                                collapsed_prop_);
        c_val = trimmed_copy_.cell_valence(vertex_to_collapse);



        if(c_val != 1 && c_val != 2){
            SS_DEBUG(" - couldn't find a vertex of valence < 3 to remove, unlocking situation and restarting RS sequence "<<full_reverse_shelling_sequence_<<std::endl<<" after "<<removed_vertices_count<<" deletions");
            //return SS_RESET_REQUESTED;

            VertexHandle mid_vertex;
            //NOTE: removing this for now, because I'm not sure if this would make the sequence re-doing more complicated
            //auto split_result = split_non_blocking_base_boundary_edge();
            int split_result = 1;
            if(split_result){

                //SS_DEBUG(" -- couldn't split an edge to unlock the situation after performing RS sequence "<<full_reverse_shelling_sequence_<<". Resetting the reverse-shelling sequence after "<<removed_vertices_count<<" deletions ");

                //SS_DEBUG(" -- looking for any base boundary edge to split and back-propagate");


                auto to_split = find_trim_base_boundary_edge_to_start_split_path();

                if(to_split.idx() ==-1){
                    SS_DEBUG(" ERROR - couldn't find base boundary edge to split"<<std::endl);
                    return SS_ERROR;
                }

                //SS_DEBUG(" -- found edge to split "<<edge(to_split)<<", backtracking to get split path");


                std::deque<EdgeHandle> path_to_split;
                auto sp_result(-1);

                sp_result = find_split_path_along_RS_sequence(to_split,
                                                              collapsed_edge_prop,
                                                              path_to_split);


                if(sp_result){
                    SS_DEBUG(" ERROR - couldn't find path to split starting from edge "<<edge(to_split)<<std::endl);
                    return SS_ERROR;
                }

                //SS_DEBUG(" -- applying "<<path_to_split.size()<<" splits path");
                for(auto to_split: path_to_split){
                    split_base_edge_and_register_split(to_split);
                    vertices_to_collapse.push_back(split_list_.back().cone_new_vertex);
                    SS_DEBUG(" -- performed split "<<split_list_.back());
                    //SS_DEBUG(" -- added vertex "<<vertices_to_collapse.back()<<" to vertices to collapse");
                }
                //SS_DEBUG(" -- last split was "<<split_list_.back());

                vertex_to_collapse = trimmed_copy_.find_removable_vertex_with_cell_valence_lower_than_3(vertices_to_collapse,
                                                                                                        collapsed_prop_);

                if(vertex_to_collapse.is_valid()){
                    SS_DEBUG(" ERROR - couldn't immediately collapse new vertex"<<std::endl);
                    return SS_ERROR;
                }

                c_val = trimmed_copy_.cell_valence(vertex_to_collapse);


                return SS_RESET_REQUESTED;
                //SS_DEBUG(" ERROR - couldn't split an edge to unlock the situation"<<std::endl);
                //return SS_EDGE_SPLIT_ERROR;
            }else{
                vertices_to_collapse.push_back(split_list_.back().cone_new_vertex);
                //SS_DEBUG(" - cell valence of new vertex: "<<trimmed_copy_.cell_valence(split_list_.back().cone_new_vertex));
            }

            //SS_DEBUG("--> updated vertices to collapse to "<<vertices_to_collapse);
        }

        if(c_val != 1 && c_val != 2){
            SS_DEBUG(" ERROR - couldn't find a valence-1 or valence-2 vertex to collapse"<<std::endl);

            SS_DEBUG(" trimmed cone details:"; trimmed_copy_.print_details();
            return SS_ERROR;
        }

        cell_valence_at_deletion_time_prop[vertex_to_collapse] = c_val;
        //SS_DEBUG(" - set deletion cell-valence of "<<vertex_to_collapse<<" at "<<c_val);


        SS_DEBUG(" - removed vertex "<<vertex_to_collapse<<" from trimmed cone after "<<removed_vertices_count<<" deletions");

        auto update_result = update_collapsed_edge_prop(vertex_to_collapse,
                                                        c_val,
                                                        collapsed_edge_prop,
                                                        collapse_depth_prop);

        if(update_result){
            SS_DEBUG(" ERROR - couldn't update target edge prop"<<std::endl);
            return SS_ERROR;
        }


       // SS_DEBUG("---------------");

        //then delete the old base vertex from the trimmed copy
        trimmed_copy_.delete_vertex(vertex_to_collapse);


        removed_vertices_count++;

    }


    SS_DEBUG(" --------------------------------------------------");
    SS_DEBUG(" ------ DONE WITH REVERSE SHELLING SEQUENCE -------");
    SS_DEBUG(" --------------------------------------------------");

    return 0;

}
#endif




int StarShapifyableExpansionCone::update_collapsed_edge_prop(const VertexHandle& vertex_to_remove,
                                                             const int cell_valence,
                                                             EdgePropertyT<EdgeHandle>& collapsed_edge_prop,
                                                             EdgePropertyT<int>& collapse_depth_prop) const{

    if(cell_valence == 1){

        SS_DEBUG(" --> removed vertex "<<vertex_to_remove<<" with cell valence = 1, looking for source and target edge");

        //because there should only be one cell incident to vh
        auto c = *trimmed_copy_.vc_iter(vertex_to_remove);

        //we only use one of the source edges for now
        EdgeHandle target_e(-1);
        EdgeHandle source_e1(-1);
        EdgeHandle source_e2(-1);

        for(auto ce_it = trimmed_copy_.ce_iter(c); ce_it.valid(); ce_it++){
            if(trimmed_copy_.is_base_boundary_edge(*ce_it)){
                if(source_e1.is_valid()){
                    source_e2 = *ce_it;
                }else{
                    source_e1 = *ce_it;
                }
            }else if(trimmed_copy_.is_base_edge(*ce_it)){
                target_e = *ce_it;
            }
        }

        if(!target_e.is_valid() ||
                !source_e1.is_valid() ||
                !source_e2.is_valid()){
            SS_DEBUG(" ERROR - couldn't find all the edges for vertex to collapse "<<vertex_to_remove<<std::endl);
            SS_DEBUG(" target edge = "<<target_e<<std::endl);
            SS_DEBUG(" source edge 1 = "<<source_e1<<std::endl);
            SS_DEBUG(" source edge 2 = "<<source_e2<<std::endl);
            return -1;
        }

        EdgeHandle source_e = source_e1;
        if(collapse_depth_prop[source_e1] > collapse_depth_prop[source_e2]){
            source_e = source_e2;
            //SS_DEBUG(" -- source edge 2 "<<edge(source_e2)<<" has lower depth ("<<collapse_depth_prop[source_e2]<<") than source edge 1 "<<edge(source_e1)<<" (depth"<<collapse_depth_prop[source_e1]<<"), using this one");
        }else{
            //SS_DEBUG(" -- source edge 1 "<<edge(source_e1)<<" has lower depth ("<<collapse_depth_prop[source_e1]<<") than source edge 2 ("<<collapse_depth_prop[source_e2]<<"), using this one");
        }

        collapsed_edge_prop[target_e] = source_e;
        collapse_depth_prop[target_e] = collapse_depth_prop[source_e] + 1;
        //SS_DEBUG(" --- set edge "<<edge(target_e)<<" as target with depth "<<collapse_depth_prop[target_e]<<" for egde "<<edge(source_e));

    }else if(cell_valence == 2){

        auto target_v = find_RS_target_vertex(vertex_to_remove);

        SS_DEBUG(" --> removed vertex "<<vertex_to_remove<<" with cell valence = 2, target vertex is "<<target_v);

        if(!target_v.is_valid()){
            SS_DEBUG(" ERROR - couldn't find target vertex for vertex to collapse "<<vertex_to_remove<<std::endl);
            return -1;
        }

        int added_edges(0);
        for(auto vc_it = trimmed_copy_.vc_iter(vertex_to_remove); vc_it.valid(); vc_it++){
            auto c_vertices = get_cell_vertices(*vc_it);
            for(auto v: c_vertices){
                if(v != vertex_to_remove &&
                        v != target_v &&
                        !tip_vertices_prop_[v]){

                    auto source_e = edge_handle(trimmed_copy_.find_halfedge(vertex_to_remove, v));
                    if(!source_e.is_valid()){
                        SS_DEBUG(" ERROR - couldn't find source edge ("<<vertex_to_remove<<", "<<v<<")"<<std::endl);
                        return -1;
                    }

                    auto target_e = edge_handle(trimmed_copy_.find_halfedge(target_v, v));
                    if(!target_e.is_valid()){
                        SS_DEBUG(" ERROR - couldn't find target edge ("<<target_v<<", "<<v<<")"<<std::endl);
                        return -1;
                    }

                    collapsed_edge_prop[target_e] = source_e;
                    collapse_depth_prop[target_e] = collapse_depth_prop[source_e] + 1;
                    //SS_DEBUG(" --- set edge "<<edge(target_e)<<" as target with depth "<<collapse_depth_prop[target_e]<<" for egde "<<edge(source_e));
                    added_edges++;
                }
            }
        }
        if(added_edges != 2){
            SS_DEBUG(" ERROR - found "<<added_edges<<" collapsed edges instead of 2"<<std::endl);
            return -1;
        }
    }

    return 0;
}

EdgeHandle StarShapifyableExpansionCone::find_trim_base_boundary_edge_to_start_split_path() const{

    EdgeHandle valid_edge(-1);

    int min_incident_vertex_cell_valence(std::numeric_limits<int>::max());

    for(auto eh: trimmed_copy_.edges()){

        if(trimmed_copy_.is_base_boundary_edge(eh) &&
                is_not_part_of_witness_vertex_1_ring(eh) /*&&
                            (witness_vertex_1_ring_base_vertex_prop_[edge(eh).from_vertex()] ||
                             witness_vertex_1_ring_base_vertex_prop_[edge(eh).to_vertex()])*/
                /* && original_edge_prop_[eh]*/){

            //SS_DEBUG(" --- trying trimmed copy base boundary, non witness 1-ring edge "<<edge(eh));

            auto base_op_vertex = opposite_trim_base_vertex(eh);

            if(!base_op_vertex.is_valid()){
                SS_DEBUG(" ERROR - couldn't recover trim base op vertex to boundary edge "<<edge(eh)<<std::endl);
                return EdgeHandle(-1);
            }

            if(trimmed_copy_.is_base_boundary_vertex(base_op_vertex)){
                //SS_DEBUG(" --- but base op vertex is boundary, skipping");
                continue;
            }


            if(this->is_base_boundary_edge(eh)){
                //SS_DEBUG(" -> found original cone base boundary edge, using this one");
                valid_edge = eh;
                break;
            }

            int incident_vertex_cell_valence(0);
            auto from_v = trimmed_copy_.edge(eh).from_vertex();
            auto to_v   = trimmed_copy_.edge(eh).to_vertex();
            if(witness_vertex_1_ring_base_vertex_prop_[from_v]){
                incident_vertex_cell_valence = trimmed_copy_.cell_valence(to_v);
                //SS_DEBUG(" --> using to-vertex "<<to_v<<" cell valence = "<<incident_vertex_cell_valence<<std::endl);
            }else if(witness_vertex_1_ring_base_vertex_prop_[to_v]){
                incident_vertex_cell_valence = trimmed_copy_.cell_valence(from_v);
                //SS_DEBUG(" --> using from-vertex "<<from_v<<" cell valence = "<<incident_vertex_cell_valence<<std::endl);
            }else{
                incident_vertex_cell_valence = std::min(trimmed_copy_.cell_valence(to_v),
                                                        trimmed_copy_.cell_valence(from_v));
                /*SS_DEBUG(" --> from-vertex "<<from_v<<" cell valence = "<<trimmed_copy_.cell_valence(from_v)<<std::endl);
                SS_DEBUG(" --> to-vertex "<<to_v<<" cell valence = "<<trimmed_copy_.cell_valence(to_v)<<std::endl);
                SS_DEBUG(" ---> using "<<incident_vertex_cell_valence<<std::endl);*/


            }

            if(incident_vertex_cell_valence < min_incident_vertex_cell_valence){
                //SS_DEBUG(" --- updated base boundary edge to "<<edge(eh)<<" with incident vertex with cell-valence = "<<incident_vertex_cell_valence);

                valid_edge = eh;
                min_incident_vertex_cell_valence = incident_vertex_cell_valence;

                if(min_incident_vertex_cell_valence < 3){
                    SS_DEBUG(" ERROR - this shouldn't happen. If there's a vertex with cell-valence < 3 it should be removable"<<std::endl);
                    SS_DEBUG(" cells around to-vertex "<<to_v<<": "<<std::endl);
                    for(auto vc_it = trimmed_copy_.vc_iter(to_v); vc_it.valid(); vc_it++){
                        SS_DEBUG(" - "<<get_cell_vertices(*vc_it)<<std::endl);
                    }
                    //SS_DEBUG(" temporarily trying this out"<<std::endl);
                    return EdgeHandle(-1);
                }else if(min_incident_vertex_cell_valence == std::numeric_limits<int>::max()){
                    SS_DEBUG(" ERROR - couldn't find a boundary edge to split, which shouldn't happen"<<std::endl);
                    return EdgeHandle(-1);
                }else if(min_incident_vertex_cell_valence == 3){
                    SS_DEBUG(" ---> min cell valence == 3, using this one");
                    break;
                }
            }
        }
    }

    return valid_edge;
}



int StarShapifyableExpansionCone::split_non_blocking_base_boundary_edge(){

    SS_DEBUG(" --------------------------------------------------------------------------");
    SS_DEBUG(" -------------- splitting non-blocking boundary edge...");

    auto cell_valence_at_deletion_time_prop = trimmed_copy_.request_vertex_property<int>("c_val_at_deletion");

    EdgeHandle to_split(-1);
    VertexHandle from_vertex(-1), to_vertex(-1);

    for(auto eh: trimmed_copy_.edges()){
        from_vertex = edge(eh).from_vertex();
        to_vertex   = edge(eh).to_vertex();

        if(trimmed_copy_.valence(eh) == 2 &&
                !tip_vertices_prop_[from_vertex] &&
                !tip_vertices_prop_[to_vertex] &&
                from_vertex != witness_vertex_ &&
                to_vertex   != witness_vertex_){
            //SS_DEBUG(" -- found edge "<<edge(eh)<<" with valence 2");

            bool found_deleted_op_vertex_with_c_val_greater_than_1(false);
            bool found_witness_vertex_as_op_vertex(false);

            for(auto ef_it = ef_iter(eh); ef_it.valid(); ef_it++){

                auto f_vertices = get_halfface_vertices(halfface_handle(*ef_it, 0));
                //SS_DEBUG(" --- checking face "<<f_vertices);
                VertexHandle op_vertex(-1);
                for(auto fv: f_vertices){
                    if(fv != from_vertex &&
                            fv != to_vertex){
                        op_vertex = fv;
                    }
                }
                if(op_vertex.idx() == -1){
                    SS_DEBUG(" ERROR - couldn't find op vertex to edge "<<edge(eh)<<
                               " on face "<<get_halfface_vertices(halfface_handle(*ef_it, 0))<<std::endl);
                    return SS_ERROR;
                }

                if(cell_valence_at_deletion_time_prop[op_vertex] > 1){
                    found_deleted_op_vertex_with_c_val_greater_than_1 = true;
                    //SS_DEBUG(" --> op vertex "<<op_vertex<<" to edge "<<edge(eh)<<" had valence > 1 when deleted");
                }

                if(op_vertex == witness_vertex_){
                    //SS_DEBUG(" --> op vertex "<<op_vertex<<" to edge "<<edge(eh)<<" is the witness vertex");
                    found_witness_vertex_as_op_vertex = true;
                }
            }

            if(!found_deleted_op_vertex_with_c_val_greater_than_1 &&
                    !found_witness_vertex_as_op_vertex){
                to_split = eh;
                //SS_DEBUG(" ----> edge "<<edge(eh)<<" is usable");
                break;
            }
        }
    }

    if(to_split.idx() == -1){
        SS_DEBUG(" --> found no valence-2 boundary edge in trimmed copy. Trimmed copy details: ");
        //IF_SS_DEBUG(trimmed_copy_.print_details(2));
        return 1;
    }


    //split both the trimmed copy and the original
    auto mid_vertex = trimmed_copy_.split_edge(to_split);
    auto original_cone_mid_vertex = split_base_edge_and_register_split(to_split);


    trimmed_copy_.set_vertex(mid_vertex, vertex(original_cone_mid_vertex));

    SS_DEBUG(" -------------- done!");
    SS_DEBUG(" --------------------------------------------------------------------------");

    return 0;

}




int StarShapifyableExpansionCone::run_collapse_sequence(const std::vector<VertexHandle>& vertices_to_collapse){

    SS_DEBUG(" ------------------------------------------"<<std::endl);
    SS_DEBUG(" ------- RUNNING COLLAPSE SEQUENCE --------"<<std::endl);
    SS_DEBUG(" ------------------------------------------"<<std::endl);

    SS_DEBUG(" full reverse-shelling sequence: "<<full_reverse_shelling_sequence_);
    //SS_DEBUG(" full reverse-shelling sequence: "<<full_reverse_shelling_sequence_<<std::endl);

    //while(collapsed_vertices_count < (int)vertices_to_collapse.size()){
    int i(0);
    for(auto vertex_to_collapse: full_reverse_shelling_sequence_){

        check_for_timeout();
        //SS_DEBUG(" ================= collapsing "<<(i+1)<<"-th vertex "<<vertex_to_collapse);

        int c_val = trimmed_copy_.cell_valence(vertex_to_collapse);

        //SS_DEBUG(" - min-valence vertex is "<<vertex_to_collapse<<" with cell-valence "<<c_val);


        switch(c_val){
        case 1: {

            auto target_face = opposite_face(vertex_to_collapse);

            SS_DEBUG(" -> collapsing vertex "<<vertex_to_collapse<<" to face "<<get_halfface_vertices(halfface_handle(target_face, 0)));

            auto v_prime = collapse_vertex(vertex_to_collapse,
                                           get_position_from_face(target_face));

            if(!v_prime.is_valid()){
                SS_DEBUG(" error while collapsing vertex"<<std::endl);
                return -1;
            }

            //attached_vertices_f_prop_[target_face] = {v_prime, collapsed_vertices_count};
            collapse_list_.push_back({v_prime, get_halfface_vertices(halfface_handle(target_face, 0))});
            //SS_DEBUG(" -> attached vertex "<<v_prime<<" to face "<<get_halfface_vertices(halfface_handle(target_face, 0)));

            break;
        }
        case 2: {

            auto target_edge = find_target_edge(vertex_to_collapse);
            if(target_edge.idx() == -1){
                SS_DEBUG( "error while trying to find target edge"<<std::endl);
                return -1;
            }
            SS_DEBUG(" -> collapsing vertex "<<vertex_to_collapse<<" to edge "<<edge(target_edge));

            auto v_prime = collapse_vertex(vertex_to_collapse,
                                           get_position_from_edge(target_edge));

            if(!v_prime.is_valid()){
                SS_DEBUG(" error while collapsing vertex"<<std::endl);
                return -1;
            }

            collapse_list_.push_back({v_prime, {edge(target_edge).from_vertex(), edge(target_edge).to_vertex()}});
            //SS_DEBUG(" -> attached vertex "<<v_prime<<" to edge "<<edge(target_edge));

            break;
        }
        default:{
            SS_DEBUG(" ERROR - next vertex in reverse-shelling sequence is cell-valence "<<c_val<<std::endl);
            SS_DEBUG(" trimmed cone details:"); trimmed_copy_.print_details();
            return 1;
        }
        }

#if ENABLE_BAD_TETS_CHECKS_IN_LOOPS
        auto bad_tets = ExactBadTetFinder::findBadTets(*this, vertex_position_prop_);
        if(bad_tets.second.size()){
            SS_DEBUG(" ERROR - cone contains flipped tets after collapse:"<<std::endl);
            for(auto bad_tet: bad_tets.second){
                SS_DEBUG(" - "<<bad_tet<<" :" <<get_cell_vertices(bad_tet)<<std::endl);
            }
            return -1;
        }
#endif

        //SS_DEBUG(" ================= done! Collapsed vertex "<<vertex_to_collapse);

        //then delete the old base vertex from the trimmed copy
        trimmed_copy_.delete_vertex(vertex_to_collapse);
        //SS_DEBUG(" - removed vertex "<<vertex_to_collapse<<" from trimmed cone");

        i++;
    }

    SS_DEBUG(" ------------------------------------------"<<std::endl);
    SS_DEBUG(" ------ DONE WITH COLLAPSE SEQUENCE -------"<<std::endl);
    SS_DEBUG(" ------------------------------------------"<<std::endl);

    return 0;
}






int StarShapifyableExpansionCone::run_contraction_sequence(bool reduce_precision){

    SS_DEBUG(" ------------------------------------------"<<std::endl);
    SS_DEBUG(" ------ RUNNING CONTRACTION SEQUENCE ------"<<std::endl);
    SS_DEBUG(" ------------------------------------------"<<std::endl);

    if(!vertices_to_contract_.empty()){
        SS_DEBUG(" ERROR - vertices to contract list is not empty"<<std::endl);
        return -1;
    }

    auto contraction_start = std::chrono::high_resolution_clock::now();

    //temp check
    /*ExpansionCone updated_cone;
    set_up_1_ring_neighborhood_as_expansion_cone(*this,
                                                 vertex_position_prop_,
                                                 /cone_tip_vertices().begin(),
                                                 witness_vertex_,
                                                 updated_cone);
    VertexPosition pos;
    auto geo_exp_result = updated_cone.is_geo_expandable(pos);

    if(geo_exp_result){
        SS_DEBUG(" ERROR - temp check for geo-exp failed before contraction"<<std::endl);
        SS_DEBUG(" updated cone: "<<std::endl);
        updated_cone.print_details();
        SS_DEBUG(" cone: "<<std::endl);
        print_details();
        return -1;
    }*/

    SS_DEBUG(" ==================================== 1-ring contraction "<<std::endl);
    //1. contract the 1-ring neighborhood and other tip vertices
    auto core_contraction_result = contract_core();
    if(core_contraction_result){
        SS_DEBUG(" ERROR - failed to contract core"<<std::endl);
        return -1;
    }


//#if ENABLE_BAD_TETS_CHECKS_OUT_OF_LOOPS
    IF_SS_DEBUG(auto initial_flipped_start = std::chrono::high_resolution_clock::now(););
    if(ExactBadTetFinder::meshContainsFlippedTets(*this, vertex_position_prop_)){
        SS_DEBUG(" ERROR - cone contains flipped tets after contracting 1-ring neighborhood: "<<std::endl);
        auto flipped_tets = ExactBadTetFinder::findBadTets(*this, vertex_position_prop_).second;
        for(auto flipped_tet: flipped_tets){
            SS_DEBUG(" - "<<flipped_tet<<": "<<get_cell_vertices(flipped_tet)<<std::endl);
        }
        return -1;
    }
//#endif

#if ENABLE_TIMINGS
    auto initial_flipped_end = std::chrono::high_resolution_clock::now();
    float initial_flipped_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(initial_flipped_end - initial_flipped_start).count() / 1000000;
    SS_DEBUG(" - initial flipped test check duration: "<<initial_flipped_duration_s<<std::endl);
#endif

    for(auto v: vertices()){
        contracted_prop_[v] = true;
        //SS_DEBUG(" - set "<<v<<" as contracted"<<std::endl);
    }

    //this allows to skip vertices that had nothing collapsed to
    //i.e. those are vertices that are at the "end" of a collapse sequence and don't need to be contracted
    auto collapsed_to_prop = request_vertex_property<std::vector<VertexHandle>>();
    for(const auto& collapse: collapse_list_){
        for(auto v: collapse.second){
            collapsed_to_prop[v].push_back(collapse.first);
        }
    }


    SS_DEBUG(" ==================================== full contraction "<<std::endl);

    for(auto& collapse: collapse_list_){
        vertices_to_contract_.push_back(collapse.first);
        contracted_prop_[collapse.first] = false;
        //SS_DEBUG(" - set "<<collapse.first<<" as uncontracted"<<std::endl);
    }

    //SS_DEBUG(" - vertices to contract (reverse order): "<<vertices_to_contract_<<std::endl);

    for(int i(collapse_list_.size()-1); i>= 0; i--){
        auto& collapse = collapse_list_[i];
        SS_DEBUG(" -----------------------------------------------------------------------------");
        SS_DEBUG(" ------- handling collapse "<<i<<": "<<collapse.first<<" -> "<<collapse.second);


        check_for_timeout();

        if(vertices_to_contract_.back() != collapse.first){
            SS_DEBUG(" ERROR - back of vertices to contract vector is "<<vertices_to_contract_.back()<<" and last contracted vertex is "<<collapse.first<<std::endl);
            return -1;
        }
        vertices_to_contract_.pop_back();
        contracted_prop_[collapse.first] = true;
        //SS_DEBUG(" - set "<<collapse.first<<" as contracted"<<std::endl);

        if(collapsed_to_prop[collapse.first].empty()){
            //SS_DEBUG(" --> nothing was collapsed to vertex "<<collapse.first<<", skipping");
            continue;
        }

        //NOTE: moving this here because we don't care about vertices that don't need to be contracted
#if ENABLE_TIMINGS
        SS_DEBUG(" -----------------------------------------------------------------------------"<<std::endl);
        SS_DEBUG(" ------- handling collapse "<<i<<": "<<collapse.first<<" -> "<<collapse.second<<std::endl);
        SS_DEBUG(" --> remaining vertices to contract: "<<(vertices_to_contract_.size()+1)<<"/"<<collapse_list_.size()<<std::endl);
        SS_DEBUG(" --> max position size = "<<new_position_max_precision_<<std::endl);
#endif

#warning TODO: put this in function
        SS_DEBUG(" -- collapse target size = "<<collapse.second.size());
        bool needs_to_be_updated(false);
        for(auto& target_v: collapse.second){
            if(secondary_vertex_prop_[target_v].idx() != -1){
                needs_to_be_updated = true;
                target_v = secondary_vertex_prop_[target_v];
            }
        }

        if(needs_to_be_updated){
            SS_DEBUG(" ---> one of the target vertices has a secondary vertex, updating to "<<collapse.second);

            IF_SS_DEBUG(auto contr_pos_start = std::chrono::high_resolution_clock::now(););
            auto contraction_pos = get_contraction_position_from_vertices(collapse.first,
                                                                          collapse.second,
                                                                          collapsed_to_prop);

            //temp check
#if 0
            auto v_dist = CGAL::squared_distance(OVMvec3ToCGALPoint3(contraction_pos), OVMvec3ToCGALPoint3(vertex(*cone_tip_vertices().begin())));
            //SS_DEBUG(" --- contraction position to tip distance = "<<CGAL::to_double(v_dist)<<std::endl);
            ExactType target_min_distance(std::numeric_limits<double>::max());
            for(auto target_v: collapse.second){
                if(is_cone_tip(target_v)){
                    continue;
                }
                auto dist = CGAL::squared_distance(OVMvec3ToCGALPoint3(vertex(target_v)), OVMvec3ToCGALPoint3(vertex(*cone_tip_vertices().begin())));
                //SS_DEBUG(" --- distance to tip for target vertex "<<target_v<<" = "<<CGAL::to_double(dist)<<std::endl);
                /*if(dist < v_dist){
                    SS_DEBUG(" ERROR - contraction position is at "<<CGAL::to_double(v_dist)<<" >= target distance "<<CGAL::to_double(dist)<<std::endl);
                    return -1;
                }*/
            }
#endif




#if ENABLE_TIMINGS
            auto contr_pos_end = std::chrono::high_resolution_clock::now();
            float contr_pos_start_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(contr_pos_end - contr_pos_start).count() / 1000000;
            SS_DEBUG(" - contraction position computation duration: "<<contr_pos_start_duration_s<<std::endl);
#endif


            IF_SS_DEBUG(auto split_spoke_start = std::chrono::high_resolution_clock::now(););
            auto secondary_v_prime = split_spoke_and_move_vertex(collapse.first,
                                                                 contraction_pos,
                                                                 collapse.second,
                                                                 reduce_precision);

#if ENABLE_TIMINGS
            auto split_spoke_end = std::chrono::high_resolution_clock::now();
            float split_spoke_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(split_spoke_end - split_spoke_start).count() / 1000000;
            SS_DEBUG(" - split spoke duration: "<<split_spoke_duration_s<<std::endl);
#endif

            if(!secondary_v_prime.is_valid()){
                SS_DEBUG(" error while trying to split spoke"<<std::endl);
                return -1;
            }
            collapse.first = secondary_v_prime;
            SS_DEBUG(" ---> updated collapse to "<<collapse.first<<" -> "<<collapse.second<<", position size = "<<byte_size(vertex(collapse.first)));

#if ENABLE_BAD_TETS_CHECKS_IN_LOOPS
            if(ExactBadTetFinder::meshContainsFlippedTets(*this, vertex_position_prop_)){
                SS_DEBUG(" ERROR - cone contains flipped tets after handling collapse "<<collapse.first<<" -> "<<collapse.second<<std::endl);
                auto flipped_tets = ExactBadTetFinder::findBadTets(*this, vertex_position_prop_).second;
                for(auto flipped_tet: flipped_tets){
                    SS_DEBUG(" - "<<flipped_tet<<": "<<get_cell_vertices(flipped_tet)<<std::endl);
                }
                return -1;
            }
#endif

            contracted_prop_[collapse.first] = true;
            //SS_DEBUG(" - set "<<collapse.first<<" as contracted"<<std::endl);

#warning TODO: remove this once validated
        }else{
            SS_DEBUG(" ERROR - A vertex to contract has no updated target vertex"<<std::endl);
            return -1;
        }

#if ENABLE_BAD_TETS_CHECKS_IN_LOOPS
        auto one_ring_check_start = std::chrono::high_resolution_clock::now();
        ExpansionCone updated_cone;
        set_up_1_ring_neighborhood_as_expansion_cone(*this,
                                                     vertex_position_prop_,
                                                     *cone_tip_vertices().begin(),
                                                     updated_cone);
        VertexPosition pos;
        auto geo_exp_result = updated_cone.is_geo_expandable(pos);

        if(geo_exp_result){
            SS_DEBUG(" ERROR - temp check for geo-exp failed after contracting vertex "<<collapse.first<<std::endl);
            updated_cone.print_details();
            return -1;
        }

#if ENABLE_TIMINGS
        auto one_ring_check_end = std::chrono::high_resolution_clock::now();
        float one_ring_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(one_ring_check_end - one_ring_check_start).count() / 1000000;
        SS_DEBUG(" - 1-ring check duration: "<<one_ring_check_duration_s<<std::endl);
#endif
#endif

#if ENABLE_BAD_TETS_CHECKS_IN_LOOPS
        auto flipped_tets_check_start = std::chrono::high_resolution_clock::now();
        auto bad_tets = ExactBadTetFinder::findBadTets(*this, vertex_position_prop_);
        if(bad_tets.second.size()){
            SS_DEBUG(" ERROR - cone contains flipped tets after contraction:"<<std::endl);
            for(auto bad_tet: bad_tets.second){
                SS_DEBUG(" - "<<bad_tet<<" :" <<get_cell_vertices(bad_tet)<<std::endl);
            }
            return -1;
        }
#if ENABLE_TIMINGS
        auto flipped_tets_check_end = std::chrono::high_resolution_clock::now();
        float flipped_tets_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(flipped_tets_check_end - flipped_tets_check_start).count() / 1000000;
        SS_DEBUG(" - flipped tets check duration: "<<flipped_tets_check_duration_s<<std::endl);
#endif
#endif

        //SS_DEBUG(" -------");
        //temporary
        //tips_1_ring_neighborhood_is_visible_from_witness_vertex(true);

    }

    IF_SS_DEBUG(auto contraction_end = std::chrono::high_resolution_clock::now();)
    IF_SS_DEBUG(float contraction_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(contraction_end - contraction_start).count() / 1000000;);
    SS_DEBUG(" - total contraction duration: "<<contraction_duration_s<<std::endl);

    SS_DEBUG(" ------------------------------------------"<<std::endl);
    SS_DEBUG(" ----- DONE WITH CONTRACTION SEQUENCE -----"<<std::endl);
    SS_DEBUG(" ------------------------------------------"<<std::endl);

    return 0;
}



int StarShapifyableExpansionCone::contract_core(){
    std::vector<VertexHandle> one_ring_neighborhood;

    SS_DEBUG(" -------------- contracting core");



    //TODO: make this a function

    //auto valence_3_spoke_prop = request_vertex_property<bool>();
    for(auto out_he: outgoing_halfedges(witness_vertex_)){
        auto neighbor = to_vertex_handle(out_he);
        if(!tip_vertices_prop_[neighbor]){

            one_ring_neighborhood.push_back(neighbor);

            //TEMP check
            auto spoke_edge = edge_handle(base_vertex_to_spoke_edge_prop_[neighbor]);
            int val = valence(spoke_edge);
            int c_val = cell_valence(spoke_edge);

            if(val == 3 && val == c_val){
                //valence_3_spoke_prop[neighbor] = true;
                SS_DEBUG(" ERROR - found valence-3 spoke in core"<<std::endl);
                return -1;
            }
        }
    }


    //TBD if this is still useful
    auto collapsed_to_prop = request_vertex_property<bool>();
    for(const auto& collapse: collapse_list_){
        for(auto v: collapse.second){
            collapsed_to_prop[v] = true;
        }
    }

    SS_DEBUG(" ---------------------------------");
    IF_SS_DEBUG(const auto witness_tip_vertex = from_vertex_handle(base_vertex_to_spoke_edge_prop_[witness_vertex_]););
    SS_DEBUG(" - contracting core around edge "<<witness_tip_vertex<<"-"<<witness_vertex_);


    for(auto neighbor: one_ring_neighborhood){
        SS_DEBUG(" ---------------------------------");
        std::vector<VertexHandle> spoke_face_vertices = {neighbor,
                                                         witness_vertex_,
                                                         from_vertex_handle(base_vertex_to_spoke_edge_prop_[neighbor])};
        SS_DEBUG(" -- moving neighbor vertex "<<neighbor<<" to center of spoke face "<<spoke_face_vertices);

        /*if(!collapsed_to_prop[neighbor]){
            SS_DEBUG(" --> nothing was collapsed to neighbor "<<neighbor);
            if(valence(edge_handle(base_vertex_to_spoke_edge_prop_[neighbor])) == 3){
                SS_DEBUG(" but spoke edge has valence 3, still contracting");
            }else{
                SS_DEBUG(" and spoke edge valence != 3, skipping");
                continue;
            }
        }*/

        auto hfh = find_halfface(spoke_face_vertices);
        if(hfh.idx() == -1){
            SS_DEBUG(" ERROR - couldn't recover halfface handle"<<std::endl);
            return -1;
        }
        //SS_DEBUG(" -- halfface handle: "<<hfh<<std::endl);
        auto op_hfh = opposite_halfface_handle(hfh);
        if(op_hfh.idx() == -1){
            SS_DEBUG(" ERROR - couldn't recover opposite halfface handle"<<std::endl);
            return -1;
        }
        //SS_DEBUG(" -- opposite halfface handle: "<<op_hfh<<std::endl);


        //hf opposite to the boundary hf
        HalfFaceHandle boundary_hf_op_hfh(-1);

        if(is_boundary(hfh)){
            boundary_hf_op_hfh = op_hfh;
        }else if(is_boundary(op_hfh)){
            boundary_hf_op_hfh = hfh;
        }

        if(boundary_hf_op_hfh.idx() != -1){
            SS_DEBUG(" -- spoke face "<<spoke_face_vertices<<" is boundary, skipping");
            continue;
        }


        auto v_prime = split_spoke_and_move_vertex(neighbor,
                                                   get_position_from_vertices(spoke_face_vertices),
                                                   {witness_vertex_, *cone_tip_vertices().begin()});

        if(!v_prime.is_valid()){
            SS_DEBUG(" error while splitting spoke"<<std::endl);
            return -1;
        }
        //temporary
        //tips_1_ring_neighborhood_is_visible_from_witness_vertex(true);
    }



    //special stuff for clusters starts here
    if(cone_tip_vertices_.size() > 1){
        SS_DEBUG(" ERROR - cluster stuff not ready"<<std::endl);
        return -1;
        //do a BFS to set-up the tip vertices layers
        auto visited_prop = request_vertex_property<bool>();
        std::vector<std::vector<VertexHandle>> tip_layers;
        std::vector<VertexHandle> to_visit = {witness_vertex_};
        while(!to_visit.empty()){
            SS_DEBUG(" -- gathering layer "<<(tip_layers.size()+1)<<" tips"<<std::endl);
            auto current_tip = to_visit.back();
            to_visit.pop_back();
            visited_prop[current_tip] = true;

            std::vector<VertexHandle> next_layer;
            for(auto out_he: outgoing_halfedges(current_tip)){
                auto neighbor = to_vertex_handle(out_he);
                if(tip_vertices_prop_[neighbor] &&
                        !visited_prop[neighbor]){
                    next_layer.push_back(neighbor);
                    visited_prop[neighbor] = true;
                }
            }

            SS_DEBUG(" -- layer "<<(tip_layers.size()+1)<<": "<<next_layer<<std::endl);
            if(!next_layer.empty()){
                tip_layers.push_back(next_layer);
            }
        }


        //then contract each layer one after the other
        //TODO

    }

    SS_DEBUG(" ...done with core contraction"<<std::endl);
    SS_DEBUG(" ---------------------------------"<<std::endl);

    return 0;
}


VertexHandle StarShapifyableExpansionCone::collapse_vertex(const VertexHandle& cone_v,
                                                           const VertexPosition& new_position){

    SS_DEBUG(" ------------ ");
    SS_DEBUG(" collapsing base vertex "<<cone_v<<" to position "<<vec2vec(new_position)<<"...");

    auto cone_v_prime = split_spoke_and_move_vertex(cone_v, new_position, {});
    if(!cone_v_prime.is_valid()){
        SS_DEBUG(" error while splitting spoke"<<std::endl);
        return VertexHandle(-1);
    }

    SS_DEBUG(" - updating collapse stack (secondary splits)");

    //then iterate through the collapse stack to do any necessary secondary split
    for(int i(collapse_list_.size()-1); i>= 0; i--){
        auto& collapse = collapse_list_[i];
        SS_DEBUG(" ------- handling collapse "<<i<<": "<<collapse.first<<" -> "<<collapse.second);

        bool needs_to_be_updated(false);
        for(auto& target_v: collapse.second){
            if(secondary_vertex_prop_[target_v].idx() != -1){
                needs_to_be_updated = true;
                target_v = secondary_vertex_prop_[target_v];
            }
        }

        if(needs_to_be_updated){
            SS_DEBUG(" ---> one of the target vertices has a secondary vertex, updating to "<<collapse.second);
            auto secondary_v_prime = split_spoke_and_move_vertex(collapse.first, get_position_from_vertices(collapse.second), {});
            if(!secondary_v_prime.is_valid()){
                SS_DEBUG(" error while splitting spoke"<<std::endl);
                return VertexHandle(-1);
            }
            collapse.first = secondary_v_prime;
            SS_DEBUG(" ---> updated collapse to "<<collapse.first<<" -> "<<collapse.second);
        }else{
            //SS_DEBUG(" ---> nothing was updated, skipping"<<std::endl);
        }

        //SS_DEBUG(" -------"<<std::endl);
    }


    SS_DEBUG(" ...done with stack for vertex "<<cone_v);
    SS_DEBUG(" ------------");


    return cone_v_prime;

}





VertexHandle StarShapifyableExpansionCone::get_last_secondary_vertex(const VertexHandle& vh) const{

    VertexHandle secondary_v(-1);
    VertexHandle next_secondary_v(secondary_vertex_prop_[vh]);

    //SS_DEBUG(" - "<<vh;
    while(next_secondary_v.idx() != -1){
        secondary_v = next_secondary_v;
        next_secondary_v = secondary_vertex_prop_[secondary_v];
        //SS_DEBUG(" -> "<<secondary_v;
    }
    //SS_DEBUG("---->"<<secondary_v<<std::endl);

    return secondary_v;
}


void StarShapifyableExpansionCone::update_vertices_with_secondaries(std::vector<VertexHandle>& vertices) const{
    for(auto& v: vertices){
        v = secondary_vertex_prop_[v];
    }
}

FaceHandle StarShapifyableExpansionCone::opposite_face(const VertexHandle& vh) const{

    auto ch = *vc_iter(vh);
    for(auto cf_it = cf_iter(ch); cf_it.valid(); cf_it++){
        auto f_vertices = get_halfface_vertices(halfface_handle(*cf_it, 0));
        if(f_vertices[0] != vh &&
                f_vertices[1] != vh &&
                f_vertices[2] != vh){
            return *cf_it;
        }
    }
    return FaceHandle(-1);
}


VertexHandle StarShapifyableExpansionCone::opposite_trim_base_vertex(const EdgeHandle& eh) const{

    for(auto ehf_it = trimmed_copy_.ehf_iter(eh); ehf_it.valid(); ehf_it++){
        if(!trimmed_copy_.is_boundary(*ehf_it)){
            auto op_v = trimmed_copy_.halfface_opposite_vertex(*ehf_it);
            if(!tip_vertices_prop_[op_v]){
                return op_v;
            }
        }
    }

    return VertexHandle(-1);

}

EdgeHandle StarShapifyableExpansionCone::find_target_edge(const VertexHandle& vh) const{
    //basically find the only face that is interior
    FaceHandle interior_face(-1);
    //SS_DEBUG(" - looking for target edge to collapse vertex "<<vh<<std::endl);
    for(auto vf_it = trimmed_copy_.vf_iter(vh); vf_it.valid(); vf_it++){
        if(!trimmed_copy_.is_boundary(*vf_it)){
            //SS_DEBUG(" -- checking boundary face "<<get_halfface_vertices(halfface_handle(*vf_it, 0))<<std::endl);
            if(interior_face.idx() != -1){
                SS_DEBUG(" ERROR - found two interior faces..."<<std::endl);
                return EdgeHandle(-1);
            }

            interior_face = *vf_it;
        }
    }
    if(interior_face.idx() == -1){
        SS_DEBUG(" ERROR - couldn't find an interior face"<<std::endl);
        return EdgeHandle(-1);
    }

    //then take the opposite edge on this face
    return opposite_edge_on_face(vh, interior_face);
}


VertexHandle StarShapifyableExpansionCone::find_RS_target_vertex(const VertexHandle& vh) const{
    VertexHandle target(-1);
    for(auto out_he: trimmed_copy_.outgoing_halfedges(vh)){
        /*SS_DEBUG(" - trying edge "<<halfedge(out_he)<<
                   ", base boundary: "<<trimmed_copy_.is_base_boundary_edge(edge_handle(out_he))<<
                   ", tip vertex: "<<tip_vertices_prop_[to_vertex_handle(out_he)]<<std::endl);*/
        if(!trimmed_copy_.is_base_boundary_edge(edge_handle(out_he)) &&
                !tip_vertices_prop_[to_vertex_handle(out_he)]){
            return to_vertex_handle(out_he);
        }
    }

    return VertexHandle(-1);
}




VertexPosition StarShapifyableExpansionCone::get_position_from_face(const FaceHandle& fh) const{
    auto f_vertices = get_halfface_vertices(halfface_handle(fh, 0));
    return get_position_from_vertices(f_vertices);
}

VertexPosition StarShapifyableExpansionCone::get_position_from_edge(const EdgeHandle& eh) const{
    return (vertex(edge(eh).from_vertex()) + vertex(edge(eh).to_vertex()))/2;
}

VertexPosition StarShapifyableExpansionCone::get_position_from_vertices(const std::vector<VertexHandle>& vertices) const{

    switch(vertices.size()){
    case 2: return (vertex(vertices[0]) + vertex(vertices[1]))/2;

    case 3: return ((vertex(vertices[0]) +
                    vertex(vertices[1]))/2 +
                vertex(vertices[2]))/2;

    case 4: return (vertex(vertices[0]) +
                vertex(vertices[1]) +
                vertex(vertices[2]) +
                vertex(vertices[3]))/4;
    default: {
        SS_DEBUG(" ERROR - unhandled case"<<std::endl);
        return {0,0,0};
    }
    }

}



VertexPosition StarShapifyableExpansionCone::get_contraction_position_from_vertices(const VertexHandle& vertex_to_contract,
                                                                                    const std::vector<VertexHandle>& target_vertices,
                                                                                    const VertexPropertyT<std::vector<VertexHandle>>& collapsed_to_prop){


    const bool use_new_version(false);

    if(use_new_version && !is_boundary(vertex_to_contract)){

        SS_DEBUG(" --> using new version for vertex "<<vertex_to_contract<<std::endl);

        ExpansionCone one_ring;
        auto tip = *(cone_tip_vertices().begin());

        auto spoke_edge = find_halfedge(vertex_to_contract, tip);
        if(!spoke_edge.is_valid()){
            SS_DEBUG(" ERROR - couldn't find spoke edge for vertex "<<vertex_to_contract<<std::endl);
            return {0,0,0};
        }

        auto added_vertex_prop = request_vertex_property<bool>();

        one_ring.add_vertex(vertex_to_contract, vertex(vertex_to_contract));
        one_ring.add_vertex(tip, vertex(tip));
        added_vertex_prop[vertex_to_contract] = true;
        added_vertex_prop[tip] = true;

        for(auto hec_it = hec_iter(spoke_edge); hec_it.valid(); hec_it++){

            //SS_DEBUG("   --- adding its vertices"<<std::endl);
            auto c_vertices = get_cell_vertices(*hec_it);
            for(auto cv: c_vertices){
                //SS_DEBUG("      ---- adding vertex "<<cv<<std::endl);
                if(!added_vertex_prop[cv]){
                    one_ring.add_vertex(cv, vertex_position_prop_[cv]);
                    added_vertex_prop[cv] = true;
                    //SS_DEBUG("        --> added vertex "<<cv<<std::endl);
                }else{
                    //SS_DEBUG("        --> already added "<<std::endl);
                }
            }

            auto cone_c = one_ring.add_cell(*hec_it, c_vertices);
            if(cone_c.idx() == -1){
                SS_DEBUG(" error while adding cell "<<c_vertices<<" incident to edge "<<halfedge(spoke_edge)<<std::endl);
                return {0,0,0};
            }
        }


        VertexPosition new_pos;
        auto exp = one_ring.is_geo_expandable(new_pos);
        if(exp){
            SS_DEBUG(" ERROR - contraction cone is not expandable"<<std::endl);
            one_ring.print_details();
            one_ring.is_geo_expandable(new_pos, true);

            return {0,0,0};
        }


        SS_DEBUG(" ----------------------------------------------"<<std::endl);
        return new_pos;

    }else{

        SS_DEBUG(" --> using old version for vertex "<<vertex_to_contract<<std::endl);

        if(target_vertices.size() == 2){

            auto target_face_vertices = target_vertices;
            target_face_vertices.push_back(vertex_to_contract);

            //SS_DEBUG(" - contracted to edge "<<target_face_vertices<<std::endl);

            return get_position_from_vertices(target_face_vertices);


            //SS_DEBUG(" WARNING - using new version for contraction with two target vertices"<<std::endl);

            auto tip_v(from_vertex_handle(base_vertex_to_spoke_edge_prop_[witness_vertex_]));
            std::vector<std::vector<VertexHandle>> halffaces;


            auto target_vertex = target_vertices[0] == tip_v ? target_vertices[1] : target_vertices[0];
            SS_DEBUG(" - contracting to edge "<<target_vertices<<std::endl);
            SS_DEBUG(" - target vertex is "<<target_vertex<<std::endl);
            /*for(auto vc_it = vc_iter(vertex_to_contract); vc_it.valid(); vc_it++){
                if(TopoHelper::cell_contains_vertex(*this, target_vertex, *vc_it)){
                    SS_DEBUG(" -- found target cell "<<
                }
            }*/

            auto mid_hf1 = find_halfface({target_vertex, vertex_to_contract, tip_v});
            auto op_vertex1 = halfface_opposite_vertex(mid_hf1);
            halffaces.push_back({target_vertex, op_vertex1, tip_v});
            SS_DEBUG(" - added side face "<<halffaces.back()<<" to constraints"<<std::endl);
            halffaces.push_back({op_vertex1, vertex_to_contract, tip_v});
            SS_DEBUG(" - added side face "<<halffaces.back()<<" to constraints"<<std::endl);
            halffaces.push_back({op_vertex1, target_vertex, vertex_to_contract});
            SS_DEBUG(" - added top face "<<halffaces.back()<<" to constraints"<<std::endl);


            auto mid_hf2 = opposite_halfface_handle(mid_hf1);
            auto op_vertex2 = halfface_opposite_vertex(mid_hf2);
            halffaces.push_back({op_vertex2, target_vertex, tip_v});
            SS_DEBUG(" - added side face "<<halffaces.back()<<" to constraints"<<std::endl);
            halffaces.push_back({vertex_to_contract, op_vertex2, tip_v});
            SS_DEBUG(" - added side face "<<halffaces.back()<<" to constraints"<<std::endl);
            halffaces.push_back({op_vertex2, vertex_to_contract, target_vertex});
            SS_DEBUG(" - added top face "<<halffaces.back()<<" to constraints"<<std::endl);

            auto v = vertex(vertex_to_contract);

            //then add the constraints from the vertices that were collapsed to this one
            for(auto out_he: outgoing_halfedges(vertex_to_contract)){
                auto neighbor = to_vertex_handle(out_he);
                if(is_cone_tip(neighbor)/* || neighbor.idx() > vertex_to_contract.idx()*/){
                    continue;
                }

                if(!secondary_vertex_prop_[neighbor].is_valid()){
                    //SS_DEBUG(" -- neighbor vertex "<<neighbor<<" has a secondary, skipping"<<std::endl);
                    //continue;
                }

                VertexPosition normal, centroid;
                triangle_normal_and_centroid({neighbor, witness_vertex_, tip_v},
                                             normal,
                                             centroid);

                auto pv = v - centroid;

                auto side = CGAL::orientation(OVMvec3ToCGALPoint3(vertex(neighbor)),
                                              OVMvec3ToCGALPoint3(vertex(witness_vertex_)),
                                              OVMvec3ToCGALPoint3(vertex(tip_v)),
                                              OVMvec3ToCGALPoint3(v));

                //if(pv.dot(normal) < 0){
                if(side == CGAL::NEGATIVE){
                    halffaces.push_back({neighbor, witness_vertex_, tip_v});
                    SS_DEBUG(" -- added face "<<halffaces.back()<<std::endl);
                    //}else if(pv.dot(normal) > 0){
                }else if(side == CGAL::POSITIVE){
                    halffaces.push_back({neighbor, tip_v, witness_vertex_});
                    SS_DEBUG(" -- added face "<<halffaces.back()<<std::endl);
                }else{
                    SS_DEBUG(" --> vertex to contract "<<vertex_to_contract<<" is on the plane spanned by vertex "<<neighbor<<", skipping"<<std::endl);
                }
            }



            VertexPosition new_pos;
            int cheb_result = find_chebyshev_center_from_halffaces(halffaces,
                                                                   new_pos);

            if(cheb_result){
                SS_DEBUG(" ERROR - couldn't find contraction Chebyshev center, details:"<<std::endl);
                find_chebyshev_center_from_halffaces(halffaces,
                                                     new_pos,
                                                     true);

                return {0,0,0};
            }

            SS_DEBUG(" ----------------------------------------------"<<std::endl);
            return new_pos;

        }else if(target_vertices.size() == 3){


            auto tip_v(from_vertex_handle(base_vertex_to_spoke_edge_prop_[witness_vertex_]));

            CGAL::Triangle_3<CGAL_ExactKernel> wtv(OVMvec3ToCGALPoint3(vertex(witness_vertex_)),
                                                   OVMvec3ToCGALPoint3(vertex(tip_v)),
                                                   OVMvec3ToCGALPoint3(vertex(vertex_to_contract)));

            CGAL::Triangle_3<CGAL_ExactKernel> target_tri(OVMvec3ToCGALPoint3(vertex(target_vertices[0])),
                    OVMvec3ToCGALPoint3(vertex(target_vertices[1])),
                    OVMvec3ToCGALPoint3(vertex(target_vertices[2])));

            //SS_DEBUG(" wtv = "<<witness_vertex_<<", "<<tip_v<<", "<<vertex_to_contract<<std::endl);
            //SS_DEBUG(" abc = "<<target_vertices<<std::endl);

            auto inter = CGAL::intersection(wtv, target_tri);

            if(!inter){
                SS_DEBUG(" ERROR - no intersection between w-t-v triangle and target face"<<std::endl);
                SS_DEBUG(" - w-t-v: "<<wtv<<std::endl);
                SS_DEBUG(" --- = "<<vec2vec(vertex(witness_vertex_))<<", "<<vec2vec(vertex(tip_v))<<", "<<vec2vec(vertex(vertex_to_contract))<<std::endl);
                SS_DEBUG(" - target face: "<<target_vertices<<std::endl);
                SS_DEBUG(" --- = "<<vec2vec(vertex(target_vertices[0]))
                        <<", "<<vec2vec(vertex(target_vertices[1]))
                        <<", "<<vec2vec(vertex(target_vertices[2]))<<std::endl);
                return {0,0,0};
            }

            using CGAL_Segment_3 = CGAL::Segment_3<CGAL_ExactKernel>;
            CGAL_Segment_3* segment_result = boost::get<CGAL_Segment_3>(&*inter);

            if(!segment_result){
                SS_DEBUG(" ERROR - intersection between w-t-v triangle and target face "<<target_vertices<<" is not a segment "<<std::endl);
                return {0,0,0};
            }

            return (vertex(vertex_to_contract) +
                    (OVMvec3ToCGALPoint3(segment_result->point(0)) +
                    OVMvec3ToCGALPoint3(segment_result->point(1)))/2)/2;

        }else{
            SS_DEBUG(" ERROR - "<<target_vertices.size()<<" target vertices, so neither a face nor an edge"<<std::endl);
            return {0,0,0};
        }
    }

}




void StarShapifyableExpansionCone::translate_to_have_tips_at_origin(){
    auto tip_pos = vertex(*cone_tip_vertices_.begin());

    for(auto v: vertices()){
        set_vertex(v, vertex(v) - tip_pos);
    }
}


void StarShapifyableExpansionCone::find_unvisible_vertices_from(const VertexHandle& witness_vertex,
                                                                std::vector<VertexHandle>& unvisible_vertices) const{
    SS_DEBUG(" - not visible: ");
    for(auto to_vertex: vertices()){
        if(to_vertex != witness_vertex){
            if(!is_visible_from(witness_vertex, to_vertex)){
                SS_DEBUG(" "<<to_vertex);
                unvisible_vertices.push_back(to_vertex);
            }
        }
    }
    SS_DEBUG(std::endl);
}


void StarShapifyableExpansionCone::gather_vertices_outside_of_1_ring_neighborhood_of_witness_vertex(std::vector<VertexHandle>& non_neighbor_vertices){

    non_neighbor_vertices.clear();
    auto neighbor_prop = request_vertex_property<bool>();
    for(auto out_he: outgoing_halfedges(witness_vertex_)){
        neighbor_prop[to_vertex_handle(out_he)] = true;
    }

    for(auto v: vertices()){
        if(v != witness_vertex_ &&
                !tip_vertices_prop_[v] &&
                !neighbor_prop[v]){
            non_neighbor_vertices.push_back(v);
        }
    }
}

bool StarShapifyableExpansionCone::is_integrity_maintained(const VertexHandle& new_mid_vertex,
                                                           const VertexPropertyT<CGAL::Orientation>& neighbors_side_prop,
                                                           const std::vector<VertexHandle>& target_vertices,
                                                           const CellPropertyT<bool>& degenerate_tet_prop,
                                                           const std::vector<std::pair<std::vector<CGAL_ExactPoint3>, CGAL::Sign>>& additional_constraints,
                                                           bool print_debug){

    if(print_debug){
        int max, total;
        float avg;
        one_ring_neighborhood_byte_size(new_mid_vertex, max, avg, total);
        SS_DEBUG(" ......................................................................"<<std::endl);
        SS_DEBUG(" --- checking integrity for vertex "<<new_mid_vertex<<std::endl);
        SS_DEBUG(" --- 1-ring max byte size = "<<max<<std::endl);
        SS_DEBUG(" --- 1-ring avg byte size = "<<avg<<std::endl);
        SS_DEBUG(" --- 1-ring total byte size = "<<total<<std::endl);
        SS_DEBUG(" --- #target vertices = "<<target_vertices.size()<<std::endl);
        SS_DEBUG(" --- #vertices to contract = "<<vertices_to_contract_.size()<<std::endl);
        SS_DEBUG(" --- vertex valence = "<<valence(new_mid_vertex)<<std::endl);
    }



    auto flipped_check_start = std::chrono::high_resolution_clock::now();
    if(ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(*this, vertex_position_prop_, new_mid_vertex)){
    //if(ExactBadTetFinder::meshContainsFlippedTets(*this, vertex_position_prop_)){

        if(print_debug){
            SS_DEBUG(" -> mesh contains flipped tets"<<std::endl);
            for(auto c: ExactBadTetFinder::findBadTets(*this, vertex_position_prop_).second){
                SS_DEBUG(" - "<<c<<": "<<get_cell_vertices(c)<<std::endl);
            }
            SS_DEBUG(" ......................................................................"<<std::endl);
        }
        return false;
    }
    if(print_debug){
        auto flipped_check_end = std::chrono::high_resolution_clock::now();
        float flipped_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(flipped_check_end - flipped_check_start).count() / 1000000;
        SS_DEBUG(" --- flipped check duration: "<<flipped_check_duration_s<<std::endl);
    }

    auto deg_check_start = std::chrono::high_resolution_clock::now();
    for(auto vc_it = vc_iter(new_mid_vertex); vc_it.valid(); vc_it++){
        auto is_deg = OVMtetToCGALtet(*this, vertex_position_prop_, *vc_it).is_degenerate();
        if(degenerate_tet_prop[*vc_it] != is_deg){
            if(print_debug){
                SS_DEBUG(" -> tet "<<get_cell_vertices(*vc_it)<<" initial deg: "<<degenerate_tet_prop[*vc_it]<<" and now: "<<is_deg<<std::endl);
                SS_DEBUG(" ......................................................................"<<std::endl);
            }
            return false;
        }
    }
    if(print_debug){
        auto deg_check_end = std::chrono::high_resolution_clock::now();
        float deg_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(deg_check_end - deg_check_start).count() / 1000000;
        SS_DEBUG(" --- deg check duration: "<<deg_check_duration_s<<std::endl);
    }

    auto tip_v = *cone_tip_vertices().begin();//from_vertex_handle(base_vertex_to_spoke_edge_prop_[new_mid_vertex]);
    if(print_debug){
        //SS_DEBUG(" tip vertex = "<<tip_v<<std::endl);
    }


    /*
    auto star_shape_check_start = std::chrono::high_resolution_clock::now();
    ExpansionCone updated_cone;
    set_up_1_ring_neighborhood_as_expansion_cone(*this,
                                                 vertex_position_prop_,
                                                 tip_v,
                                                 updated_cone);
    VertexPosition pos;
    auto geo_exp_result = updated_cone.is_geo_expandable(pos);

    if(geo_exp_result){

        if(print_debug){
            SS_DEBUG(" -> EC is no longer star-shaped"<<std::endl);
            SS_DEBUG(" ......................................................................"<<std::endl);
        }
        //updated_cone.print_details();
        //SS_DEBUG(" -> no longer star-shaped"<<std::endl);
        return false;
    }
    if(print_debug){
        auto star_shape_check_end = std::chrono::high_resolution_clock::now();
        float star_shape_check_end_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(star_shape_check_end - star_shape_check_start).count() / 1000000;
        SS_DEBUG(" --- star-shape check duration: "<<star_shape_check_end_duration_s<<std::endl);
    }
    */

    auto v = vertex(new_mid_vertex);

    //SS_DEBUG(" - v = "<<v<<std::endl);
    //then add the constraints from the vertices that were collapsed to this one
    //for(auto out_he: outgoing_halfedges(new_mid_vertex)){
    // auto vh = to_vertex_handle(out_he);
    auto sides_check_start = std::chrono::high_resolution_clock::now();
    for(auto vh: vertices_to_contract_){
        if(is_cone_tip(vh) || vh == witness_vertex_/* || neighbor.idx() > vertex_to_contract.idx()*/){
            continue;
        }

        if(secondary_vertex_prop_[vh].is_valid()){
            //SS_DEBUG(" -- vertex "<<vh<<" has a secondary, skipping"<<std::endl);
            continue;
        }

        //SS_DEBUG(" -- vertex "<<vh<<" initial side: "<<neighbors_side_prop[vh]<<std::endl);
        if(!neighbors_side_prop[vh]){
            //SS_DEBUG(" -- neighbor "<<vh<<" was aligned, skipping"<<std::endl);
            continue;
        }

        /*VertexPosition normal, centroid;
        triangle_normal_and_centroid({vh, witness_vertex_, tip_v},
                                     normal,
                                     centroid);

        const auto& pv = v - centroid;
        const auto& dot = pv.dot(normal);
        int side = dot < 0 ? -1 : (dot > 0 ? 1 : 0);*/

        auto side = CGAL::orientation(OVMvec3ToCGALPoint3(vertex(vh)),
                                      OVMvec3ToCGALPoint3(vertex(witness_vertex_)),
                                      OVMvec3ToCGALPoint3(vertex(tip_v)),
                                      OVMvec3ToCGALPoint3(v));

        //SS_DEBUG(" -- side for base vertex "<<vh<<": "<<side<<std::endl);

        if(print_debug){
            //SS_DEBUG(" -- side of vertex "<<vh<<": "<<side<<std::endl);
        }
        //SS_DEBUG("  -- "<<vh<<" normal = "<<normal<<", centroid = "<<centroid<<", dot = "<<dot<<", side: "<<side<<std::endl);
        if(side != neighbors_side_prop[vh]){

            if(print_debug){
                SS_DEBUG(" --> vertex "<<new_mid_vertex<<" is not on the same side of vertex "<<vh<<" anymore"<<std::endl);
                SS_DEBUG(" ......................................................................"<<std::endl);
            }
            return false;
        }
    }
    if(print_debug){
        auto sides_check_end = std::chrono::high_resolution_clock::now();
        float sides_check_end_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(sides_check_end - sides_check_start).count() / 1000000;
        SS_DEBUG(" --- sides check duration: "<<sides_check_end_duration_s<<std::endl);
    }

    auto target_sides_check_start = std::chrono::high_resolution_clock::now();
    for(auto vh: target_vertices){
        if(is_cone_tip(vh)){
            continue;
        }
        /*VertexPosition normal, centroid;
        triangle_normal_and_centroid({vh, witness_vertex_, tip_v},
                                     normal,
                                     centroid);

        const auto& pv = v - centroid;
        const auto& dot = pv.dot(normal);
        int side = dot < 0 ? -1 : (dot > 0 ? 1 : 0);*/
        auto side = CGAL::orientation(OVMvec3ToCGALPoint3(vertex(vh)),
                                      OVMvec3ToCGALPoint3(vertex(witness_vertex_)),
                                      OVMvec3ToCGALPoint3(vertex(tip_v)),
                                      OVMvec3ToCGALPoint3(v));

        //SS_DEBUG(" -- side for base target vertex "<<vh<<": "<<side<<std::endl);
        if(side != neighbors_side_prop[vh]){
            if(print_debug){
                SS_DEBUG(" --> target vertex "<<new_mid_vertex<<" is not on the same side of vertex "<<vh<<" anymore"<<std::endl);
                SS_DEBUG(" ......................................................................"<<std::endl);
            }
            return false;
        }
    }
    if(print_debug){
        auto target_sides_check_end = std::chrono::high_resolution_clock::now();
        float target_sides_check_end_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(target_sides_check_end - target_sides_check_start).count() / 1000000;
        SS_DEBUG(" --- target sides check duration: "<<target_sides_check_end_duration_s<<std::endl);
    }


    int i(0);
    auto additional_sides_check_start = std::chrono::high_resolution_clock::now();
    for(const auto& constr: additional_constraints){

        auto side = CGAL::orientation(constr.first[0], constr.first[1], constr.first[2], OVMvec3ToCGALPoint3(v));

        if(side != constr.second){

            if(print_debug){
                SS_DEBUG(" --> switched side from "<<constr.second<<" to "<<side<<" for additional constraint n°"<<i<<std::endl);
                SS_DEBUG(" ......................................................................"<<std::endl);
            }
            /*SS_DEBUG("     points used: "<<std::endl);
            for(const auto& p: constr.first){
                SS_DEBUG("     - "<<p<<std::endl);
            }
            std::string ans;
            std::cin>>ans;*/
            return false;
        }
        i++;
    }

    if(print_debug){
        auto additional_sides_check_end = std::chrono::high_resolution_clock::now();
        float additional_sides_check_end_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(additional_sides_check_end - additional_sides_check_start).count() / 1000000;
        SS_DEBUG(" --- additional sides check duration: "<<additional_sides_check_end_duration_s<<std::endl);
    }

    if(print_debug){
        SS_DEBUG(" --> integrity maintained"<<std::endl);
        SS_DEBUG(" ......................................................................"<<std::endl);
    }
    return true;
}

void StarShapifyableExpansionCone::one_ring_neighborhood_byte_size(const VertexHandle& vertex_to_contract,
                                                                   int& max_byte_size,
                                                                   float& average_byte_size,
                                                                   int& total_byte_size) const{
    max_byte_size = 0;
    average_byte_size = 0;
    total_byte_size = 0;

    float n(0);
    for(auto vv_it = vv_iter(vertex_to_contract); vv_it.valid(); vv_it++){
        auto size = byte_size(vertex(*vv_it));
        max_byte_size = std::max(max_byte_size, size);
        total_byte_size += size;
        n++;
    }

    average_byte_size = (float)total_byte_size / n;
}

VertexPosition StarShapifyableExpansionCone::binary_search_min_precision_position(const VertexHandle& new_mid_vertex,
                                                                                  const std::vector<VertexHandle>& target_vertices){

    //return vertex(new_mid_vertex);


    SS_DEBUG(" =======================================================");
    SS_DEBUG("  looking for minimum precision position for new mid-vertex "<<new_mid_vertex<<" initially at "<<std::setprecision(20)<<vec2vec(vertex(new_mid_vertex))<<" with binary search...");

#if ENABLE_TIMINGS
    SS_DEBUG(" ======================================================="<<std::endl);
    SS_DEBUG("  looking for minimum precision position for new mid-vertex "<<new_mid_vertex<<" initially at "<<vec2vec(vertex(new_mid_vertex))<<" with binary search..."<<std::endl);
#endif

    const int byte_size_lower_bound(400);
    const int min_shift_factor_diff(16);
    const int min_byte_size_diff(byte_size_lower_bound/2);
    const int max_iterations(10);

    auto new_vertex_initial_position = vertex(new_mid_vertex);

    auto initial_byte_size = byte_size(new_vertex_initial_position);
    if(initial_byte_size < byte_size_lower_bound){

#if ENABLE_TIMINGS
        SS_DEBUG(" --> initial precision lower than "<<byte_size_lower_bound<<", using that"<<std::endl);
        SS_DEBUG(" ======================================================="<<std::endl);
#endif

        SS_DEBUG(" --> initial precision enough lower than "<<byte_size_lower_bound<<", using that");
        SS_DEBUG(" =======================================================");
        return new_vertex_initial_position;
    }

    auto setup_start = std::chrono::high_resolution_clock::now();

    //SS_DEBUG(" - initial mid-vertex position = "<<new_vertex_initial_position<<std::endl);
    //SS_DEBUG(" - initial sides: "<<std::endl);
    //auto neighbors_side_prop = request_vertex_property<int>();
    auto exact_neighbors_side_prop = request_vertex_property<CGAL::Orientation>();

    auto tip_v = *cone_tip_vertices().begin();
    auto v = OVMvec3ToCGALPoint3(new_vertex_initial_position);
    auto w = OVMvec3ToCGALPoint3(vertex(witness_vertex_));
    auto t = OVMvec3ToCGALPoint3(vertex(tip_v));

    //then add the constraints from the vertices that were collapsed to this one
    for(auto vh: vertices_to_contract_){
        if(is_cone_tip(vh) || vh == witness_vertex_/* || neighbor.idx() > vertex_to_contract.idx()*/){
            continue;
        }

        if(secondary_vertex_prop_[vh].is_valid()){
            //SS_DEBUG(" -- vertex "<<vh<<" has a secondary, skipping"<<std::endl);
            continue;
        }

        exact_neighbors_side_prop[vh] = CGAL::orientation(OVMvec3ToCGALPoint3(vertex(vh)), w, t, v);

        //SS_DEBUG(" -- side for base vertex "<<vh<<": "<<exact_neighbors_side_prop[vh]<<std::endl);
    }

    //SS_DEBUG(" - target vertices: "<<target_vertices<<std::endl);
    int target_side_sum(0);
    for(auto vh: target_vertices){
        if(!is_cone_tip(vh)){

            exact_neighbors_side_prop[vh] = CGAL::orientation(OVMvec3ToCGALPoint3(vertex(vh)), w, t, v);
            //SS_DEBUG(" -- side for base target vertex "<<vh<<": "<<exact_neighbors_side_prop[vh]<<std::endl);

        }else{
            //by convention, for next check
            //neighbors_side_prop[vh] = 0;
            exact_neighbors_side_prop[vh] = CGAL::Sign(0);
        }
        //target_side_sum += neighbors_side_prop[vh];
        target_side_sum += exact_neighbors_side_prop[vh];
    }

    if(target_vertices.size() == 3 && target_side_sum){
        SS_DEBUG(" ERROR - target sides sum is not zero. Sides:"<<std::endl);
        for(auto target:target_vertices){
            //SS_DEBUG(" - "<<target<<" :" <<neighbors_side_prop[target]<<std::endl);
            SS_DEBUG(" - "<<target<<" :" <<exact_neighbors_side_prop[target]<<std::endl);

        }
        return {0,0,0};
    }

    auto degenerate_tet_prop = request_cell_property<bool>();
    for(auto vc_it = vc_iter(new_mid_vertex); vc_it.valid(); vc_it++){
        degenerate_tet_prop[*vc_it] = OVMtetToCGALtet(*this, vertex_position_prop_, *vc_it).is_degenerate();
    }


    std::vector<std::pair<std::vector<CGAL_ExactPoint3>, CGAL::Sign>> additional_constraints;
    int i = 0;
    for(auto vv_it = vv_iter(new_mid_vertex); vv_it.valid(); vv_it++){
        if(is_cone_tip(*vv_it)){
            continue;
        }
        for(auto target_v: target_vertices){

            if(target_v == *vv_it || is_cone_tip(target_v)){
                continue;
            }

            auto target = OVMvec3ToCGALPoint3(vertex(target_v));
            auto s = OVMvec3ToCGALPoint3(vertex(*vv_it));
            auto side = CGAL::orientation(s, target, t, v);

            if(!side){
                //SS_DEBUG(" --- skipping additional constraint "<<(*vv_it)<<" - "<<target_v<<" with side = "<<side<<std::endl);
                continue;
            }
            additional_constraints.push_back({{s, target, t}, side});
            //SS_DEBUG(" --- added additional constraint n°"<<i<<": "<<(*vv_it)<<" - "<<target_v<<" with side = "<<side<<std::endl);
            /*SS_DEBUG("     points used: "<<std::endl);
            for(const auto& p: additional_constraints.back().first){
                SS_DEBUG("     - "<<p<<std::endl);
            }*/
            i++;
        }

    }


#if ENABLE_TIMINGS
    auto setup_end = std::chrono::high_resolution_clock::now();
    float setup_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start).count() / 1000000;
    SS_DEBUG(" -- set-up duration: "<<setup_duration_s<<std::endl);
    auto first_integrity_check_start = std::chrono::high_resolution_clock::now();
#endif


    if(!is_integrity_maintained(new_mid_vertex, exact_neighbors_side_prop, target_vertices, degenerate_tet_prop, additional_constraints, ENABLE_TIMINGS)){
        SS_DEBUG(" ERROR - integrity already not maintained with initial position "<<vec2vec(new_vertex_initial_position)<<" for vertex "<<new_mid_vertex<<std::endl);
        SS_DEBUG(" - list of remaining vertices to contract: "<<vertices_to_contract_<<std::endl);
        is_integrity_maintained(new_mid_vertex, exact_neighbors_side_prop, target_vertices, degenerate_tet_prop, additional_constraints, true);

        //found_deg_tets_ = true;
        //return new_vertex_initial_position;
        return {0,0,0};
    }

#if ENABLE_TIMINGS
    auto first_integrity_check_end = std::chrono::high_resolution_clock::now();
    float first_integrity_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(first_integrity_check_end - first_integrity_check_start).count() / 1000000;
    SS_DEBUG(" -- initial position integrity check duration: "<<first_integrity_check_duration_s<<std::endl);
#endif

    //SS_DEBUG(" WARNING - precision reduction diabled"<<std::endl);
    //return new_vertex_initial_position;


    //first, try by simply using double precision
    const auto initial_pos = this->vertex(new_mid_vertex);
    auto current_pos = vec2vec(vec2vec(new_vertex_initial_position));
    this->set_vertex(new_mid_vertex, current_pos, false);

    auto double_integrity_check_start = std::chrono::high_resolution_clock::now();
    if(is_integrity_maintained(new_mid_vertex, exact_neighbors_side_prop, target_vertices, degenerate_tet_prop, additional_constraints, ENABLE_TIMINGS)){

#if ENABLE_TIMINGS
        SS_DEBUG(" --> double precision enough, using that"<<std::endl);
        SS_DEBUG(" ======================================================="<<std::endl);
#endif
        SS_DEBUG(" --> double precision enough, using that");
        SS_DEBUG(" =======================================================");
        return current_pos;
    }

#if ENABLE_TIMINGS
    auto double_integrity_check_end = std::chrono::high_resolution_clock::now();
    float double_integrity_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(double_integrity_check_end - double_integrity_check_start).count() / 1000000;
    SS_DEBUG(" -- double integrity check duration: "<<double_integrity_check_duration_s<<std::endl);
#endif

    current_pos = new_vertex_initial_position;

    //SS_DEBUG(" - initial byte size = "<<initial_byte_size<<std::endl);


    int max_shift_factor = PEHelpers::find_maximum_shift_factor(new_vertex_initial_position);
    int min_shift_factor = 0;

    int shift_factor = max_shift_factor;

    bool found_best_factor(false);
    int last_byte_size(initial_byte_size);
    int last_shift_factor(shift_factor);

    i = 0;
    while(!found_best_factor && i < max_iterations){

        current_pos = new_vertex_initial_position;
        PEHelpers::lower_precision(shift_factor, current_pos);
        int current_byte_size = byte_size(current_pos);

#if ENABLE_TIMINGS
        SS_DEBUG(" ------ "<<std::endl);
        SS_DEBUG(" ---- byte size at iteration "<<i<<"/"<<max_iterations<<" : "<<byte_size(current_pos)<<std::endl);
        SS_DEBUG("          min shift factor = "<<min_shift_factor<<std::endl);
        SS_DEBUG("          max shift factor = "<<max_shift_factor<<std::endl);
        SS_DEBUG("              shift factor = "<<shift_factor<<std::endl);
#endif

        this->set_vertex(new_mid_vertex, current_pos, false);

        auto integrity_check_start = std::chrono::high_resolution_clock::now();
        //if the precision is too low, we use it as lower bound
        if(!is_integrity_maintained(new_mid_vertex, exact_neighbors_side_prop, target_vertices, degenerate_tet_prop, additional_constraints, ENABLE_TIMINGS)){
            max_shift_factor = shift_factor;

            //if it's too high, we use it as UPPER bound
        }else{
            min_shift_factor = shift_factor;
        }
#if ENABLE_TIMINGS
        auto integrity_check_end = std::chrono::high_resolution_clock::now();
        float integrity_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(integrity_check_end - integrity_check_start).count() / 1000000;
        SS_DEBUG(" -- "<<(i+1)<<"-th integrity check duration: "<<integrity_check_duration_s<<std::endl);
#endif

        if((shift_factor <= min_shift_factor && current_byte_size < byte_size_lower_bound) ||
                std::abs(current_byte_size - last_byte_size) < min_byte_size_diff){
            //SS_DEBUG(" --> ok"<<std::endl);
            found_best_factor = true;
        }

        last_byte_size = current_byte_size;
        last_shift_factor = shift_factor;
        shift_factor = (min_shift_factor + max_shift_factor) / 2;
        if(std::abs(last_shift_factor - shift_factor) < min_shift_factor_diff){
            //SS_DEBUG(" --> shift factor change < "<<min_shift_factor_diff<<" -> ok"<<std::endl);
            found_best_factor = true;
        }
        //SS_DEBUG(" - updated shift factor = "<<shift_factor<<std::endl);

        i++;
    }

    //SS_DEBUG(" final min shift factor = "<<min_shift_factor<<std::endl);

    //re-compute best valid position
    current_pos = new_vertex_initial_position;
    PEHelpers::lower_precision(min_shift_factor, current_pos);

    //SS_DEBUG(" final byte size = "<<byte_size(current_pos)<<std::endl);

    this->set_vertex(new_mid_vertex, initial_pos);


    SS_DEBUG(" ...done, reduced precision from "<<initial_byte_size<<" bytes to "<<byte_size(current_pos));
    SS_DEBUG(" =======================================================");

    /*if(byte_size(current_pos) < initial_byte_size){
        SS_DEBUG(" --> actually reduced here"<<std::endl);
        std::string ans;
        std::cin>>ans;
    }*/

#if ENABLE_TIMINGS
    SS_DEBUG(" ...done, reduced precision from "<<initial_byte_size<<" bytes to "<<byte_size(current_pos)<<" after "<<i<<" iterations"<<std::endl);
    SS_DEBUG(" ======================================================="<<std::endl);
#endif
    //std::string ans;
    //std::cin>>ans;

    return current_pos;
}




VertexHandle StarShapifyableExpansionCone::split_spoke_and_move_vertex(const VertexHandle& cone_v,
                                                                       const VertexPosition& new_position,
                                                                       const std::vector<VertexHandle>& target_vertices,
                                                                       bool reduce_precision){


    //SS_DEBUG(" ------------ ");
    //SS_DEBUG(" splitting spoke of base vertex "<<cone_v<<" to position "<<vec2vec(new_position)<<"...");

    if(new_position == VertexPosition(0,0,0)){
        SS_DEBUG(" ERROR - trying to move spoke mid-vertex to origin. This should never happen"<<std::endl);
        return VertexHandle(-1);
    }



    const auto cone_spoke_heh = base_vertex_to_spoke_edge_prop_[cone_v];

    //then, split the spoke (both on cone and mesh)
    const auto cone_tip       = from_vertex_handle(cone_spoke_heh);
    const auto cone_v_prime   = split_edge(cone_spoke_heh);
    base_vertex_to_spoke_edge_prop_[cone_v_prime] = find_halfedge(cone_tip, cone_v_prime);
    collapsed_prop_[cone_v] = true;
    secondary_vertex_prop_[cone_v] = cone_v_prime;

    original_base_vertex_prop_[cone_v_prime] = original_base_vertex_prop_[to_vertex_handle(cone_spoke_heh)];
    //SS_DEBUG(" - set original base vertex of mid-vertex "<<cone_v_prime<<" as "<<original_base_vertex_prop_[cone_v_prime]<<std::endl);
    mid_vertices_positions_prop_[original_base_vertex_prop_[cone_v_prime]].push_back(new_position);
    mid_vertices_prop_[original_base_vertex_prop_[cone_v_prime]].push_back(cone_v_prime);
    //SS_DEBUG(" - and added mid-vertex position "<<vec2vec(min_prec_new_pos)<<" to list of mid-vertices positions for base vertex "<<original_base_vertex_prop_[cone_v_prime]<<std::endl);
    //SS_DEBUG(" - mid-vertices for base vertex "<<original_base_vertex_prop_[cone_v_prime]<<": "<<mid_vertices_prop_[original_base_vertex_prop_[cone_v_prime]]<<std::endl);

    //SS_DEBUG(" - split edge "<<halfedge(cone_spoke_heh)<<" -> "<<cone_v_prime<<std::endl);

    set_vertex(cone_v_prime, new_position);

    if(ENABLE_PRECISION_REDUCTION &&
            reduce_precision &&
            !is_boundary(cone_v_prime)){
        auto min_prec_new_pos = binary_search_min_precision_position(cone_v_prime, target_vertices);
        if(min_prec_new_pos == VertexPosition(0,0,0)){
            SS_DEBUG(" ERROR - min-precision position is the origin. This should never happen"<<std::endl);
            SS_DEBUG(" failed to split spoke for vertex "<<cone_v<<" with target vertices "<<target_vertices<<std::endl);
            return VertexHandle(-1);
        }
        set_vertex(cone_v_prime, min_prec_new_pos);

        saved_position_bytes_ += (byte_size(new_position) - byte_size(min_prec_new_pos));
    }else {
        SS_DEBUG(" - skipping precision reduction");
        SS_DEBUG("     --            macro: "<<ENABLE_PRECISION_REDUCTION);
        SS_DEBUG("     -- reduce precision: "<<reduce_precision);
        SS_DEBUG("     --         boundary: "<<is_boundary(cone_v_prime));
        //SS_DEBUG(" - skipping precision reduction (boundary: "<<is_boundary(cone_v_prime)<<")"<<std::endl);
    }



    //add this new split to the list
    split_list_.push_back({cone_tip,
                           cone_v,
                           cone_v_prime,
                           vertex(cone_v_prime)});

    //if(max_byte_size(new_position) > 100){
        //SS_DEBUG(" --> moved vertex "<<cone_v_prime<<" to position with byte size = "<<max_byte_size(new_position)<<std::endl);
    //}
    SS_DEBUG(" - "<<cone_v<<"' = "<<cone_v_prime<<" at "<<vec2vec(vertex(cone_v_prime)));
    //SS_DEBUG(" ------------");
    return cone_v_prime;
}


//NOTE: the details are done in the HalfEdge version because it can sometimes be necessary to
//know which of the from- and to-vertex is which
VertexHandle StarShapifyableExpansionCone::split_base_edge_and_register_split(const HalfEdgeHandle& to_split){

    auto cone_from_v(from_vertex_handle(to_split));
    auto cone_to_v(to_vertex_handle(to_split));
    auto mid_vertex = split_edge(to_split);
    auto mid_pos = (vertex(cone_from_v) + vertex(cone_to_v))/2;
    //to update the exact vertex position prop
    this->set_vertex(mid_vertex, mid_pos);


    split_list_.push_back({cone_from_v,
                          cone_to_v,
                          mid_vertex,
                          mid_pos});

    SS_DEBUG(" --> performed split "<<halfedge(to_split)<<" -> "<<mid_vertex<<" at "<<vec2vec(mid_pos)<<" and added it to list");

    original_base_vertex_prop_[mid_vertex] = mid_vertex;
    //SS_DEBUG(" --> set "<<mid_vertex<<" as original base vertex of itself");

    if(original_edge_prop_[edge_handle(to_split)]){
        auto first_half_e = edge_handle(find_halfedge(cone_from_v, mid_vertex));
        auto second_half_e = edge_handle(find_halfedge(cone_to_v,   mid_vertex));
        original_edge_prop_[first_half_e] = true;
        original_edge_prop_[second_half_e] = true;
        SS_DEBUG(" --> marked new edges "<<edge(first_half_e)<<" and "<<edge(second_half_e)<<" as pieces of and original edge");
    }else{
        SS_DEBUG(" --> not an original edge");
    }

    return mid_vertex;
}


bool StarShapifyableExpansionCone::is_not_part_of_witness_vertex_1_ring(const EdgeHandle& eh) const{
    return edge(eh).from_vertex() != witness_vertex_ &&
           edge(eh).to_vertex()   != witness_vertex_ &&
            (!witness_vertex_1_ring_base_vertex_prop_[edge(eh).from_vertex()]
            || !witness_vertex_1_ring_base_vertex_prop_[edge(eh).to_vertex()]);

}



void StarShapifyableExpansionCone::count_base_edges(int& original_edges_count,
                                                    int& new_edges_count) const{
    original_edges_count = 0;
    new_edges_count = 0;
    for(auto e: edges()){
        original_edges_count += original_edge_prop_[e] && is_base_edge(e);
        new_edges_count += !original_edge_prop_[e] && is_base_edge(e);
    }
}


VertexHandle StarShapifyableExpansionCone::split_base_edge_and_register_split(const EdgeHandle& edge){
    return split_base_edge_and_register_split(halfedge_handle(edge, 0));
}


int StarShapifyableExpansionCone::add_visible_neighbors_to_unremoved_candidates_to_candidates_list(std::vector<VertexHandle>& candidate_vertices){

    SS_DEBUG(" --------------------------------------------------"<<std::endl);
    SS_DEBUG(" - adding visible neighbors to unremoved candidates to candidates list..."<<std::endl);
    SS_DEBUG(" - initial size: "<<candidate_vertices.size()<<std::endl);

    auto already_candidate_prop = request_vertex_property<bool>();
    for(auto candidate_v: candidate_vertices){
        already_candidate_prop[candidate_v] = true;
    }

    std::vector<VertexHandle> new_candidates;

    bool added_at_least_one(false);
    for(auto candidate_v: candidate_vertices){
        SS_DEBUG(" -- checking candidate "<<candidate_v<<std::endl);
        if(!collapsed_prop_[candidate_v]){
            for(auto out_he: outgoing_halfedges(candidate_v)){
                auto neighbor = to_vertex_handle(out_he);
                if(!collapsed_prop_[neighbor] &&
                        neighbor != witness_vertex_ &&
                        !tip_vertices_prop_[neighbor]){
                    new_candidates.push_back(neighbor);
                    already_candidate_prop[neighbor] = true;
                    added_at_least_one = true;
                    SS_DEBUG(" ---> added candidate neighbor "<<neighbor<<std::endl);
                }
            }
        }
    }

    for(auto v: new_candidates){
        candidate_vertices.push_back(v);
    }

    SS_DEBUG(" ...done, final size: "<<candidate_vertices.size()<<std::endl);
    SS_DEBUG(" --------------------------------------------------"<<std::endl);

    return !added_at_least_one;
}



bool StarShapifyableExpansionCone::tips_1_ring_neighborhood_is_visible_from_witness_vertex(bool print_details) const{


    if(print_details){
        SS_DEBUG(" =========================="<<std::endl);
        SS_DEBUG(" VISIBILITY CHECK:"<<std::endl);
    }
    bool all_visible(true);
    for(auto cone_tip: cone_tip_vertices_){
        //SS_DEBUG(" - cone tip "<<cone_tip<<std::endl);
        for(auto out_he: outgoing_halfedges(cone_tip)){
            auto neighbor = to_vertex_handle(out_he);
            //SS_DEBUG(" -- neighbor "<<neighbor<<std::endl);
            if(!tip_vertices_prop_[neighbor] && neighbor != witness_vertex_){
                bool visible = is_visible_from(witness_vertex_, neighbor);

                all_visible &= visible;
                if(print_details){
                    SS_DEBUG(" - "<<neighbor<<": "<<visible<<std::endl);
                }
            }
        }
    }

    if(print_details){
        SS_DEBUG(" =========================="<<std::endl);
    }

    return all_visible;
}

void StarShapifyableExpansionCone::check_for_timeout() const{

    auto current_time = std::chrono::high_resolution_clock::now();
    int current_duration_s = (float)std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time_).count() /1000;

    if(current_duration_s > max_allocated_time_s_){
        SS_DEBUG(" --> star-shapification timeout"<<std::endl);
        throw TimeOutException();
    }

}

}
