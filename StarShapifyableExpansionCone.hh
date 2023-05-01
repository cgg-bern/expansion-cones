#pragma once

#include "ExpansionCone.hh"
#include <VertexPriorityQueue.hh>
#include <deque>


namespace OpenVolumeMesh{


typedef enum {SS_ERROR = -1,
              SS_SUCCESS = 0,
              SS_NO_VERTICES_TO_COLLAPSE,
#warning the next one is deprecated and should eventually be removed
              SS_MIN_VALENCE_GREATER_THAN_TWO,
              SS_EDGE_SPLIT_ERROR,
              SS_RESET_REQUESTED,
              SS_GENERAL_FAILURE} STAR_SHAPIFICATION_RESULT;


class StarShapifyableExpansionCone : public ExpansionCone{

public:


    StarShapifyableExpansionCone(const ExpansionCone& other_cone,
                                 const int max_allocated_time_s);

    STAR_SHAPIFICATION_RESULT is_star_shapifyable();

    STAR_SHAPIFICATION_RESULT star_shapify(bool use_faster_but_unsafe_RS,
                                           bool reduce_precision);

    const SplitList& get_split_list() const{
        return split_list_;
    }

    const VertexPropertyT<std::vector<VertexPosition>> get_mid_vertices_position_prop() const{
        return mid_vertices_positions_prop_;
    }

    const VertexPropertyT<std::vector<VertexHandle>> get_mid_vertices_prop() const{
        return mid_vertices_prop_;
    }

    const VertexPropertyT<VertexHandle> get_original_base_vertex_prop(){
        return original_base_vertex_prop_;
    }


    int get_saved_position_bytes() const{
        return saved_position_bytes_;
    }

    void set_vertex(const VertexHandle& vh,
                    const VertexPosition& pos,
                    bool update_max_size = true);


#warning TODO eventually make private passed this point

    void operator=(const StarShapifyableExpansionCone& other);


    /* helper function for tests */
    void translate_to_have_tips_at_origin();

    void find_unvisible_vertices_from(const VertexHandle& witness_vertex,
                                      std::vector<VertexHandle>& unvisible_vertices) const;

    void gather_vertices_outside_of_1_ring_neighborhood_of_witness_vertex(std::vector<VertexHandle>& non_neighbor_vertices);

    //int uncollapsed_cell_valence(const VertexHandle& vh) const;

    /* might be deprecated */
    int add_visible_neighbors_to_unremoved_candidates_to_candidates_list(std::vector<VertexHandle>& candidate_vertices);

    VertexHandle get_last_secondary_vertex(const VertexHandle& vh) const;

    bool tips_1_ring_neighborhood_is_visible_from_witness_vertex(bool print_details = false) const;

    VertexHandle split_base_edge_and_register_split(const EdgeHandle& edge);
    VertexHandle split_base_edge_and_register_split(const HalfEdgeHandle& edge);

    bool is_not_part_of_witness_vertex_1_ring(const EdgeHandle& eh) const;

private:

    //void reset();


    EdgeHandle find_trim_base_boundary_edge_to_start_split_path() const;

    int update_collapsed_edge_prop(const VertexHandle& removed_vertex,
                                   const int cell_valence,
                                   EdgePropertyT<EdgeHandle>& collapsed_edge_prop,
                                   EdgePropertyT<int>& collapse_depth_prop) const;

    int find_split_path_along_RS_sequence(const EdgeHandle& start,
                                          const EdgePropertyT<EdgeHandle>& collapse_to_prop,
                                          std::deque<EdgeHandle>& path_to_split);

    int find_shortest_dual_path_to_boundary(const EdgeHandle& start,
                                            std::deque<EdgeHandle>& path_to_split);

    int split_base_non_link_edges_not_connected_to_witness();

    int split_base_1_2_ring_edges();

    int split_base_non_extended_link_edges_not_connected_to_witness();

    //int unlock_2_ring_vertices_collapse(std::vector<VertexHandle>& vertices_to_collapse);

    //IMPORTANT NOTE: using shortest path is faster but has no guarantee
    int run_reverse_shelling_sequence(const int max_iteration_count,
                                      bool backtrack_splits_along_shortest_path,
                                      std::vector<VertexHandle>& vertices_to_collapse);

    int run_reverse_shelling_sequence_internal(bool backtrack_splits_along_shortest_path,
                                               std::vector<VertexHandle>& vertices_to_collapse);

    //int run_direct_reverse_shelling_sequence(std::vector<VertexHandle>& vertices_to_collapse);

    int split_non_blocking_base_boundary_edge();

    int run_collapse_sequence(const std::vector<VertexHandle>& vertices_to_collapse);

    int run_contraction_sequence(bool reduce_precision);



    /**
     * IMPORTANT: assumes that vh is on a SINGLE tet (as it should be in this context) */
    FaceHandle opposite_face(const VertexHandle& vh) const;

    //assuming that this is a boundary edge, otherwise there could be two opposite vertices
    VertexHandle opposite_trim_base_vertex(const EdgeHandle& eh) const;

    EdgeHandle common_base_edge(const CellHandle& c1,
                                const CellHandle& c2) const;

    /** IMPORTANT: assumes that the vertex is on exactly TWO tets (as it should be in this context)
     *  NOTE: looks into the trimmed copy */
    EdgeHandle find_target_edge(const VertexHandle& vh) const;


    VertexHandle find_RS_target_vertex(const VertexHandle& vh) const;

    VertexPosition get_position_from_face(const FaceHandle& fh) const;
    VertexPosition get_position_from_edge(const EdgeHandle& eh) const;
    VertexPosition get_position_from_vertices(const std::vector<VertexHandle>& vertices) const;
    VertexPosition get_contraction_position_from_vertices(const VertexHandle& vertex_to_contract,
                                                          const std::vector<VertexHandle>& target_vertices,
                                                          const VertexPropertyT<std::vector<VertexHandle>>& collapsed_to_prop);

    VertexHandle collapse_vertex(const VertexHandle& vertex_to_collapse,
                                 const VertexPosition& new_position);

    VertexHandle split_spoke_and_move_vertex(const VertexHandle& vertex_to_collapse,
                                             const VertexPosition& new_position,
                                             const std::vector<VertexHandle>& target_vertices,
                                             bool reduce_precision = false);

    int contract_core();

    void update_vertices_with_secondaries(std::vector<VertexHandle>& vertices) const;

    void count_base_edges(int& original_edges_count,
                          int& new_edges_count) const;

    /** IMPORTANT: assumes that that vh is already located at its exact precision, safe, position*/
    VertexPosition binary_search_min_precision_position(const VertexHandle& vh,
                                                        const std::vector<VertexHandle>& target_vertices);

    /**
     * checks that
     * 1. there are no flipped tets
     * 2. the cone is still star-shaped
     * 3. that all collapsed-to vertices are still "visible" from the witness.
     *    This is done by checking that the neighbors don't switch "sides" after precision reduction
     *    The argument prop stores the side for each neighbor (i.e. the sign based on the dot product */
    bool is_integrity_maintained(const VertexHandle& vertex_to_contract,
                                 const VertexPropertyT<CGAL::Orientation>& neighbors_side_prop,
                                 const std::vector<VertexHandle>& target_vertices,
                                 const CellPropertyT<bool>& degenerate_tet_prop,
                                 const std::vector<std::pair<std::vector<CGAL_ExactPoint3>, CGAL::Sign>>& additional_constraints,
                                 bool print_debug = false);

    void one_ring_neighborhood_byte_size(const VertexHandle& vertex_to_contract,
                                         int& max_byte_size,
                                         float& average_byte_size,
                                         int& total_byte_size) const;


    //stores the same cone but without the collapsed vertices and
    //the vertices whose spoke has been split
    ExpansionCone trimmed_copy_;

    VertexHandle witness_vertex_;

    VertexPropertyT<bool> collapsed_prop_;

    VertexPropertyT<VertexHandle> secondary_vertex_prop_;

    CellPropertyT<bool> witness_vertex_1_ring_cell_prop_;
    VertexPropertyT<bool> witness_vertex_1_ring_base_vertex_prop_;
    EdgePropertyT<bool> original_edge_prop_;
    VertexPropertyT<std::vector<VertexPosition>> mid_vertices_positions_prop_;
    VertexPropertyT<std::vector<VertexHandle>>   mid_vertices_prop_;
    VertexPropertyT<VertexHandle> original_base_vertex_prop_;

    //second value is the shelling order
    //FacePropertyT<std::pair<VertexHandle, int>> attached_vertices_f_prop_;
    //EdgePropertyT<std::pair<VertexHandle, int>> attached_vertices_e_prop_;

    //each halfedge is tip->base
    VertexPropertyT<HalfEdgeHandle> base_vertex_to_spoke_edge_prop_;
    /*EdgePropertyT<FaceHandle>   base_edge_to_spoke_face_prop_;*/

    //first element is the vertex and the second is the vertices of either
    //the edge, face or tet it's been collapsed to
    std::vector<std::pair<VertexHandle, std::vector<VertexHandle>>> collapse_list_;

    std::vector<VertexHandle> primary_reverse_shelling_sequence_;
    std::vector<VertexHandle> full_reverse_shelling_sequence_;
    int remaining_trimmed_copy_vertices_at_last_RS_iteration_;
    int max_vertex_index_at_last_RS_iterations_ = -1;

    std::vector<VertexHandle> vertices_to_contract_;
    VertexPropertyT<bool> contracted_prop_;
    bool found_deg_tets_ = false;

    SplitList split_list_;

    //diagnostics stuff
    int splits_count_ = 0;
    int backtracking_splits_count_ = 0;


    int new_position_max_precision_ = 0;
    std::vector<int> new_position_precision_record_;
    int saved_position_bytes_ = 0;


    //-------------------------- timeout stuff
    void check_for_timeout() const;

    const int max_allocated_time_s_;
    const std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;


};

struct TimeOutException : public std::exception{
    const char * what () const throw ()
    {
        return "timeout time reached";
    }
};

}

