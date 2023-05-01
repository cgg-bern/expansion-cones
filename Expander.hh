#pragma once


#include "ExpansionDataLogger.hh"
#include "ConnectedVertexSubsetIterator.hh"
#include "StarShapifyableExpansionCone.hh"
#include "ProgEmbeddingHelpers.hh"

#include <OpenVolumeMesh/FileManager/FileManager.hh>


#define DEFAULT_EXPANSION_TIMEOUT_S (1*12*3600)


namespace OpenVolumeMesh{

    template<typename _A, typename _B>
    class ordered_pair : public std::pair<_A, _B>{

    public:
        ordered_pair() : std::pair<_A, _B>(){}

        ordered_pair(const _A& a, const _B& b) : std::pair<_A, _B>(std::min(a, b), std::max(a,b)){}
    };


    typedef enum{EXPANSION_ERROR=-1,
        EXPANSION_SUCCESS,
        EXPANSION_FAILURE,   //1
        NO_MIN_EXPANSION,    //2
        UNEXP_MIN_CONE,      //3
        MAX_CONE_EQUALS_MIN, //4
        MAX_CONE_SMALLER_THAN_4_CELLS, //5
        FAILURE_CODES_COUNT}
            EXPANSION_RESULT_STATUS;


    struct VertexExpanse{
#warning TODO: rename this to tip vertex
        VertexHandle center_vertex = VertexHandle(-1);

        ExpansionCone expansion_cone;

        VertexPosition new_tip_vertex_pos;

        //probably deprecated
        EXPANSION_RESULT_STATUS result_status = EXPANSION_FAILURE;

        EXPANSION_CHECK_RESULT exp_result = IS_NOT_TOPO_EXPANDABLE;

        STAR_SHAPIFICATION_RESULT ss_result = SS_GENERAL_FAILURE;

    };

    typedef enum{FOUR_STAGE_EXPANSION,
                 FIVE_STAGE_EXPANSION,
                 SIX_STAGE_EXPANSION,
                 EIGHT_STAGE_EXPANSION,
                 EXPANSION_SCHEME_COUNT}
                EXPANSION_SCHEME;

#define DEFAULT_EXPANSION_SCHEME SIX_STAGE_EXPANSION


    class Expander{
    public:


        //static constexpr double default_radius_epsilon = 1e-7;
        //static constexpr double default_distance_epsilon = 1e-6;

        /** return value:
         *  -1 Error
         *   0 Success
         *   1 reached max iterations
         *   2 remaining degenerate tets after full expansion
         *   3 stuck
         *   4 timeout
         *   5 max growth reached
         */
        static int fully_expand_mesh(TetrahedralMesh& mesh,
                                     const TetrahedralMesh& input_mesh,
                                     const std::string& mesh_name,
                                     const std::string& output_file_path,
                                     double& expansion_time_s,
                                     bool silent_mode = false,
                                     EXPANSION_SCHEME expansion_scheme = DEFAULT_EXPANSION_SCHEME,
                                     VertexPropertyT<VertexPosition>* exact_vertices_positions = nullptr);





        /**
         * @brief Expander
         * @param mesh
         * @param mesh_name
         * @param exact_vertices_positions
         * @param expanding_cluster
         * @param allow_flipped_initial_flipped_tets */
        Expander(TetrahedralMesh& mesh,
                 const TetrahedralMesh& input_mesh,
                 const std::string& mesh_name = "",
                 bool silent_mode = false,
                 VertexPropertyT<VertexPosition>* exact_vertices_positions = nullptr,
                 bool expanding_cluster = false,
                 bool allow_flipped_initial_flipped_tets = false,
                 const int max_allocated_time_s = DEFAULT_EXPANSION_TIMEOUT_S);



        [[nodiscard]] bool is_fully_expanded() const;

#warning TODO: use ints
        [[nodiscard]] int to_expand_count() const;

        [[nodiscard]] int left_to_expand_count() const;

        [[nodiscard]] int expanded_count() const;

        int star_shapification_count() const;
        int cluster_expansion_count() const;
        int cluster_star_shapification_count() const;

        [[nodiscard]] const ExpansionDataLogger& data_logger() const;

        //result based on the expanded_prop_ property
        //useful to know if the mesh is already partially expanded
        static bool at_least_one_expanded_vertex(TetrahedralMesh& mesh);

        /* Returns true if single_cc is false or the euler charac. is not 1 */
        int check_main_cluster_topology(bool& single_connected_component,
                                         int& euler_characteristic);


        const SplitList& get_split_list() const{
            return split_list_;
        }

        const VertexPropertyT<VertexPosition>& get_vertex_position_prop() const{
            return vertex_position_prop_;
        }

        VertexPropertyT<bool>& TEMP_expanded_prop_ = expanded_prop_;
        //VertexPropertyT<VertexExpanse>& TEMP_full_cone_at_last_iteration_prop_ = expanse_at_last_iteration_prop_;


        [[nodiscard]] VertexPosition vertex(const VertexHandle& v) const;

        int iteration_count() const;


#warning those should eventually be private


        /** return value:
         *  -1 Error
         *   0 Success
         *   1 reached max iterations
         *   2 remaining degenerate tets after full expansion
         *   3 stuck
         *   4 timeout
         *   5 max growth ratio
         *   6 unknown error */
        int fully_expand_mesh(const std::string& output_file_path,
                              EXPANSION_SCHEME expansion_scheme,
                              double& expansion_time_s);


        int four_stage_full_expansion(const std::string& output_file_path);

        int five_stage_full_expansion(const std::string& output_file_path);

        int six_stage_full_expansion(const std::string& output_file_path);

        int eight_stage_full_expansion(const std::string& output_file_path);


        EXPANSION_RESULT_STATUS expand_all_expandable_vertices(VertexExpanse& expanse,
                                                               bool enable_star_shapification);


        /* \arg stop_at_first_feasible_max_cone passed to find_maximal_expansion_cone */
        int break_best_unexpandable_cone_with_memory(VertexExpanse& expanse,
                                                     bool stop_at_first_feasible_max_cone);

        int brute_force_find_expandable_cluster_and_expand();

        int find_expandable_cluster_and_expand(bool enable_star_shapification,
                                               int max_k = std::numeric_limits<int>::max(),
                                               int max_n_vertices = std::numeric_limits<int>::max());

        //NOTE: only handling edges connected to the primary cluster for now
        int cluster_interface_expansion(const ExpansionCone& cluster_ec);

        /** NOTE: reference to prop so we can add the mapping between the submesh
         *  mid-vertices to the original mesh's */
        int apply_split_list(const SplitList& split_list,
                             VertexPropertyT<VertexHandle>& submesh_to_mesh_prop,
                             bool star_shapifying_cluster = false);

        void reduce_mid_vertices_precision(int split_list_start_index,
                                           int position_byte_size_threshold);

        void smooth_full_interior(int max_iterations,
                                  double delta_epsilon,
                                  int position_byte_size_threshold);

        void smooth_neighborhood_of_unexpanded_vertices(int ring_k,
                                                        int max_iterations,
                                                        double delta_epsilon,
                                                        int position_byte_size_threshold);

        void smooth_mid_vertices_1_ring_neighborhood(int split_list_start_index,
                                                     int max_iterations,
                                                     double delta_epsilon,
                                                     int position_byte_size_threshold);

        void smooth_vertices(const std::vector<VertexHandle>& vertices_to_smooth,
                             int max_iterations,
                             double delta_epsilon,
                             int position_byte_size_threshold);

        int collapse_new_edges(const int split_list_start_index);

        int collapse_new_edges_opposite_to_unexpanded_vertices();

        int collapse_as_many_edges_as_possible(const std::vector<std::pair<VertexHandle, VertexHandle>>& to_collapse);


        //NOTE: cluster_ec only const because we need to add a property for the cluster EC
        int star_shapify_cluster(const ExpansionCone& cluster_ec);


        int gather_full_expansion_cone(const VertexHandle& center_vertex,
                                       ExpansionCone& cone,
                                       bool accumulating_cluster = false) const;

        int split_edges_between_cone_vertices_but_not_part_of_the_cone(ExpansionCone& ss_cone);

        int split_edges_of_tets_with_three_faces_on_cone_base(ExpansionCone& ss_cone);


        void export_cone(ExpansionCone& cone_before_SS,
                         const ExpansionCone& cone_after_SS,
                         int sv_increase_threshold = 200,
                         int sv_increase_ratio_threshold = 10) const;

    private:


        int n_boundary_vertices() const;

        //will expand the first expandable vertex that has unexpansion valence <= target_expansion_valence <= max_expansion_valence
        //OR the expandable vertex with the lowest unexpansion valence if this minimum is greater than target_expansion_valence
        EXPANSION_RESULT_STATUS expand_best_expandable_vertex_with_memory(bool enable_star_shapification,
                                                                          int target_expansion_valence,
                                                                          int max_expansion_valence = std::numeric_limits<int>::max());

        int merge_expansion_cones(const std::vector<VertexHandle>& unexpanded_vertices,
                                  const std::vector<int>& cluster_indices,
                                  ExpansionCone& cluster);





        /********************************************************* CONE BREAKING ***/


        /* starts by finding the max cone so don't use this one if it's already computed */
        int break_expansion_cone_and_expand(const VertexHandle& vertex,
                                            ExpansionCone& maximal_cone,
                                            bool stop_at_first_feasible_max_cone);


        int find_maximal_expansion_cone(const VertexHandle& vertex,
                                        ExpansionCone& maximal_cone,
                                        bool stop_at_first_feasible_max_cone);


        int find_best_minimal_expansion_cone(const VertexHandle& center_vertex,
                                             ExpansionCone& minimal_cone);


        int extract_expansion_cone_around_spoke_edge(const HalfEdgeHandle& spoke,
                                                     ExpansionCone& cone);


        int grow_maximal_expansion_cone(const VertexHandle& center_vertex,
                                        const ExpansionCone& minimal_cone,
                                        ExpansionCone& maximal_cone);

        /** This is needed because sometimes, adding a cell makes it so that
         *  all vertices of an element (face or tet) are inside the cone,
         *  but not the element itself.
         *  This can cause errors when computing the geo-expandability
         *  e.g., with a full cone of three cells and one face
         *  (The top triangle here being "empty")
         *   ____
         *  |\xx/|
         *  | \/ |
         *  | /\ |
         *  |/__\|
         *  If we add the three cells to the cone, it will appear as expandable,
         *  but it's not, because of the lone face.
         *  This lone face should be added */
        int add_side_elements_to_expansion_cone(const VertexHandle& center_vertex,
                                                ExpansionCone& cone);


        int perform_expansion_with_max_cone(const VertexHandle& center_vertex,
                                            const ExpansionCone& maximal_cone);

        int expand_cluster(const ExpansionCone& cluster,
                           const VertexPosition& pos);



        int split_expanded_spokes_out_of_expansion_cone(const VertexHandle& center_vertex,
                                                        const ExpansionCone& cone,
                                                        std::vector<VertexHandle>& steiner_vertices);



        int post_op_fix_degenerate_cells(const std::vector<CellHandle>& degenerate_cells);


        /**************************************************** CONE HELPERS ***/

        CellHandle add_cell_to_cone(const CellHandle& ch,
                                    ExpansionCone& cone) const;

        FaceHandle add_face_to_cone(const FaceHandle& fh,
                                    ExpansionCone& cone) const;

        [[nodiscard]] int expanded_vertices_count(const CellHandle& ch) const;

        [[nodiscard]] int expanded_vertices_count(const FaceHandle& ch) const;


        /********************************************************* HELPERS ***/

        /** Finalize expanse by
         *  1. moving the vertex
         *  2. updating the expanded vertices count (increment)
         *  3. marking the vertex as expanded
         *  4. logging the expanse
         *  5. assembling the VertexExpanse
         *
         *  NOTE: new_vertices_to_expand_count only there for logging process */
        VertexExpanse finalize_expanse(const VertexHandle& center_vertex,
                                       const ExpansionCone& cone, //NOT SURE THIS IS STILL RELEVANT
                                       const VertexPosition& new_vertex_position,
                                       int to_expand_count_increase = 0);

        void move_vertex_and_mark_as_expanded(const VertexHandle& tip_vertex,
                                              const VertexPosition& new_vertex_position);

        VertexPosition find_minimum_precision_position(const VertexHandle& center_vertex,
                                                       const VertexPosition& new_vertex_position);

        VertexPosition binary_search_minimum_precision_position(const VertexHandle& center_vertex,
                                                                ExpansionCone& cone,
                                                                const VertexPosition& new_vertex_position,
                                                                int position_byte_size_threshold);


        void reset_EC_memory();

        bool all_unexpanded_vertices_are_at_the_center(const VertexHandle& center_vertex) const;

        [[nodiscard]] int expansion_valence(const VertexHandle& center_vertex) const;

        [[nodiscard]] int unexpansion_valence(const VertexHandle& center_vertex) const;


        void set_vertex(const VertexHandle& v,
                        const VertexPosition& pos);

        void set_expanded_prop(const VertexHandle& v,
                               const bool status);

        bool expanded_prop(const VertexHandle& v) const;

        VertexHandle split_edge(const EdgeHandle& e);

        VertexHandle split_edge(const HalfEdgeHandle& e);

        //marks halfedges (mid-from) and (mid-to)
        void mark_halfedges_as_new_and_collapsible(const VertexHandle& from_vertex,
                                                   const VertexHandle& mid_vertex,
                                                   const VertexHandle& to_vertex);


        void check_for_timeout() const;

        int remaining_seconds_before_timeout() const;



        /********************************************************* FIELDS ***/

        int iteration_count_ = 0;

        const int initial_n_vertices_;

        //const int initial_non_link_edges_count_;

        const std::chrono::time_point<std::chrono::high_resolution_clock> global_start_time_;
        const int max_allocated_time_s_;

        TetrahedralMesh& mesh_;
        VertexPropertyT<bool> expanded_prop_;
        VertexPropertyT<VertexPosition> vertex_position_prop_;
        int expanded_count_ = 0;
        int to_expand_count_ = 0;
        TopoHelper topo_helper_;
        bool fully_expanded_ = false;
        int cone_break_count_ = 0;
        bool got_stuck_ = false;

        const bool expanding_cluster_;

        VertexPropertyT<VertexPosition> domain_vertex_position_prop_;

        VertexPropertyT<VertexExpanse>             expanse_at_last_iteration_prop_;
        VertexPropertyT<std::vector<VertexHandle>> unexpanded_neighbors_at_expanse_computation_;
        //VertexPropertyT<std::vector<VertexHandle>> unexpanded_neighbors_at_last_cluster_expansion_iteration_;
        VertexPropertyT<bool>                      moved_during_smoothing_prop_;
        //contains the {exp-valence, unexp-valence} pair for each vertex
        VertexPropertyT<std::pair<int,int>>        valence_at_expanse_computation_;

        HalfEdgePropertyT<bool> collapsible_new_halfedge_prop_;

        //VertexPropertyT<std::pair<ExpansionCone, int>> full_cone_at_last_iteration_prop_;


        VertexPropertyT<int> cluster_index_prop_;
        EdgePropertyT<std::pair<int, int>> connecting_cluster_indices_prop_;
        int latest_cluster_index_ = 0;

        ExpansionDataLogger data_logger_;
        int last_single_sv_ec_count_ = 0;
        int last_unexp_valence_ = 0;
        //NOTE: _B as in "Bytes"
        int max_vertex_precision_B_ = 0;
        int total_saved_position_bytes_ = 0;

        //quite the mouthful names, let's try to find a better one
        int single_SV_topo_expandable_but_not_expandable_SV_count_ = 0;
        int single_SV_topo_expandable_and_expandable_SV_count_ = 0;

        double total_simple_expansions_time_s_ = 0;

        double total_smoothing_time_s_ = 0;

        int star_shapification_successes_count_ = 0;        
        double total_ss_time_s_ = 0;

        std::vector<int> cluster_size_count_;
        double total_cluster_exp_time_s_ = 0;

        int cluster_star_shapifications_count_ = 0;
        double total_css_time_s_ = 0;

        int total_collapsed_edges_count_ = 0;
        double total_edge_collapsing_time_s_ = 0;

        //mostly used when a cluster needs some star-shapification and
        //the splits need to be reproduced on the base mesh.
        SplitList split_list_;

         bool silent_mode_;

    };

}

