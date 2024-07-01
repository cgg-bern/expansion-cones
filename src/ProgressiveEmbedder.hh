#pragma once

#include "Expander.hh"
#include "TetMapper.hh"
#include "SphereMapper.hh"
#include "LightWeightJsonExporter.hh"


namespace OpenVolumeMesh{

typedef enum {TET_ONLY=1, STIFF_TET, TET_TO_SPHERE, TET_TO_RANDOM_STAR_SHAPE} BOUNDARY_MAPPING_METHOD;



class ProgressiveEmbedder{

public:


    /** return value:
     * -3 initial interior mapping sufficient
     * -2 expansion error
     * -1 pre-processing/input error
     *  0 expansion succeded
     *  1 shrinkage failed
     *  2 expansion failed
     */
    static int shrinkAndExpand(TetrahedralMesh& domain_mesh,
                               TetrahedralMesh& codomain_mesh,
                               const std::string& mesh_name,
                               const std::string& output_file_path,
                               int debug_expander);



    static int map_to_unit_tet(TetrahedralMesh& mesh);

    //NOTE: "stiff" in the sense that one of its face will be a single triangle
    static int map_to_stiff_unit_tet(TetrahedralMesh& mesh);

    static int map_to_unit_ball_using_tet(TetrahedralMesh& mesh);

    static int map_to_random_star_shape_using_tet(TetrahedralMesh& mesh);

    static int mesh_contains_flipped_tets_or_degenerate_boundary_triangles(TetrahedralMesh& mesh);




    //interior mappings

    static int interior_tutte_mapping(TetrahedralMesh& mesh);

    static void interior_uniform_smoothing(TetrahedralMesh& mesh,
                                           const double delta_eps = 1e-7,
                                           const int max_iterations = 1e6);



    bool split_boundary_connecting_edges_in_submesh(TetrahedralMesh& mesh,
                                                    ExpansionCone& submesh);

    /*void set_vertex_exactly(const VertexHandle& vh,
                    const VertexPosition& pos);

    VertexPosition exact_vertex(const VertexHandle& vh) const;*/

private:

    ProgressiveEmbedder(TetrahedralMesh& domain_mesh,
                        TetrahedralMesh& codomain_mesh,
                        const std::string& mesh_name,
                        const std::string& output_file_path);

    int shrink_and_expand(int debug_expander);


    int count_interior_vertices() const;

    int count_boundary_vertices() const;

    int count_interior_edges() const;

    int count_boundary_edges() const;


    bool is_interior(const EdgeHandle& eh) const;


    //------ JSON export stuff
    using JsonExporter = LightWeightJsonExporter<std::stringstream>;

    /**
     * @brief setup_json_exporter writes basic mesh stuff:
     *          * mesh name
     *          * mesh elements (V, E, F, C)
     *          * #interior vertices
     *          * #interior edges
     *          * #boudnary edges
     *          * initial #flipped tets
     *          * vertex valence histogram
     *          * edge valence histogram
     * @param exporter
     */
    void setup_json_exporter(JsonExporter& exporter);


    void export_shrinkage_data_to_json(JsonExporter& exporter,
                                       const float shrinkage_time_s,
                                       const int shrinkage_result);


    void export_expansion_data_to_json(JsonExporter& exporter,
                                       const float expansion_time_s,
                                       const int expansion_result,
                                       const int left_to_expand_count,
                                       const int expansion_iteration_count,
                                       const int star_shapifications_count,
                                       const int cluster_expansions_count,
                                       const int cluster_star_shapifications_count);


    TetrahedralMesh& domain_mesh_;
    TetrahedralMesh& codomain_mesh_;
    const std::string mesh_name_;
    const std::string output_file_path_;


#if 0
    //int initial_flipped_tets_count_;

    std::vector<std::list<VertexHandle>> clusters_;
    VertexPropertyT<int> cluster_prop_;

    std::vector<bool> expanded_cluster_;
    int expanded_cluster_count_ = 0;

    std::deque<EdgeHandle> edges_to_shrink_;
    int edge_queue_reset_count_ = 0;

    //std::queue<EdgeHandle> candidate_queue_;

    //int edges_to_shrink_count_;
    EdgePropertyT<bool> to_shrink_prop_;
    EdgePropertyT<bool> shrunk_prop_;
    int shrunk_edges_count_ = 0;

    int LPs_solved_count_ = 0;
    int back_up_positions_count_ = 0;

    //used to skip them when computing negative volumes in 1-ring
    CellPropertyT<bool> zero_volume_cells_prop_;

    //also true if the cluster itself changed
    VertexPropertyT<bool> changed_neighboring_cluster_at_last_iteration_prop_;
    EdgePropertyT<double> total_negative_volume_in_1_ring_prop_;
    EdgePropertyT<bool> computed_1_ring_volume_at_least_once_prop_;
    EdgePropertyT<VertexPosition> best_position_prop_;

    VertexPropertyT<VertexPosition> vertex_position_prop_;


    bool handling_flipped_tets_ = false;

    float avg_cluster_merge_duration_s_ = 0;
    float avg_cc_edges_update_duration_s_ = 0;
    float avg_edges_list_reset_duration_s_ = 0;

    int skipped_flipped_tet_edges_count_ = 0;

    const int n_initial_interior_edges_;
#endif

    //VolumePriorityQueue edge_selection_priority_queue_;
};



}
