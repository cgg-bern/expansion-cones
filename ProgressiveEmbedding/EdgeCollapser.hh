#pragma once

#include "TopologicalLink.hh"
#include "BadTetFinder.hh"
#include "EdgeCycle.hh"
#include "MeshTopologicalComparator.hh"
#include "SubTriangleMap.hh"
#include "ProgEmbeddingHelpers.hh"

#include <chrono>
#include <random>
#include <algorithm>
#include <queue>
#include <list>


namespace OpenVolumeMesh{

typedef enum{COLLAPSE_ERROR=-1, COLLAPSE_SUCCESS, COLLAPSE_FAILURE, NON_INTERIOR_EDGE_FAILURE, NON_LINK_EDGE_FAILURE, NON_INTERIOR_AND_NON_LINK_EDGE_FAILURE} COLLAPSE_RESULT_STATUS;


/* contains the edge collapsed (represented by its source and target vertices)
 * and the list of equatorial vertices around the original edge.
 * NOTE: to_vertex is the surviving vertex of the collapse */
struct EdgeCollapse{
    VertexHandle from_vertex;
    VertexHandle to_vertex;
    VertexHandle replacement_vertex;

    /** equatorial vertices are in counter-clockwise order around edge (from, to)*/
    std::vector<VertexHandle> equatorial_vertices;

    std::vector<std::vector<VertexHandle>> north_hemispherical_hfs;

    COLLAPSE_RESULT_STATUS result_status = COLLAPSE_SUCCESS;

    bool TEMP_after_split = false;
};

/* contains the edge collapsed (represented by its source and target vertices)
 * and the resulting middle vertex */
struct EdgeSplit{
    VertexHandle from_vertex;
    VertexHandle to_vertex;
    VertexHandle middle_vertex;

    //not sure this is useful
    std::vector<VertexHandle> equatorial_vertices;

    bool is_valid() const {
        return from_vertex.idx() != -1 &&
                middle_vertex.idx() != -1 &&
                to_vertex.idx() != -1;
    }
};


struct SplitAndCollapse{
    EdgeSplit split;
    EdgeCollapse collapse;
};


using SACSequence = std::vector<SplitAndCollapse>;
using CollapseSequence = std::vector<EdgeCollapse>;


class EdgeCollapser
{
public:


    EdgeCollapser(TetrahedralMesh& mesh);


    typedef enum{RANDOM_SELECTION, MAX_VALENCE_SELECTION, MIN_VALENCE_SELECTION, COLLAPSE_HEURISTIC_COUNT} COLLAPSE_SELECTION_HEURISTIC;

    /** ------ STATIC COLLAPSING FUNCTIONS ------- **/

    /* a CollapseFunction takes a tet-mesh as input, collapses it,
     * and returns how many interior vertices remain after the full collapse */
    using CollapseFunction =  int (*)(TetrahedralMesh&);

    using CollapseFunctionWithFileOutput =  int (*)(TetrahedralMesh&, const std::string& mesh_name, const std::string& output_file_path, int option1, int option2);


    static int twoStepCollapse(TetrahedralMesh& mesh);

    /** return value: */
    static int threeStepCollapse(TetrahedralMesh& mesh);


    CollapseSequence getCollapseSequence() const;

    SACSequence getSplitAndCollapseSequence() const;

    const SubTriangleMap& getSubTriangleMap() const;


    /** ------ COLLAPSING ALGORITHMS ------- **/

    /* return value: remaining interior vertices */
    COLLAPSE_RESULT_STATUS twoStepCollapse(bool random_order = false);

    /** \return whether there's a single interior vertex remaining or not */
    COLLAPSE_RESULT_STATUS threeStepCollapse();

    /* return value: success count */
    int twoStepCollapseInRandomOrder(int nb_runs = 100);


    COLLAPSE_RESULT_STATUS collapseUntilPositiveVolume();

    COLLAPSE_RESULT_STATUS collapseAllCollapsibleEdges();

    COLLAPSE_RESULT_STATUS collapseAllCollapsibleBoundaryEdges();

    COLLAPSE_RESULT_STATUS collapseAllCollapsibleEdgesInRandomOrder();


    COLLAPSE_RESULT_STATUS breakTriTets(bool random_order = false);

    COLLAPSE_RESULT_STATUS breakFanTets();

    COLLAPSE_RESULT_STATUS splitAndCollapseCycles(COLLAPSE_SELECTION_HEURISTIC = RANDOM_SELECTION);



    /** ------ EXPERIMENTAL ------- **/


    void splitAllBoundaryEdges();

    // -------- end of experimental

    /** ------- INDIVIDUAL SPLIT AND COLLAPSE ------- **/

    bool split(const EdgeSplit& split);

    bool collapse(const EdgeCollapse& collapse);


    /** ------- COLLAPSING SUB-ALGORITHMS ------- **/

    COLLAPSE_RESULT_STATUS splitAndCollapseNewNeighbor(EdgeHandle eh,
                                                       COLLAPSE_SELECTION_HEURISTIC = RANDOM_SELECTION);

    COLLAPSE_RESULT_STATUS breakCycle(const EdgeCycle& cycle);

    HalfEdgeHandle findNonBoundaryCollapsibleEdgeOnBadTet(const BadTetFinder::BadTetList& bad_tets);

    HalfEdgeHandle findAnyNonBoundaryCollapsibleEdge();

    HalfEdgeHandle findAnyBoundaryCollapsibleEdge();

    std::set<EdgeHandle> findLoneNonLinkEdges();

    std::set<EdgeHandle> findFaceValenceNotEqualCellValenceEdges();

    std::set<EdgeCycle> findNonLinkInteriorEdgeCycles(bool include_face_cycles = true);
    std::set<EdgeCycle> findNonLinkEdgeCyclesWithAtLeastOneInteriorEdge();
    std::set<EdgeCycle> findNonLinkInteriorEdgeCyclesAroundVertex(const VertexHandle& vertex,
                                                                  bool count_special_config_cycles = true);

    EdgeSplit splitEdge(const EdgeHandle& eh);

    EdgeSplit splitEdge(const HalfEdgeHandle& eh);

    EdgeCollapse collapseEdge(const HalfEdgeHandle& eh);

    bool checkThatCyclesAreIndeedNotFaces();

    bool containsNonLinkInteriorCycle();

    void storeLinkConditionAsProperty(bool print_out_status = true);
    bool linkConditionUpToDate() const;

    static void shuffleEdgeVector(std::vector<EdgeHandle>& edge_vector);

    void register_split_in_subtriangle_map(const EdgeSplit& split);

private:


    void setLinkConditionPropStatus(bool up_to_date);

    TetrahedralMesh& mesh_;

    TetrahedralMesh initial_mesh_;

    int total_collapse_count_ = 0;

    EdgePropertyT<bool> link_condition_prop_;
    bool link_condition_is_up_to_date_ = false;

    std::vector<std::pair<int,int>> bad_cycles_count_;

    std::vector<std::vector<int>> split_edge_valence_histogram;

    SACSequence split_and_collapse_sequence_;
    CollapseSequence collapse_sequence_;

    SubTriangleMap sub_triangle_map_;

    TopoHelper topo_helper_;

};


std::ostream& operator<<(std::ostream& os,
                         const EdgeCollapse& collapse);

std::ostream& operator<<(std::ostream& os,
                         const EdgeSplit& split);

std::ostream& operator<<(std::ostream& os,
                         const SplitAndCollapse& sac);




}
