#pragma once

#include "EdgeCollapser.hh"
#include "GeometryHelper.hh"

namespace OpenVolumeMesh{

typedef enum{UNSPLIT_SUCCESS, UNSPLIT_GEO_FAILURE, UNSPLIT_TOPO_FAILURE, UNSPLIT_ERROR} UNSPLIT_RESULT;


class Uncollapser{

public:


    Uncollapser(TetrahedralMesh& mesh,
                const SubTriangleMap& sub_triangle_map);

    bool uncollapse_sequence(CollapseSequence collapse_sequence);


    bool uncollapse_sequence(SACSequence split_and_collapse_sequence);

    VertexHandle uncollapse(const EdgeCollapse& collapse);


    UNSPLIT_RESULT unsplit(const EdgeSplit& split);

#if 0
    int experimental_unsplit(const EdgeSplit& split);
#endif

    const SubTriangleMap& get_sub_triangle_map() const;

private:





    /** ------- UNCOLLAPSE HELPERS ------- **/


    VertexHandle uncollapse_helper(const EdgeCollapse&      collapse,
                                   const SubFaceSet&        equatorial_neighborhood,
                                   VertexPropertyT<VertexHandle>& split_result_prop);



    bool insert_spoke_secondary_uncollapses_in_queue(const EdgeCollapse& collapse,
                                                     std::queue<EdgeCollapse>& queue);


    /* This works by extracting the equatorial disc given as a triangle mesh
     * and then mark the vertices that are on the boundary of this tri-mesh.
     * This can thus be used for a complete equatorial disc or any subset of
     * (e.g. an equatorial triangle) */
    bool mark_equatorial_disc_boundary_vertices(const EquatorialDisc& equatorial_disc,
                                                VertexPropertyT<bool>& equatorial_disc_boundary_prop) const;


    EquatorialDisc recover_equatorial_disc(const EdgeCollapse& collapse);


    std::set<VertexHandle> recover_subface_set_vertices(const SubFaceSet& subface_set);



    SubFaceSet find_equatorial_neighborhood(const EdgeCollapse& collapse,
                                            const EquatorialDisc& equatorial_disc) const;





    bool add_cells_from_vertex_and_opposite_halffaces(const VertexHandle& center_vertex,
                                                      const std::vector<std::vector<VertexHandle>>& opposite_halffaces);


    bool extract_sorted_spoke_vertices(const EdgeCollapse& collapse,
                                       const int& spoke_index,
                                       std::vector<VertexHandle>& sorted_spoke_vertices) const;


    EdgeCollapse make_secondary_collapse(const VertexHandle& secondary_to_vertex,
                                         const EdgeCollapse& primary_collapse) const;




    /** ------- UNSPLIT HELPERS ------- **/

    bool is_valid(const EdgeSplit& split) const;

    bool is_unsplittable(const EdgeSplit& split,
                         const HalfEdgeHandle& heh) const;

    bool register_unsplit_in_subtriangle_map(const EdgeSplit& split);


    /** ------- EXPERIMENTAL AND INCOMPLETE UNSPLIT STUFF ------- **/

#if 0
    bool gather_split_center_subhalfedges(const EdgeSplit& split,
                                          std::vector<std::pair<HalfEdgeHandle, int>>& center_subhalfedges);

    bool gather_triangle_halfedges(const SubTriangleFace& triangle,
                                   const int& triangle_index,
                                   std::vector<std::pair<HalfEdgeHandle, int>>& halfedges_on_triangle);

    bool is_unsplittable(const std::pair<HalfEdgeHandle, int>& heh,
                         const EdgeSplit& split) const;


    bool register_unsplit_in_subtriangle_map(const std::pair<HalfEdgeHandle, int>& heh,
                                             const EdgeSplit& split);

#endif

    /** ------- SUBFACE MAP UPDATE ------- **/

    bool update_subface_map(const EdgeCollapse& collapse,
                            const VertexPropertyT<VertexHandle>& uncollapse_result_prop);


    /*update all triangles (i,j,k) s.t.
     * - i = collapse.to_vertex
     * - j is either an equatorial vertex or on the north hemisphere
     * - k is on the north hemisphere
     * by replacing i with collapse.from_vertex */
    bool update_north_hemisphere_triangles(const EdgeCollapse& collapse,
                                           const VertexPropertyT<VertexHandle>& uncollapse_result_prop);


    bool update_equatorial_triangles_and_spokes(const EdgeCollapse& collapse,
                                                const VertexPropertyT<VertexHandle>& uncollapse_result_prop);


    bool update_bypass_triangles(const EdgeCollapse& collapse,
                                 const VertexPropertyT<VertexHandle>& uncollapse_result_prop);

    /**\arg spoke_index the i-th spoke is between the (i-1)-th and i-th triangle
    * i.e. it splits the spoke face connecting the from_vertex to the i-th equatorial vertex
    * NOTE: the orientation of the spoke is (from_vertex, to_vertex, eq_vertex) */
    bool partition_spoke_triangle(const EdgeCollapse& collapse,
                              int spoke_index,
                              const VertexPropertyT<VertexHandle>& split_result_prop);

    bool copy_subfaces_and_replace_vertices_with_uncollapse_result(const SubTriangleFace& source_face,
                                                                   const SubTriangleFace& target_face,
                                                                   const VertexPropertyT<VertexHandle>& uncollapse_result_prop);

    void print_subfaces(std::vector<VertexHandle> face);




    TetrahedralMesh& mesh_;

    SubTriangleMap sub_triangle_map_;

    TopoHelper topo_helper_;
};

}

