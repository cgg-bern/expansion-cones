#pragma once

#include "CommonMeshDefinitions.hh"
#include "PrecisionDiagnostics.hh"

//temp
#include <boost/optional/optional_io.hpp>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/intersections.h>


#include <map>
#include <iostream>




namespace OpenVolumeMesh{


    using VertexPosition = ExactVertexPosition;


    typedef enum{
        EXPANDABILITY_CHECK_ERROR = -1,
        IS_EXPANDABLE = 0,               //0
        IS_NOT_TOPO_EXPANDABLE,          //1
        IS_NOT_GEO_EXPANDABLE,           //2
        NO_CELL_IN_EC,                   //3
        EULER_CHARAC_IS_NOT_ONE,         //4
        CONE_BASE_EULER_CHARAC_IS_NOT_ONE, //5
        MULTIPLE_CONNECTED_COMPONENTS,   //6
#warning TODO: check that antenna detection really works
        ANTENNA,                          //7
        SPOKES_NOT_CONTINUOUSLY_EXPANDED, //8
        CLUSTER_TIPS_EULER_CHARACTERISTIC_IS_NOT_ONE, //9
        CLUSTER_TIPS_ARE_NOT_A_SINGLE_CONNECTED_COMPONENT, //10
        CLUSTER_NOT_ALL_FACES_ARE_INCIDENT_TO_AT_LEAST_ONE_TET, //11
        CLUSTER_TIP_IS_NOT_BOUNDARY_2_MANIFOLD, //12
        EXPANSION_CHECK_RESULT_COUNT
    } EXPANSION_CHECK_RESULT;


#warning TODO: rename those for "submesh_"
    struct Split{
        VertexHandle cone_from_vertex;
        VertexHandle cone_to_vertex;

        VertexHandle cone_new_vertex;
        VertexPosition new_vertex_position;

        //only used when a collapse follows the split
        //(cone_new_vertex is collapsed to collapse_to_vertex)
        VertexHandle cone_collapsed_to_vertex = VertexHandle(-1);

        //std::vector<VertexHandle> TEMP_neighbors_at_split_time;

        //keep track of neighbors when split, and check when it diverges when applying split list
    };

    using SplitList = std::vector<Split>;

    std::ostream& operator<<(std::ostream& ss, const Split& split);


    class ExpansionCone : public TetrahedralMesh{
    public:

        /*default lower bound to avoid non-terminating LPs
         * For some unknown reason, setting the lower bound to exactly zero
         * makes the LP solving not terminate...
         * As a hotfix, we'll set this lower bound to something extremely low */
        /*const CGAL::Gmpq chebyshev_radius_lower_bound = ExactType(CGAL::Gmpz(1),
                                                                  CGAL::Gmpz(std::numeric_limits<double>::max()));
*/

        ExpansionCone();

        ExpansionCone(const ExpansionCone& other_cone);

        /**
         * @brief ExpansionCone creates an expansion cone from a full tet mesh
         */
        ExpansionCone(const TetrahedralMesh& mesh);

        void operator=(const ExpansionCone& other_cone);

        EXPANSION_CHECK_RESULT is_expandable(VertexPosition& pos,
                                             bool print_debug = false);


        void set_vertex(const VertexHandle& v,
                        const VertexPosition& pos);

        const std::set<VertexHandle>& cone_tip_vertices() const;

        ///NOTE: only pass CONE vertices, not mesh ones.
        bool is_cone_tip(const VertexHandle& cone_tip) const{
            return tip_vertices_prop_[cone_tip];
        }



        const VertexPropertyT<VertexPosition>& vertex_position_prop() const{
            return vertex_position_prop_;
        }

        int euler_characteristic() const{
            return (int)n_logical_vertices()
                 - (int)n_logical_edges()
                 + (int)n_logical_faces()
                 - (int)n_logical_cells();
        }

        //WARNIN: iterates through all neighboring edges to find a base boundary one
        //so a bit costly
        bool is_base_boundary_vertex(const VertexHandle& vh) const;

        bool is_base_boundary_edge(const EdgeHandle& eh) const;

        bool is_base_edge(const EdgeHandle& eh) const;


#warning handle override stuff

        VertexHandle add_vertex(const VertexHandle& mesh_vertex,
                                const VertexPosition& pos,
                                bool tip_vertex = false);

        virtual CellHandle add_cell(const CellHandle& mesh_cell,
                                    const std::vector<VertexHandle>& original_indices);

        FaceHandle add_face(const std::vector<VertexHandle>& original_indices);

        EdgeHandle add_edge(const VertexHandle& from_vertex,
                            const VertexHandle& to_vertex);

        VertexHandle split_edge(const HalfEdgeHandle& cone_e);

        VertexHandle split_edge(const EdgeHandle& cone_e);

        VertexHandle split_edge(const EdgeHandle& cone_e,
                                const VertexHandle& mesh_mid_v);

        void set_as_tip(const VertexHandle& mesh_vertex);


        const VertexPropertyT<VertexHandle>& cone_to_mesh_v_handle_prop(){
            return cone_to_mesh_v_handle_prop_;
        }


        VertexHandle cone_to_mesh_handle(const VertexHandle& cone_handle) const;
        VertexHandle mesh_to_cone_handle(const VertexHandle& mesh_handle) const;

        CellHandle cone_to_mesh_handle(const CellHandle& cone_handle) const;


        void clear();

        bool contains_flipped_tets() const;
        bool contains_degenerate_tets() const;

        bool base_is_a_single_connected_component();

        int tips_euler_characteristic() const;

        EdgeHandle opposite_edge_on_face(const VertexHandle& vh,
                                         const FaceHandle& fh) const;

        int cell_valence(const VertexHandle& v) const;
        int cell_valence(const EdgeHandle& v) const;


        void remove_tips();

        //replaces all tips with a single one, connected to all base vertices
        //the return value is the new tip
        //WARNING: does NOT call collect_garbage so the deleted entities are still accessible afterwards
        VertexHandle merge_tips();

        VertexHandle find_removable_vertex_with_cell_valence_lower_than_3(const std::vector<VertexHandle>& candidate_vertices,
                                                                          const VertexPropertyT<bool>& to_ignore_prop,
                                                                          const int min_index = -1);




        /** NOTE: mesh is not const because we need a property */
        static int set_up_1_ring_neighborhood_as_expansion_cone(TetrahedralMesh& mesh,
                                                                const VertexPropertyT<VertexPosition>& exact_vertex_position_prop,
                                                                const VertexHandle& vh,
                                                                ExpansionCone& one_ring_EC);


        /** NOTE: mesh is not const because we need a property */
        static int set_up_1_ring_neighborhood_as_expansion_cone(TetrahedralMesh& mesh,
                                                                const VertexPropertyT<VertexPosition>& exact_vertex_position_prop,
                                                                const std::vector<VertexHandle>& vs,
                                                                ExpansionCone& one_ring_EC);


        static int set_up_1_ring_neighborhood_as_expansion_cone(TetrahedralMesh& mesh,
                                                                const std::vector<VertexHandle>& vs,
                                                                ExpansionCone& one_ring_EC);


        /**
         * checks whether a line segment going from one vertex to the other only stays
         * inside of the cone's volume.
         *
         * returns true if there's no other BOUNDARY triangle containing a tip vertex
         * that intersects with the segment between the two vertices passed in argument
         * (other than the ones containing either vertices, which intersect by definition)
         * AND either there's an INTERIOR triangle intersected by the segment
         * or the two vertices are neighbors
         * */
        bool is_visible_from(const VertexHandle& cone_vertex_from,
                             const VertexHandle& cone_vertex_to) const;


        /** computes the north pole as the average of all side faces' normals */
        VertexPosition north_pole() const;

        void project_base_to_sphere();

        /**  projection of the cone base to the plane tangential to -north_pole **/
        void project_base_stereographically(const VertexPosition& north_pole);

        void translate_to_move_tip_to_origin();

#warning TODO eventually make protected passed this point

        EXPANSION_CHECK_RESULT is_topo_expandable(bool print_debug = false);

        EXPANSION_CHECK_RESULT is_geo_expandable(VertexPosition& new_position,
                                                 bool print_debug = false) const;

        bool same_entities_count(const ExpansionCone& other) const;

        int count_interior_vertices() const{
            int count(0);
            for(auto v: vertices()){
                count += !is_boundary(v);
            }
            return count;
        }

        void print_details(int detail_level = 2) const;

        VertexPosition vertex(const VertexHandle& v) const;


        /*int find_max_min_volume_center(VertexPosition& new_position,
                                       int library_to_use,
                                       double eps = 1e-7,
                                       bool print_debug = false);*/

        int find_max_min_volume_center_with_CGAL(VertexPosition& new_position,
                                                 double eps = 1e-7,
                                                 bool print_debug = false) const;


        int total_position_byte_size() const;


        int triangle_normal_and_centroid(const HalfFaceHandle& face,
                                         VertexPosition& normal,
                                         VertexPosition& centroid) const;

        int triangle_normal_and_centroid(const std::vector<VertexHandle>& face_vertices,
                                         VertexPosition& normal,
                                         VertexPosition& centroid) const;

    protected:

        void set_as_tip_internal(const VertexHandle& cone_vertex);


        int find_chebyshev_center(VertexPosition& new_position,
                                  bool print_debug = false) const;


        //NOTE: the center will be on the negative side of each halfface
        //i.e. the halffaces' normal points outward the feasible space
        int find_chebyshev_center_from_halffaces(const std::vector<std::vector<VertexHandle>>& halffaces,
                                                 VertexPosition& new_position,
                                                 bool print_debug = false) const;



        //NOTE: only used when merging multiple ECs. Kind of a hack...
        VertexPropertyT<bool> tip_vertices_prop_;
        std::set<VertexHandle> cone_tip_vertices_;

        VertexPropertyT<VertexHandle> cone_to_mesh_v_handle_prop_;
        CellPropertyT<CellHandle>     cone_to_mesh_c_handle_prop_;
        std::map<VertexHandle, VertexHandle> mesh_to_cone_v_handle_map_;

        VertexPropertyT<VertexPosition> vertex_position_prop_;


    };


    ExactType norm(const VertexPosition& pos);

    ExactType halley_norm(const VertexPosition& pos);

    ExactType l1_norm(const VertexPosition& pos);

    ExactType linf_norm(const VertexPosition& pos);

    ExactType quake_norm(const VertexPosition& pos);

//temporarily passed through double precision for Eigen
    ExactType eigen_norm(const VertexPosition& pos);

    VertexPosition normalize(VertexPosition& pos);

#if ENABLE_EXACT_REPRESENTATION
    VertexPosition vec2vec(const Vec3d& v);

    Vec3d vec2vec(const VertexPosition& v);
#else
    VertexPosition vec2vec(const VertexPosition& v);

#endif


    std::ostream& operator<<(std::ostream& os,
                             const ExpansionCone& cone);



}
