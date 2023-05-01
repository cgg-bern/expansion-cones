#pragma once

#include "CommonMeshDefinitions.hh"
#include <random>

//#include <CGAL/GMP

namespace OpenVolumeMesh{

class TetMapper{

public:

    /** NOTE: regular tetrahedron centered at the origin by default
        * */
       static bool map_boundary_to_tet(TetrahedralMesh& mesh,
                                       bool restrict_single_triangle_tet_face = false,
                                       std::vector<Vec3d> corner_positions = {{std::sqrt(8.)/3.,0,-1./3.},
                                                                                   {-std::sqrt(2.)/3., std::sqrt(2./3.), -1./3.},
                                                                                   {-std::sqrt(2.)/3., -std::sqrt(2./3.), -1./3.},
                                                                                   {0,0,1}});
private:

        TetMapper(TetrahedralMesh& mesh,
                  std::vector<Vec3d> corner_positions,
                  int rng_seed);

        bool map_boundary_to_tet(bool restrict_single_triangle_tet_face,
                                 std::vector<VertexHandle>* corner_vertices = nullptr);


        bool embed_and_map_to_tet_edges(bool restrict_single_triangle_tet_face,
                                        std::vector<VertexHandle>* corner_vertices);

        bool find_edge_between_corners(int from_idx, int to_idx);

        bool embed_and_map_to_tet_faces();

        //void find_next_face(std::vector<VertexHandle>& face_vertices);

        bool map_to_tet_faces(const std::vector<VertexHandle>& face_vertices);


        VertexHandle find_next_corner(bool restrict_single_triangle_tet_face);

        bool shortest_geodesic_path(VertexHandle from_v,
                                    VertexHandle to_v,
                                    std::vector<VertexHandle>& path) const;

        bool is_tet_embedded(const VertexHandle& vh) const;

        bool all_boundary_vertices_are_tet_embedded() const;

        TetrahedralMesh& mesh_;

        std::vector<Vec3d> tet_corners_positions_;

        std::mt19937 rng_;

        std::vector<VertexHandle> tet_corners_;

        VertexPropertyT<bool> tet_corner_prop_;
        VertexPropertyT<bool> tet_edge_prop_;
        VertexPropertyT<bool> tet_face_prop_;

};


}
