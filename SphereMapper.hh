#pragma once

#include "CommonMeshDefinitions.hh"
#include "BadTetFinder.hh"
#include "ProgEmbeddingHelpers.hh"

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace OpenVolumeMesh{


class SphereMapper
{
public:


    static bool project_boundary_to_sphere(TetrahedralMesh& mesh,
                                           double sphere_radius = 1.f);

    /* Basically map the interior to the center and checks that all
     * tets are positive */
    static bool is_boundary_star_shaped_and_non_degenerate(TetrahedralMesh mesh);


    static bool center_mesh(TetrahedralMesh& mesh,
                            const TetrahedralMesh::PointT& centroid = {0,0,0});

    static Vec3d boundary_centroid(const TetrahedralMesh& mesh);

    static double encompassing_sphere_radius(TetrahedralMesh& mesh);


    static void perturbate_vertices_positions(TetrahedralMesh& mesh);

    static void scale_mesh_to(TetrahedralMesh& mesh,
                              double new_scale);


    static void move_centroid_and_scale_mesh_boundary_surface(TetrahedralMesh& mesh,
                                                              const TetrahedralMesh::PointT& target_centroid,
                                                              double target_surface_area,
                                                              bool move_and_scale_interior = false);

    static void map_interior_to_centroid(TetrahedralMesh& mesh,
                                         const TetrahedralMesh::PointT& centroid = {0,0,0});

    static double boundary_area(const TetrahedralMesh& mesh);



    static double triangle_face_area(const TetrahedralMesh& mesh,
                                     const FaceHandle& face);

};


}
