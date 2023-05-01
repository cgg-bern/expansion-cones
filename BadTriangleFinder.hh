#pragma once

#include "CommonMeshDefinitions.hh"


class BadTriangleFinder{

public:

    static std::vector<OpenMesh::FaceHandle> find_flipped_triangles(TriMesh& mesh);


    /* Assumes all points to be in the XY-plane (i.e. zero-Z)*/
    static double signed_triangle_area(TriMesh& mesh,
                                       const OpenMesh::FaceHandle& tri);

private:


};

