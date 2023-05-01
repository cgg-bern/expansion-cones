#pragma once


#include "CommonMeshDefinitions.hh"
#include <chrono>
#include <queue>




class TetMeshBoundarySplitter
{
public:


    /* \brief splits:
     * - tets with three (or more) boundary faces
     * - interior edges that connect two boundary vertices
     *
     * \arg split_segment_count see splitInteriorEdgesConnectingBoundaryVertices
     * This makes it possible to then map any feature to a plane without creating zero-volume (planar) tets */
    static void preProcessProblematicRegions(TetrahedralMesh& mesh,
                                             int split_segment_count = 2);

    /* \arg split_segment_count the number of segments those edges should be split into
     * 1 means they're not split, n means those edges are split (n-1) times */
    static void splitInteriorEdgesConnectingBoundaryVertices(TetrahedralMesh& mesh,
                                                             int split_segment_count = 2);

    static void handleDBCIvertices(TetrahedralMesh& mesh);

    static void splitTetsWithMoreThanTwoBoundaryFaces(TetrahedralMesh& mesh);

    static void splitInteriorFacesWithBoundaryOnlyEdges(TetrahedralMesh& mesh);

};

