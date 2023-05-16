#ifndef FACESET_HH
#define FACESET_HH

#include <initializer_list>

#include "CommonMeshDefinitions.hh"


using VertexSet = std::set<OpenVolumeMesh::VertexHandle>;
using EdgeSet   = std::set<OpenVolumeMesh::EdgeHandle>;
using FaceSet   = std::set<OpenVolumeMesh::FaceHandle>;
using CellSet   = std::set<OpenVolumeMesh::CellHandle>;


/** \brief  a set of Faces in the topological sense
 * i.e. vertices, edges, triangles and tetrahedra are all 'faces'
 * This, however, only stores the first three because it's meant to be used to represent the
 * "Link" of a vertex or an edge (see Links.hh)
 */
class TopologicalFaceSet
{
public:

    TopologicalFaceSet(){};

    TopologicalFaceSet(const std::initializer_list<OpenVolumeMesh::VertexHandle>& vertices,
                       const std::initializer_list<OpenVolumeMesh::EdgeHandle>& edges,
                       const std::initializer_list<OpenVolumeMesh::FaceHandle>& faces,
                       const std::initializer_list<OpenVolumeMesh::CellHandle>& cells);


    TopologicalFaceSet(const VertexSet& vertices,
                       const EdgeSet& edges,
                       const FaceSet& faces,
                       const CellSet& cells);

    bool operator==(const TopologicalFaceSet& other) const;
    bool operator!=(const TopologicalFaceSet& other) const;

    const VertexSet& vertices() const;
    const EdgeSet& edges() const;
    const FaceSet& faces() const;
    const CellSet& cells() const;

    /* Computes the per-face-type intersection between two sets.*/
    TopologicalFaceSet intersection(const TopologicalFaceSet& other) const;

    TopologicalFaceSet subtract(const TopologicalFaceSet& other) const;

    void print() const;


private:

    VertexSet vertices_;
    EdgeSet edges_;
    FaceSet faces_;
    CellSet cells_;


};

#endif // FACESET_HH
