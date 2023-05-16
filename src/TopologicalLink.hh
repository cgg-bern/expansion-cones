#pragma once

#include "TopologicalFaceSet.hh"


/* Computes the 'Link' of a vertex in the topological sense */
TopologicalFaceSet link(const TetrahedralMesh& mesh,
                               const OpenVolumeMesh::VertexHandle& vertex);



/* Computes the 'Link' of an edge in the topological sense */
TopologicalFaceSet link(const TetrahedralMesh& mesh,
                             const OpenVolumeMesh::EdgeHandle& edge);


/* Returns the set of topological faces that are in the intersection of
the edge's vertices' links but not in the edge's link */
TopologicalFaceSet link_outsiders(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::EdgeHandle& edge);


/* Checks that the 'link condition' for a particular edge holds.
i.e. If edge goes from vertex a to vertex b, it checks whether

Lk(a) [intersection] Lk(b) = Lk(ab)

is true or not.

This can typically be used to check if this edge can be safely collapsed
without hurting the mesh's topology. */
bool link_condition(const TetrahedralMesh& mesh,
                    const OpenVolumeMesh::EdgeHandle& eh);


bool link_condition(const TetrahedralMesh& mesh,
                    const OpenVolumeMesh::HalfEdgeHandle& heh);
