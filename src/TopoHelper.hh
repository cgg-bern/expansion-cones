#pragma once


#include "BadTetFinder.hh"
#include "TopologicalLink.hh"

namespace OpenVolumeMesh{


class TopoHelper{
public:

    /** ------------ STATIC HELPERS ------------ **/

    static CellHandle cell_exists(const TetrahedralMesh& mesh,
                                  const std::vector<VertexHandle>& cell_vertices);



    /** \return HalfEdgeHandle(-1) if not found */
    static HalfEdgeHandle halfedge_exists(const TetrahedralMesh& mesh,
                                          VertexHandle from_vertex,
                                          VertexHandle to_vertex,
                                          bool look_for_opposite_halfedge_too = false);

    static bool face_contains_vertex(const TetrahedralMesh& mesh,
                                     const VertexHandle& vertex,
                                     const FaceHandle& face);

    static bool cell_contains_vertex(const TetrahedralMesh& mesh,
                                     const VertexHandle& vertex,
                                     const CellHandle& cell);

    /** ------- MESH PROPERTIES ANALYSIS ------- **/

    //WARNING: const ref static members make a copy of the mesh so
    //they can modify the mesh argument to create properties
    static std::set<std::pair<std::set<VertexHandle>, bool>> findNonCellTets(const TetrahedralMesh&);
    static bool singleConnectedComponent(const TetrahedralMesh&);
    static bool containsVoid(const TetrahedralMesh&);
    static std::vector<VertexHandle> nonManifoldBoundaryVertices(const TetrahedralMesh&);
    static bool manifoldVertex(const TetrahedralMesh&,
                               const VertexHandle&);

    static std::set<std::pair<std::set<VertexHandle>, bool>> findNonCellTets(TetrahedralMesh&);
    static bool singleConnectedComponent(TetrahedralMesh&);
    static bool containsVoid(TetrahedralMesh&);
    static std::vector<VertexHandle> nonManifoldBoundaryVertices(TetrahedralMesh&);
    static bool manifoldVertex(TetrahedralMesh&,
                               const VertexHandle&);

    static std::vector<EdgeHandle> findNonLinkEdges(TetrahedralMesh&);

    static void connectivity_matrix_condition_number(TetrahedralMesh& mesh, double& kappa);

    static void compute_valence_gradient(TetrahedralMesh& mesh,
                                         VertexPropertyT<double>& valence_gradient_prop,
                                         double& avg);

    static void regularize_valence(TetrahedralMesh& mesh);

    /** ------- CONSISTENCY CHECKS ------- **/

    static bool noDoubleEdges(const TetrahedralMesh&);
    static bool interiorVertexLinksAreTwoManifold(const TetrahedralMesh&);

    static bool noDoubleEdges(TetrahedralMesh&);
    static bool interiorVertexLinksAreTwoManifold(TetrahedralMesh&);



    /** ------- TEXTUAL ANALYSIS ------- **/

    void printOutMesh();
    void printOutHalfedgeInfo(HalfEdgeHandle) const;
    void printOutMeshInformation();



    //------------ MEMBER HELPERS ------------//

    //those let you use the functions without passing the mesh everytime

    TopoHelper(const TetrahedralMesh& mesh);


    CellHandle cell_exists(const std::vector<VertexHandle>& cell_vertices) const;


    //std::vector<VertexHandle> halfface_vertices(const)


    /** \return HalfEdgeHandle(-1) if not found */
    HalfEdgeHandle halfedge_exists(VertexHandle from_vertex,
                                   VertexHandle to_vertex,
                                   bool look_for_opposite_halfedge_too = false) const;


    /** ------- HELPERS ------- **/

    /* True if the edge satisfies the link condition and both its vertices are non-boundary*/
    bool isCollapsible(const EdgeHandle& edge) const;

    bool isInterior(const EdgeHandle& edge) const;

    int uncollapsibleCount() const;

    int valence4InteriorVerticesCount() const;

    int valence3InteriorEdgesCount() const;

    int interiorVerticesCount() const;

    int boundaryVerticesCount() const;

    int edgeValence(EdgeHandle edge) const;

    int vertexValence(VertexHandle vertex) const;

private:

    const TetrahedralMesh& mesh_;
};





}

