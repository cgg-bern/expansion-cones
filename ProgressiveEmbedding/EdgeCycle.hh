#ifndef EDGECYCLE_HH
#define EDGECYCLE_HH

#include <vector>


#include "CommonMeshDefinitions.hh"


class EdgeCycle : public std::vector<OpenVolumeMesh::EdgeHandle>{
public:

    //default constructor -> invalid EdgeCycle
    EdgeCycle();

    EdgeCycle(OpenVolumeMesh::EdgeHandle,
              OpenVolumeMesh::EdgeHandle,
              OpenVolumeMesh::EdgeHandle);

    std::set<OpenVolumeMesh::VertexHandle> toVertexSet(const TetrahedralMesh& mesh) const;

    void print(bool endline = true) const;

    void print(const TetrahedralMesh&,
               bool  endline = true) const;

    void printVertices(const TetrahedralMesh&,
                       bool  endline = true) const;


private:
};




#endif // EDGECYCLE_HH
