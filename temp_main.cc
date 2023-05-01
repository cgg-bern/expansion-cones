#include <iostream>
#include "CommonMeshDefinitions.hh"

using namespace OpenVolumeMesh;


int main(){
    
    TetrahedralMesh mesh;
    mesh.add_vertex({0,0,0});
    
    std::cout<<" caca"<<std::endl;
    std::cout<<" mesh: "<<mesh.n_vertices()<<std::endl;
    return 0;
}
