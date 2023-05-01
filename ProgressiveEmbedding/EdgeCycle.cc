#include "EdgeCycle.hh"

using namespace OpenVolumeMesh;

EdgeCycle::EdgeCycle(){}

EdgeCycle::EdgeCycle(OpenVolumeMesh:: EdgeHandle first,
                     OpenVolumeMesh:: EdgeHandle second,
                     OpenVolumeMesh:: EdgeHandle third){

    int min_idx = std::min(first.idx(), std::min(second.idx(), third.idx()));
    int max_idx = std::max(first.idx(), std::max(second.idx(), third.idx()));
    int middle_idx(0);

    if(first.idx() > min_idx && first.idx() < max_idx){
        middle_idx = first.idx();
    }else if(second.idx() > min_idx && second.idx() < max_idx){
        middle_idx = second.idx();
    }else{
        middle_idx = third.idx();
    }

    push_back(EdgeHandle(min_idx));
    push_back(EdgeHandle(middle_idx));
    push_back(EdgeHandle(max_idx));
}


std::set<VertexHandle> EdgeCycle::toVertexSet(const TetrahedralMesh& mesh) const{

    std::set<VertexHandle> cycle_vertices;
    for(auto e: (*this)){
        cycle_vertices.insert(mesh.edge(e).from_vertex());
        cycle_vertices.insert(mesh.edge(e).to_vertex());
    }
    return cycle_vertices;
}

void EdgeCycle::print(bool endline) const{
    std::cout<<" {"<<(*this)[0]<<" -> "<<(*this)[1]<<" -> "<<(*this)[2]<<"} ";
    if(endline) std::cout<<std::endl;
}



void EdgeCycle::printVertices(const TetrahedralMesh& mesh, bool endline) const{
    std::set<VertexHandle> vertices;

    for(auto e: (*this)){
        vertices.insert(mesh.edge(e).from_vertex());
        vertices.insert(mesh.edge(e).to_vertex());
    }

    std::cout<<" ( ";
    for(auto v: vertices){
        std::cout<<v<<" ";
    }
    std::cout<<")";
    if(endline) std::cout<<std::endl;
}


void EdgeCycle::print(const TetrahedralMesh& mesh,
                      bool  endline) const{
    print(false);
    std::cout<<" -> ";
    printVertices(mesh, endline);

}
