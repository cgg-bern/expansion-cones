#include "TopologicalFaceSet.hh"



TopologicalFaceSet::TopologicalFaceSet(const std::initializer_list<OpenVolumeMesh::VertexHandle>& vertices,
                                       const std::initializer_list<OpenVolumeMesh::EdgeHandle>& edges,
                                       const std::initializer_list<OpenVolumeMesh::FaceHandle>& faces,
                                       const std::initializer_list<OpenVolumeMesh::CellHandle>& cells)
    : vertices_(vertices), edges_(edges), faces_(faces), cells_(cells){}


TopologicalFaceSet::TopologicalFaceSet(const std::set<OpenVolumeMesh::VertexHandle>& vertices,
                                       const std::set<OpenVolumeMesh::EdgeHandle>& edges,
                                       const std::set<OpenVolumeMesh::FaceHandle>& faces,
                                       const std::set<OpenVolumeMesh::CellHandle>& cells)
    : vertices_(vertices), edges_(edges), faces_(faces), cells_(cells){}


bool TopologicalFaceSet::operator==(const TopologicalFaceSet& other) const{

    return vertices_ == other.vertices_ &&
            edges_ == other.edges_ &&
            faces_ == other.faces_ &&
            cells_ == other.cells_;

}


bool TopologicalFaceSet::operator!=(const TopologicalFaceSet& other) const{

    return !((*this) == other);

}


TopologicalFaceSet TopologicalFaceSet::intersection(const TopologicalFaceSet& other) const{
    VertexSet vertices;
    for(auto& v: vertices_){
        if(other.vertices_.find(v) != other.vertices_.end()){
            vertices.insert(v);
        }
    }

    EdgeSet edges;
    for(auto& e: edges_){
        if(other.edges_.find(e) != other.edges_.end()){
            edges.insert(e);
        }
    }

    FaceSet faces;
    for(auto& f: faces_){
        if(other.faces_.find(f) != other.faces_.end()){
            faces.insert(f);
        }
    }

    CellSet cells;
    for(auto& c: cells_){
        if(other.cells_.find(c) != other.cells_.end()){
            cells.insert(c);
        }
    }

    return {vertices, edges, faces, cells};
}


TopologicalFaceSet TopologicalFaceSet::subtract(const TopologicalFaceSet& other) const{

    VertexSet vertices;
    for(auto& v: vertices_){
        if(other.vertices_.find(v) == other.vertices_.end()){
            vertices.insert(v);
        }
    }

    EdgeSet edges;
    for(auto& e: edges_){
        if(other.edges_.find(e) == other.edges_.end()){
            edges.insert(e);
        }
    }

    FaceSet faces;
    for(auto& f: faces_){
        if(other.faces_.find(f) == other.faces_.end()){
            faces.insert(f);
        }
    }

    CellSet cells;
    for(auto& c: cells_){
        if(other.cells_.find(c) == other.cells_.end()){
            cells.insert(c);
        }
    }

    return {vertices, edges, faces, cells};
}


const VertexSet& TopologicalFaceSet::vertices() const{ return vertices_;}

const EdgeSet& TopologicalFaceSet::edges() const{ return edges_; }

const FaceSet& TopologicalFaceSet::faces() const{ return faces_; }

const CellSet& TopologicalFaceSet::cells() const{ return cells_; }


void TopologicalFaceSet::print() const{
    std::cout<<" Edges: {"<<std::endl;
    for(auto e: edges_){
        std::cout<<e<<" ";
    }
    std::cout<<"}"<<std::endl;

    std::cout<<" Vertices: {"<<std::endl;
    for(auto v: vertices_){
        std::cout<<v<<" ";
    }
    std::cout<<"}"<<std::endl;
}
