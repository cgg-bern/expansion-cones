#include "TopologicalLink.hh"


using namespace OpenVolumeMesh;

TopologicalFaceSet link(const TetrahedralMesh& mesh,
                               const VertexHandle& vertex){

    VertexSet vertices;
    EdgeSet edges;
    FaceSet faces;

    for(auto v = mesh.vv_iter(vertex); v.valid(); v++){
        vertices.insert(*v);
    }

    for(auto c = mesh.vc_iter(vertex); c.valid(); c++){
        //std::cout<<" - checking neighbor cell "<<*c<<std::endl;

        for(auto hf = mesh.chf_iter(*c); hf.valid(); hf++){
            //std::cout<<" -- checking halfface "<<*hf<<std::endl;
            bool found(false);
            for(auto v = mesh.hfv_iter(*hf); v.valid(); v++){
                //std::cout<<" --- checking vertex "<<*v<<std::endl;
                if(*v == vertex){
                    found = true;
                }
            }
            if(!found){
                //std::cout<<" --> face "<<mesh.face_handle(*hf)<<" is opposite"<<std::endl;
                faces.insert(mesh.face_handle(*hf));
                for(auto e = mesh.hfe_iter(*hf); e.valid(); e++){
                   // std::cout<<" ----> inserting edge "<<*e<<std::endl;
                    edges.insert(*e);
                    //vertices.insert(mesh.edge(*e).from_vertex());
                    //vertices.insert(mesh.edge(*e).to_vertex());
                }
            }
        }
    }

    for(auto f = mesh.vf_iter(vertex); f.valid(); f++){
        for(auto e = mesh.fe_iter(*f); e.valid(); e++){
            bool found(false);
            if(mesh.edge(*e).to_vertex() == vertex ||
                    mesh.edge(*e).from_vertex() == vertex){
                found = true;
            }
            if(!found){
                edges.insert(*e);
            }
        }
    }

    return {vertices, edges, faces, {}};
}




TopologicalFaceSet link(const TetrahedralMesh& mesh,
                             const EdgeHandle& edge){


    VertexSet vertices;
    EdgeSet edges;

    auto from_vertex = mesh.edge(edge).from_vertex();
    auto to_vertex = mesh.edge(edge).to_vertex();

    //std::cout<<" building edge link for edge "<<edge<<" : "<<from_vertex<<" -> "<<to_vertex<<std::endl;

    //std::cout<<" cell count : "<<mesh.n_cells()<<std::endl;

    for(auto c = mesh.ec_iter(edge); c.valid(); c++){
        //std::cout<<" -- checking cell "<<*c<<std::endl;

        for(auto e = mesh.ce_iter(*c); e.valid(); e++){

            auto other_from_vertex = mesh.edge(*e).from_vertex();
            auto other_to_vertex = mesh.edge(*e).to_vertex();

            //std::cout<<" ---- checking edge "<<*e<<" : "<<other_from_vertex<<" -> "<<other_to_vertex<<std::endl;

            if(from_vertex != other_from_vertex &&
                    from_vertex != other_to_vertex &&
                    to_vertex != other_from_vertex &&
                    to_vertex != other_to_vertex){

                //std::cout<<" ----> found link edge "<<std::endl;

                edges.insert(*e);
                vertices.insert(other_from_vertex);
                vertices.insert(other_to_vertex);
            }
        }
    }


    for(auto f = mesh.ef_iter(edge); f.valid(); f++){
        for(auto v = mesh.fv_iter(*f); v++;){
            if(*v != mesh.edge(edge).to_vertex() &&
                    *v != mesh.edge(edge).from_vertex()){
                vertices.insert(*v);
            }
        }
    }


    return {vertices, edges, {}, {}};

}



TopologicalFaceSet link_outsiders(const TetrahedralMesh& mesh,
                                  const OpenVolumeMesh::EdgeHandle& edge){

    auto edge_link = link(mesh, edge);
    auto from_vertex_link = link(mesh, mesh.edge(edge).from_vertex());
    auto to_vertex_link = link(mesh, mesh.edge(edge).to_vertex());

    return from_vertex_link.intersection(to_vertex_link).subtract(edge_link);
}



bool link_condition(const TetrahedralMesh& mesh,
                    const EdgeHandle& edge){

    auto from_vertex = mesh.edge(edge).from_vertex();
    auto to_vertex = mesh.edge(edge).to_vertex();

    auto vertex_link_intersection = link(mesh, from_vertex).intersection(link(mesh, to_vertex));

    if(vertex_link_intersection.faces().size()){
        return false;
    }

    return vertex_link_intersection == link(mesh, edge);
}




bool link_condition(const TetrahedralMesh& mesh,
                    const OpenVolumeMesh::HalfEdgeHandle& heh){

    return link_condition(mesh, mesh.edge_handle(heh));
}
