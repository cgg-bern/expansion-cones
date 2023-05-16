#include "TetMeshBoundarySplitter.hh"

using namespace OpenVolumeMesh;


void TetMeshBoundarySplitter::preProcessProblematicRegions(TetrahedralMesh& mesh,
                                                           int split_segment_count){


    std::cout<<" ============================================================== "<<std::endl;
    std::cout<<" preprocessing... "<<std::endl;
    std::cout<<"    initial stats: "<<std::endl;
    std::cout<<"    vertices: "<<mesh.n_vertices()<<std::endl;
    std::cout<<"       edges: "<<mesh.n_edges()<<std::endl;
    std::cout<<"       faces: "<<mesh.n_faces()<<std::endl;
    std::cout<<"       cells: "<<mesh.n_cells()<<std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();


    splitInteriorFacesWithBoundaryOnlyEdges(mesh);

    splitInteriorEdgesConnectingBoundaryVertices(mesh, split_segment_count);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout<<" FULLY PRE-PROCESSED PROBLEMATIC REGIONS IN "<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
    std::cout<<"    vertices: "<<mesh.n_vertices()<<std::endl;
    std::cout<<"       edges: "<<mesh.n_edges()<<std::endl;
    std::cout<<"       faces: "<<mesh.n_faces()<<std::endl;
    std::cout<<"       cells: "<<mesh.n_cells()<<std::endl;
    std::cout<<" ============================================================== "<<std::endl;


}


    void TetMeshBoundarySplitter::handleDBCIvertices(TetrahedralMesh& mesh){

        std::cout<<" ============================================================== "<<std::endl;
        std::cout<<" splitting edges connecting Double-Boundary-Connected Interior vertices to boundary..."<<std::endl;
        int count(0);

        int n_vertices = mesh.n_vertices();

        for(int i(0); i<n_vertices; i++){

            VertexHandle v(i);

            if(!mesh.is_boundary(v)){
                int boundary_neighbors_count(0);
                for(auto out_he: mesh.outgoing_halfedges(v)){
                    if(boundary_neighbors_count >= 2){
                        break;
                    }
                    auto neighbor = mesh.to_vertex_handle(out_he);

                    if(mesh.is_boundary(neighbor)){
                        boundary_neighbors_count++;
                    }
                }

                if(boundary_neighbors_count >= 2){
                    bool first_BC_edge(true);

                    for(auto out_he: mesh.outgoing_halfedges(v)){

                        auto neighbor = mesh.to_vertex_handle(out_he);

                        if(mesh.is_boundary(neighbor)){
                            if(first_BC_edge){
                                first_BC_edge = false;
                            }else{
                                mesh.split_edge(out_he);
                                count++;
                            }
                        }
                    }
                }
            }
        }

        std::cout<<"...done! Split "<<count<<" edges "<<std::endl;
        std::cout<<" ============================================================== "<<std::endl;
    }


void TetMeshBoundarySplitter::splitInteriorEdgesConnectingBoundaryVertices(TetrahedralMesh& mesh,
                                                                           int split_segment_count){

    if(split_segment_count < 2){
        std::cout<<" nothing to split"<<std::endl;
        return;
    }

    std::cout<<" split segment count = "<<split_segment_count<<std::endl;

    //looking for problematic edges
    int split_count = 0;

    std::vector<HalfEdgeHandle> halfedges_to_split;

    for(auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++){
        auto edge_vertices = mesh.edge_vertices(*e_it);
        if(mesh.is_boundary(edge_vertices[0]) &&
                mesh.is_boundary(edge_vertices[1]) &&
                !mesh.is_boundary(*e_it)){

            auto mid_vertex = mesh.split_edge(*e_it);

            if(split_segment_count > 2){
                auto segment = mesh.halfedge(mid_vertex, edge_vertices[0]);
                //std::cout<<" - split he "<<mesh.edge(*e_it)<<", adding he "<<segment<<std::endl;
                halfedges_to_split.push_back(segment);
            }

            split_count++;
        }
    }
    split_segment_count -= 2;

    std::cout<<" done first "<<split_count<<" splits, remaining segments to create: "<<split_segment_count<<std::endl;

    while(split_segment_count){

        split_segment_count--;

        std::vector<HalfEdgeHandle> next_round_to_split;

        for(auto to_split: halfedges_to_split){
            auto from_vertex = mesh.from_vertex_handle(to_split);

            auto mid_vertex = mesh.split_edge(to_split);
            split_count++;

            auto segment = mesh.halfedge(from_vertex, mid_vertex);

            //std::cout<<" - split he "<<mesh.halfedge(to_split);

            if(split_segment_count){
                //std::cout<<", adding he "<<segment;
                next_round_to_split.push_back(segment);
            }
            //std::cout<<std::endl;

        }

        halfedges_to_split = next_round_to_split;

        std::cout<<" --> done round of splits, remaining segments to create: "<<split_segment_count<<std::endl;
        std::cout<<"                           remaining    edges to  split: "<<halfedges_to_split.size()<<std::endl;
        std::cout<<"                                            split count: "<<split_count<<std::endl;
        std::cout<<" -------------------------------------------------------------------------"<<std::endl;
    }



    mesh.collect_garbage();
    std::cout<<" ---------------------"<<std::endl;
    std::cout<<"    resulting stats after "<<split_count<<" interior edge splits: "<<std::endl;
    std::cout<<"    vertices: "<<mesh.n_vertices()<<std::endl;
    std::cout<<"       edges: "<<mesh.n_edges()<<std::endl;
    std::cout<<"       faces: "<<mesh.n_faces()<<std::endl;

}



void TetMeshBoundarySplitter::splitTetsWithMoreThanTwoBoundaryFaces(TetrahedralMesh& mesh){
    int split_count(0);

    std::cout<<" ---------------------"<<std::endl;
    std::cout<<" looking for tets with more than two boundary faces"<<std::endl;
    std::vector<OpenVolumeMesh::FaceHandle> faces_to_split;


    for(auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); c_it++){

        auto cf_it = mesh.cf_iter(*c_it);
        int boundary_face_count(0);

        while(cf_it.valid()){

            if(mesh.is_boundary(*cf_it)){
                boundary_face_count++;
            }

            cf_it++;
        }

        if(boundary_face_count > 2){

            //std::cout<<"cell "<<*c_it<<" has "<<boundary_face_count<<" boundary faces. "<<std::endl;
            split_count++;

            cf_it = mesh.cf_iter(*c_it);

            while(cf_it.valid()){
                if(!mesh.is_boundary(*cf_it)){
                    faces_to_split.push_back(*cf_it);
                    //std::cout<<"    face "<<*cf_it<<" added to the split list"<<std::endl;
                }
                cf_it++;
            }
        }
    }

    for(auto f_it = faces_to_split.begin(); f_it != faces_to_split.end(); f_it++){
        mesh.split_face(*f_it);
    }

    mesh.collect_garbage();
    std::cout<<"    resulting stats after "<<split_count<<" tet splits: "<<std::endl;
    std::cout<<"    vertices: "<<mesh.n_vertices()<<std::endl;
    std::cout<<"       edges: "<<mesh.n_edges()<<std::endl;
    std::cout<<"       faces: "<<mesh.n_faces()<<std::endl;
    std::cout<<" ---------------------"<<std::endl;



}

void TetMeshBoundarySplitter::splitInteriorFacesWithBoundaryOnlyEdges(TetrahedralMesh& mesh){

    int split_count(0);

    std::vector<OpenVolumeMesh::FaceHandle> faces_to_split;

    std::cout<<" ---------------------"<<std::endl;
    std::cout<<" looking for non boundary faces with boundary only edges..."<<std::endl;

    for(auto f: mesh.faces()){

        if(!mesh.is_boundary(f)){

            int boundary_edges_count(0);

            for(auto fe_it = mesh.fe_iter(f); fe_it.valid(); fe_it++){
                if(mesh.is_boundary(*fe_it)){
                    boundary_edges_count++;
                }
            }

            if(boundary_edges_count == 3){
                faces_to_split.push_back(f);
            }
        }
    }

    std::cout<<" faces to split: "<<faces_to_split.size()<<std::endl;

    for(auto f: faces_to_split){
        if(!mesh.is_deleted(f)){
            mesh.split_face(f);
            split_count++;
        }
    }

    mesh.collect_garbage();
    std::cout<<"    resulting stats after "<<split_count<<" tet splits: "<<std::endl;
    std::cout<<"    vertices: "<<mesh.n_vertices()<<std::endl;
    std::cout<<"       edges: "<<mesh.n_edges()<<std::endl;
    std::cout<<"       faces: "<<mesh.n_faces()<<std::endl;
    std::cout<<" ---------------------"<<std::endl;

}

