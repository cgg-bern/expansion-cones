#include "TopoHelper.hh"

#include <Eigen/Sparse>
#include <Eigen/SVD>


namespace OpenVolumeMesh{




CellHandle TopoHelper::cell_exists(const TetrahedralMesh& mesh,
                                   const std::vector<VertexHandle>& cell_vertices){
    return TopoHelper(mesh).cell_exists(cell_vertices);
}



/** \return HalfEdgeHandle(-1) if not found */
HalfEdgeHandle TopoHelper::halfedge_exists(const TetrahedralMesh& mesh,
                                           VertexHandle from_vertex,
                                           VertexHandle to_vertex,
                                           bool look_for_opposite_halfedge_too){
    return TopoHelper(mesh).halfedge_exists(from_vertex, to_vertex, look_for_opposite_halfedge_too);
}


bool TopoHelper::face_contains_vertex(const TetrahedralMesh& mesh,
                                      const VertexHandle& vertex,
                                      const FaceHandle& face){
    auto hf_vertices = mesh.get_halfface_vertices(mesh.halfface_handle(face, 0));

    return  hf_vertices[0] == vertex ||
            hf_vertices[1] == vertex ||
            hf_vertices[2] == vertex;
}


bool TopoHelper::cell_contains_vertex(const TetrahedralMesh& mesh,
                                      const VertexHandle& vertex,
                                      const CellHandle& cell){
    auto c_vertices = mesh.get_cell_vertices(cell);


    return  c_vertices[0] == vertex ||
            c_vertices[1] == vertex ||
            c_vertices[2] == vertex ||
            c_vertices[3] == vertex;
}

//------------ MEMBER HELPERS ------------//


TopoHelper::TopoHelper(const TetrahedralMesh& mesh) : mesh_(mesh) {}



CellHandle TopoHelper::cell_exists(const std::vector<VertexHandle>& cell_vertices) const{
    if(cell_vertices.size() == 4){
        for(auto vc_it = mesh_.vc_iter(cell_vertices[0]); vc_it.valid(); vc_it++){
            int found_vertices_count(0);
            for(auto cv_it = mesh_.cv_iter(*vc_it); cv_it.valid(); cv_it++){
                if(std::find(cell_vertices.begin(), cell_vertices.end(), *cv_it) != cell_vertices.end()){
                    found_vertices_count++;
                }
            }

            if(found_vertices_count == 4){
                return *vc_it;
            }
        }
    }
    return CellHandle(-1);
}



OpenVolumeMesh::HalfEdgeHandle TopoHelper::halfedge_exists(OpenVolumeMesh::VertexHandle from_vertex,
                                                           OpenVolumeMesh::VertexHandle to_vertex,
                                                           bool look_for_opposite_halfedge_too) const{

    HalfEdgeHandle heh(-1);

    for(auto out_he: mesh_.outgoing_halfedges(from_vertex)){
        if(mesh_.to_vertex_handle(out_he) == to_vertex){
            return out_he;
        }
    }
    if(look_for_opposite_halfedge_too){
        for(auto out_he: mesh_.outgoing_halfedges(to_vertex)){
            if(mesh_.to_vertex_handle(out_he) == from_vertex){
                return out_he;
            }
        }
    }

    return heh;
}




void TopoHelper::printOutHalfedgeInfo(HalfEdgeHandle edge) const{
    auto halfedge = mesh_.halfedge(edge);


    std::cout<<" HalfEdge "<<edge<<": "<<halfedge<<std::endl;
    std::cout<<" --           valence: "<<edgeValence(mesh_.edge_handle(edge))<<std::endl;
    std::cout<<" --          boundary: "<<mesh_.is_boundary(edge)<<std::endl;
    std::cout<<" -- boundary vertices: ("<<mesh_.is_boundary(halfedge.from_vertex())<<", "<<mesh_.is_boundary(halfedge.to_vertex())<<")"<<std::endl;
    std::cout<<" --    link condition: "<<link_condition(mesh_, mesh_.edge_handle(edge))<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
}



void TopoHelper::printOutMesh(){

    std::cout<<" ====================== "<<std::endl;
    std::cout<<" VERTICES : "<<std::endl;
    for(const auto& v: mesh_.vertices()){
        std::cout<<v<<" : "<<mesh_.vertex(v)<<std::endl;
    }
    std::cout<<" ====================== "<<std::endl;
    std::cout<<" CELLS : "<<std::endl;
    for(const auto& c: mesh_.cells()){
        std::cout<<c<<" : ";
        for(auto cv_it = mesh_.cv_iter(c); cv_it.valid(); cv_it++){
            std::cout<<*cv_it<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<" ====================== "<<std::endl;

}

std::set<std::pair<std::set<VertexHandle>, bool>> TopoHelper::findNonCellTets(const TetrahedralMesh& mesh){
    auto mesh_copy = mesh;
    return findNonCellTets(mesh_copy);
}



std::set<std::pair<std::set<VertexHandle>, bool>> TopoHelper::findNonCellTets(TetrahedralMesh& mesh_){

    std::set<std::pair<std::set<VertexHandle>, bool>> non_cell_tets;

    for(auto v: mesh_.vertices()){

        auto neighbor = mesh_.request_vertex_property<bool>("neighbor");

        for(auto vv_it = mesh_.vv_iter(v); vv_it.valid(); vv_it++){
            neighbor[*vv_it] = true;
        }

        for(auto vf_it = mesh_.vf_iter(v); vf_it.valid(); vf_it++){

            std::set<VertexHandle> tet_vertices;
            tet_vertices.insert(v);

            //find opposite edge to v
            EdgeHandle opposite_edge(-1);

            for(auto fe_it = mesh_.fe_iter(*vf_it); fe_it.valid(); fe_it++){
                if(mesh_.edge(*fe_it).from_vertex() != v &&
                        mesh_.edge(*fe_it).to_vertex() != v){
                    tet_vertices.insert(mesh_.edge(*fe_it).from_vertex());
                    tet_vertices.insert(mesh_.edge(*fe_it).to_vertex());

                    opposite_edge = *fe_it;
                }
            }

            //iterate through the opposite edge's faces to try to find another neighbor vertex
            for(auto ef_it = mesh_.ef_iter(opposite_edge); ef_it.valid(); ef_it++){
                std::set<VertexHandle> partial_cell = tet_vertices;

                //find the opposite vertex
                for(auto fv_it = mesh_.fv_iter(*ef_it); fv_it.valid(); fv_it++){
                    if(mesh_.edge(opposite_edge).from_vertex() != *fv_it &&
                            mesh_.edge(opposite_edge).to_vertex() != *fv_it &&
                            neighbor[*fv_it]){

                        partial_cell.insert(*fv_it);

                        bool found_cell(false);

                        //check if it's a cell or not
                        for(auto vc_it = mesh_.vc_iter(v); vc_it.valid(); vc_it++){
                            int found_vertices(0);
                            for(auto cv_it = mesh_.cv_iter(*vc_it); cv_it.valid(); cv_it++){
                                if(partial_cell.find(*cv_it) != partial_cell.end()){
                                    found_vertices++;
                                }
                            }
                            if(found_vertices == 4){
                                found_cell = true;
                            }
                        }

                        if(!found_cell){

                            if(partial_cell.size() != 4){
                                std::cout<<" ERROR - tet does not contain 4 vertices"<<std::endl;
                                return {};
                            }
                            //check that at least two vertices are non-boundary
                            int boundary_count(0);
                            for(auto v: partial_cell){
                                boundary_count += mesh_.is_boundary(v);
                            }
                            if(boundary_count <=2){

                                bool all_edge_non_collapsible(true);
                                std::set<EdgeHandle> tet_edges;
                                for(auto v: partial_cell){
                                    for(auto out_he: mesh_.outgoing_halfedges(v)){
                                        auto edge = mesh_.edge_handle(out_he);

                                        if(partial_cell.find(mesh_.to_vertex_handle(out_he)) != partial_cell.end()){
                                            tet_edges.insert(edge);
                                            all_edge_non_collapsible &= !link_condition(mesh_, edge);
                                        }
                                    }
                                }

                                if(tet_edges.size() != 6){
                                    std::cout<<" ERROR - tet does not contain 6 edges"<<std::endl;
                                    return {};
                                }

                                std::set<FaceHandle> tet_faces;
                                for(auto e: tet_edges){
                                    for(auto ef_it = mesh_.ef_iter(e); ef_it.valid(); ef_it++){
                                        for(auto fv_it = mesh_.fv_iter(*ef_it); fv_it.valid(); fv_it++){
                                            if(*fv_it != mesh_.edge(e).from_vertex() &&
                                                    *fv_it != mesh_.edge(e).to_vertex() &&
                                                    partial_cell.find(*fv_it) != partial_cell.end()){
                                                tet_faces.insert(*ef_it);
                                            }
                                        }
                                    }
                                }

                                if(tet_faces.size() == 4){
                                    non_cell_tets.insert({partial_cell, all_edge_non_collapsible});
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return non_cell_tets;
}




bool TopoHelper::singleConnectedComponent(const TetrahedralMesh& mesh){
    auto mesh_copy = mesh;
    return singleConnectedComponent(mesh_copy);
}

bool TopoHelper::singleConnectedComponent(TetrahedralMesh&  mesh){

    auto visited_prop = mesh.request_vertex_property<bool>("visited");

    std::set<VertexHandle> to_visit;
    to_visit.insert(*mesh.vertices().first);

    while(!to_visit.empty()){
        auto next = *to_visit.begin();
        to_visit.erase(to_visit.begin());

        visited_prop[next] = true;

        for(auto out_he_it = mesh.outgoing_halfedges(next).first; out_he_it.valid(); out_he_it++){
            auto to_vertex = mesh.to_vertex_handle(*out_he_it);
            if(!visited_prop[to_vertex]){
                to_visit.insert(to_vertex);
            }
        }
    }

    bool all_visited(true);

    for(auto v: mesh.vertices()){
        if(!visited_prop[v]){
            all_visited = false;
        }
    }


    return all_visited;
}


bool TopoHelper::containsVoid(const TetrahedralMesh& mesh){
    auto mesh_copy = mesh;
    return containsVoid(mesh_copy);
}


bool TopoHelper::containsVoid(TetrahedralMesh&  mesh){

    std::vector<FaceHandle> to_visit;

    for(auto f: mesh.faces()){
        if(mesh.is_boundary(f)){
            to_visit.push_back(f);
            break;
        }
    }

    if(!to_visit.size()){
        std::cerr<<" couldn't find a boundary face, which is very weird"<<std::endl;
        return false;
    }

    auto visited = mesh.request_face_property<bool>("visited");

    int iteration_count(0);

    while(to_visit.size()){
        //std::cout<<" stack size at iteration "<<iteration_count<<": "<<to_visit.size()<<std::endl;

        auto current_face = to_visit.back();
        to_visit.pop_back();

        if(!visited[current_face]){

            visited[current_face] = true;

            for(auto fe_it = mesh.fe_iter(current_face); fe_it.valid(); fe_it++){

                for(auto ef_it = mesh.ef_iter(*fe_it); ef_it.valid(); ef_it++){
                    if(mesh.is_boundary(*ef_it) && !visited[*ef_it]){
                        to_visit.push_back(*ef_it);
                    }
                }
            }
        }

        iteration_count++;
    }


    for(auto f: mesh.faces()){
        if(mesh.is_boundary(f) && !visited[f]){
            std::cout<<" found a face not visited"<<std::endl;
            return true;
        }
    }



    return false;
}


std::vector<VertexHandle> TopoHelper::nonManifoldBoundaryVertices(const TetrahedralMesh& mesh){
    auto mesh_copy = mesh;
    return nonManifoldBoundaryVertices(mesh_copy);
}

std::vector<VertexHandle> TopoHelper::nonManifoldBoundaryVertices(TetrahedralMesh& mesh){

    std::vector<VertexHandle> non_manifold_boundary_vertices;

    for(auto v: mesh.vertices()){
        if(mesh.is_boundary(v)){

            if(!manifoldVertex(mesh, v)){
                non_manifold_boundary_vertices.push_back(v);
            }

        }
    }

    return non_manifold_boundary_vertices;
}





bool TopoHelper::manifoldVertex(const TetrahedralMesh& mesh,
                                const VertexHandle& vertex){
    auto mesh_copy = mesh;
    return manifoldVertex(mesh_copy, vertex);
}

bool TopoHelper::manifoldVertex(TetrahedralMesh& mesh,
                                const VertexHandle& vertex){

    if(!mesh.is_boundary(vertex)){
        return true;
    }

    //find an initial boundary face
    FaceHandle initial_boundary_face(-1);
    for(auto vf_it = mesh.vf_iter(vertex); vf_it.valid(); vf_it++){
        if(mesh.is_boundary(*vf_it)){
            initial_boundary_face = *vf_it;
            break;
        }
    }

    if(initial_boundary_face.idx() == -1){
        std::cerr<<" couldn't find boundary face adjacent to boundary vertex. returning false"<<std::endl;
        return false;
    }

    //std::cout<<" ----------------------------"<<std::endl;
    //std::cout<<" checking vertex "<<vertex<<std::endl;

    //std::cout<<" -- initial boundary face : "<<initial_boundary_face<<": "<<mesh.get_halfface_vertices(mesh.halfface_handle(initial_boundary_face,0))<<std::endl;


    auto visited_edge = mesh.request_edge_property<bool>("visited");
    auto visited_face = mesh.request_face_property<bool>("visited");

    EdgeHandle next_edge(-1);

    //visit all edges of the initial face
    for(auto fe_it = mesh.fe_iter(initial_boundary_face); fe_it.valid(); fe_it++){
        if(mesh.edge(*fe_it).to_vertex() == vertex ||
                mesh.edge(*fe_it).from_vertex() == vertex){
            visited_edge[*fe_it] = true;
            //std::cout<<" ---- visited edge "<<mesh.edge(*fe_it)<<" and marked it as current edge"<<std::endl;
            next_edge = *fe_it;
        }
    }
    visited_face[initial_boundary_face] = true;

    while(next_edge.idx() != -1){
        //std::cout<<" --------- "<<std::endl;

        auto current_edge = next_edge;
        next_edge = EdgeHandle(-1);

        //circulate through boundary faces around edge to find the next one
        for(auto ef_it = mesh.ef_iter(current_edge); ef_it.valid(); ef_it++){
            //std::cout<<" -- checking face "<<(*ef_it)<<": "<<mesh.get_halfface_vertices(mesh.halfface_handle(*ef_it,0))<<std::endl;

            if(mesh.is_boundary(*ef_it) && !visited_face[*ef_it]){
                //found a neighboring face to both current_edge and vertex
                //so visit the next edge and use it as next current_edge
                visited_face[*ef_it] = true;

                //std::cout<<" ---- neighbor face: "<<*ef_it<<": "<<mesh.get_halfface_vertices(mesh.halfface_handle(*ef_it,0))<<std::endl;

                for(auto fe_it = mesh.fe_iter(*ef_it); fe_it.valid(); fe_it++){
                    if(mesh.is_boundary(*fe_it)&&
                            (mesh.edge(*fe_it).to_vertex() == vertex ||
                             mesh.edge(*fe_it).from_vertex() == vertex) &&
                            !visited_edge[*fe_it]){
                        visited_edge[*fe_it] = true;
                        next_edge = *fe_it;
                       //std::cout<<" ---- next neighbor edge: "<<*fe_it<<": "<<mesh.edge(*fe_it)<<std::endl;
                       break;
                    }
                }
                break;
            }
            visited_face[*ef_it] = true;
        }
    }

    //and finally we check that all boundary faces were visited
    bool all_visited(true);
    for(auto vf_it = mesh.vf_iter(vertex); vf_it.valid(); vf_it++){
        if(mesh.is_boundary(*vf_it)){
            all_visited &= visited_face[*vf_it];
        }
    }
    if(!all_visited){
        //std::cout<<" vertex "<<vertex<<" is non-manifold"<<std::endl;
    }

    return all_visited;
}



std::vector<EdgeHandle> TopoHelper::findNonLinkEdges(TetrahedralMesh& mesh){
    std::vector<EdgeHandle> non_link_edges;

    for(auto e: mesh.edges()){
        if(!link_condition(mesh, e)){
            non_link_edges.push_back(e);
        }
    }
    return non_link_edges;
}


void TopoHelper::connectivity_matrix_condition_number(TetrahedralMesh& mesh, double& kappa){


    using Matrix = Eigen::MatrixXf;
    Matrix A((int)mesh.n_vertices(), (int)mesh.n_vertices());

    for(auto v: mesh.vertices()){
        int val(0);
        for(auto vv_it = mesh.vv_iter(v); vv_it.valid(); vv_it++){
            A(v.idx(), vv_it->idx()) = 1;
            val++;
        }
        A(v.idx(), v.idx()) = val;
    }


   // A.setFromTriplets(triplets.begin(), triplets.end());


    std::cout<<" computing singular values..."<<std::endl;

    Eigen::BDCSVD<Matrix> svd(A);
    auto sigma = svd.singularValues();
    kappa = sigma(0) / sigma(sigma.size()-1);
}



void TopoHelper::compute_valence_gradient(TetrahedralMesh& mesh,
                                          VertexPropertyT<double>& valence_gradient_prop,
                                          double& avg){
    //std::cout<<" ------------------------"<<std::endl;
    //std::cout<<" valence gradients: "<<std::endl;
    avg = 0;
    for(auto v: mesh.vertices()){
        if(mesh.is_deleted(v)){
            continue;
        }
        double variation(0);
        double v_val(mesh.valence(v));
        for(auto vv_it = mesh.vv_iter(v); vv_it.valid(); vv_it++){
            variation += std::abs((double)mesh.valence(*vv_it) - v_val);
        }
        valence_gradient_prop[v] = variation/v_val;
        avg += valence_gradient_prop[v];
        //std::cout<<"  -- vertex "<<v<<" of valence "<<std::setw(4)<<v_val<<" has valence gradient "<<std::setw(10)<<std::setprecision(2)<<valence_gradient_prop[v]<<std::endl;
    }

    avg /= (double)mesh.n_logical_vertices();

    //std::cout<<" ------------------------"<<std::endl;
    //std::cout<<" average = "<<avg<<std::endl;
    //std::cout<<" ------------------------"<<std::endl;

}


void TopoHelper::regularize_valence(TetrahedralMesh& mesh){

    std::cout<<" ------------------------------------"<<std::endl;
    std::cout<<" regularizing interior valence..."<<std::endl;

    const int from_to_val_max_diff(2);
    const int min_diff(2);
    const int max_iterations(10);
    int i(0);
    std::vector<std::pair<EdgeHandle, int>> to_split;

    int initial_n_edges(mesh.n_edges());

    do{
        std::cout<<" ------------------------------------"<<std::endl;

        to_split.clear();

        for(auto eh: mesh.edges()){

            if(mesh.is_boundary(eh) || eh.idx() > initial_n_edges){
                continue;
            }

            auto from_v = mesh.edge(eh).from_vertex();
            auto to_v   = mesh.edge(eh).to_vertex();

            int from_v_val = mesh.valence(from_v);
            int to_v_val   = mesh.valence(to_v);
            int from_to_min_val = std::min(from_v_val, to_v_val);

            if(std::abs(from_v_val - to_v_val) <= from_to_val_max_diff){
                //std::cout<<" - edge "<<mesh.edge(eh)<<" with min valence "<<from_to_min_val<<" and valence diff = "<<std::abs(from_v_val - to_v_val)<<std::endl;

                int max_val(0);
                int min_val(std::numeric_limits<int>::max());

                auto heh = mesh.halfedge_handle(eh, 0);

                for(auto ehf_it = mesh.ehf_iter(eh); ehf_it.valid(); ehf_it++){
                    auto op_v = mesh.to_vertex_handle(mesh.next_halfedge_in_halfface(heh, *ehf_it));

                    int val = mesh.valence(op_v);
                    max_val = std::max(max_val, val);

                    //std::cout<<"  -- neighbor vertex "<<op_v<<" has valence "<<val<<std::endl;
                }

                if(max_val < from_to_min_val){

                    int diff = std::abs(max_val - from_to_min_val);


                    //std::cout<<"    --> from-to min val to max neighbor val diff = "<<diff<<std::endl;

                    if(diff > min_diff){
                        //std::cout<<" - edge "<<mesh.edge(eh)<<" with min valence "<<from_to_min_val<<" and valence diff = "<<std::abs(from_v_val - to_v_val)<<std::endl;
                        //std::cout<<"   ---> greater than min diff, adding for split"<<std::endl;
                        to_split.push_back({eh, diff});
                    }
                }
            }
        }

        std::cout<<" - found "<<to_split.size()<<" edges to split at iteration "<<(i+1)<<std::endl;


        for(const auto& pair: to_split){

            VertexHandle from_v = mesh.edge(pair.first).from_vertex();
            VertexHandle to_v   = mesh.edge(pair.first).to_vertex();
            std::cout<<" splitting edge "<<mesh.edge(pair.first)<<" "<<pair.second<<" times"<<std::endl;
            for(int j(0); j<pair.second; j++){
                auto heh = mesh.halfedge(from_v, to_v);
                auto mid_v = mesh.split_edge(heh);
                from_v = mid_v;
            }
        }

        i++;

    }while(i < max_iterations &&
           !to_split.empty());


    std::cout<<" ...done"<<std::endl;
    std::cout<<" ------------------------"<<std::endl;
}


bool TopoHelper::noDoubleEdges(const TetrahedralMesh& mesh){
    auto mesh_copy = mesh;
    return noDoubleEdges(mesh_copy);
}


bool TopoHelper::noDoubleEdges(TetrahedralMesh& mesh){

    for(auto v: mesh.vertices()){
        auto visited = mesh.request_vertex_property<HalfEdgeHandle>("visited through", HalfEdgeHandle(-1));

        for(auto out_he: mesh.outgoing_halfedges(v)){
            auto visited_he = visited[mesh.to_vertex_handle(out_he)];
            if(visited_he.idx() == -1){
                visited[mesh.to_vertex_handle(out_he)] = out_he;
            }else{
                std::cout<<" FOUND DOUBLE EDGE : "<<std::endl;
                std::cout<<"    "<<visited_he<<": "<<mesh.halfedge(visited_he)<<std::endl;
                std::cout<<"    "<<out_he<<": "<<mesh.halfedge(out_he)<<std::endl;
                return false;
            }
        }
    }

    return true;
}





bool TopoHelper::interiorVertexLinksAreTwoManifold(TetrahedralMesh& mesh){

    for(auto v: mesh.vertices()){
        std::cout<<" checking vertex "<<v<<std::endl;
        if(mesh.is_boundary(v)){
            continue;
        }
        auto vertex_link_mesh = mesh;
        auto vertex_link = link(mesh, v);

        for(auto v2: vertex_link_mesh.vertices()){
            if(vertex_link.vertices().find(v2) == vertex_link.vertices().end()){
                vertex_link_mesh.delete_vertex(v2);
            }
        }
        for(auto e: vertex_link_mesh.edges()){
            if(vertex_link.edges().find(e) == vertex_link.edges().end()){
                vertex_link_mesh.delete_edge(e);
            }
        }
        for(auto f: vertex_link_mesh.faces()){
            if(vertex_link.faces().find(f) == vertex_link.faces().end()){
                vertex_link_mesh.delete_face(f);
            }
        }

        //vertex_link_mesh.delete_vertex(v);

        for(auto e: vertex_link_mesh.edges()){
            if(vertex_link_mesh.valence(e) != 2){
                std::cout<<" FOUND VERTEX WITH NON-MANIFOLD LINK"<<std::endl;
                std::cout<<" edge "<<e<<": "<<vertex_link_mesh.edge(e)<<" is not valence 2 "<<std::endl;
                mesh = vertex_link_mesh;
                return false;
            }
        }
    }
    return true;
}


bool TopoHelper::isCollapsible(const OpenVolumeMesh::EdgeHandle& edge) const{
    return isInterior(edge) && link_condition(mesh_, edge);
}


bool TopoHelper::isInterior(const OpenVolumeMesh::EdgeHandle& edge) const{
    return !mesh_.is_boundary(mesh_.edge(edge).from_vertex()) &&
            !mesh_.is_boundary(mesh_.edge(edge).to_vertex());
}


int TopoHelper::uncollapsibleCount() const{
    int count(0);
    for(auto e: mesh_.edges()){
        count += !link_condition(mesh_, e);
    }
    return count;
}

int TopoHelper::valence4InteriorVerticesCount() const{
    int count(0);
    for(auto v: mesh_.vertices()){
        count += !mesh_.is_boundary(v) && mesh_.valence(v) == 4;
    }
    return count;
}


int TopoHelper::valence3InteriorEdgesCount() const{
    int count(0);
    for(auto e: mesh_.edges()){
        count += isInterior(e) && edgeValence(e) == 3;
    }
    return count;
}

int TopoHelper::interiorVerticesCount() const{
    int interior_vertices_count(0);
    //std::cout<<" interior vertices: ";
    for(auto v: mesh_.vertices()){
        interior_vertices_count += !mesh_.is_boundary(v);
        /*if(!mesh_.is_boundary(v)){
            std::cout<<v<<" ";
        }*/
    }
    //std::cout<<std::endl;

    return interior_vertices_count;
}


int TopoHelper::boundaryVerticesCount() const{
    int interior_vertices_count(0);
    for(auto v: mesh_.vertices()){
        interior_vertices_count += mesh_.is_boundary(v);
    }
    return interior_vertices_count;
}


int TopoHelper::edgeValence(OpenVolumeMesh::EdgeHandle edge) const{
    int valence(0);
    for(auto c_it = mesh_.ec_iter(edge); c_it.valid(); c_it++){
        valence++;
    }
    return valence;
}


int TopoHelper::vertexValence(OpenVolumeMesh::VertexHandle vertex) const{
    int valence(0);
    for(auto vc_it = mesh_.vc_iter(vertex); vc_it.valid(); vc_it++){
        valence++;
    }
    return valence;
}



void TopoHelper::printOutMeshInformation() {

    auto bad_tets = BadTetFinder::findBadTets(mesh_);

    auto non_manifold_vertices = nonManifoldBoundaryVertices(mesh_);

    std::cout<<"    vertices: "<<mesh_.n_logical_vertices()<<std::endl;
    std::cout<<"       edges: "<<mesh_.n_logical_edges()<<std::endl;
    std::cout<<"       faces: "<<mesh_.n_logical_faces()<<std::endl;
    std::cout<<"       cells: "<<mesh_.n_logical_cells()<<std::endl;

    std::cout<<" ==============================================================="<<std::endl;
    std::cout<<" - CREATED CUBE MAPPING ENVIRONMENT. INITIAL BAD TETS : "<<std::endl;
    std::cout<<" -----                    degenerates: "<<bad_tets.first.size()<<std::endl;
    std::cout<<" -----                        flipped: "<<bad_tets.second.size()<<std::endl;
    std::cout<<" -----                          genus: "<<mesh_.genus()<<std::endl;
    std::cout<<" -----                      single cc: "<<singleConnectedComponent(mesh_)<<std::endl;
    std::cout<<" -----                  contains void: "<<containsVoid(mesh_)<<std::endl;
    std::cout<<" -----                no double edges: "<<noDoubleEdges(mesh_)<<std::endl;
    //std::cout<<" -----    vertex links are 2-manifold: "<<EdgeCollapser::interiorVertexLinksAreTwoManifold(mesh_)<<std::endl;
    std::cout<<" ----- non-manifold boundary vertices: "<<non_manifold_vertices.size()<<std::endl;
    std::cout<<" -----                  non-cell tets: "<<findNonCellTets(mesh_).size()<<std::endl;
    std::cout<<" -----       valence-3 interior edges: "<<valence3InteriorEdgesCount()<<std::endl;
    std::cout<<" -----    valence-4 interior vertices: "<<valence4InteriorVerticesCount()<<std::endl;
    //std::cout<<" -----            lone non-link edges: "<<edge_collapser.findLoneNonLinkEdges().size()<<std::endl;
    std::cout<<" ==============================================================="<<std::endl;

    std::cout<<"    vertices: "<<mesh_.n_logical_vertices()<<std::endl;
    std::cout<<"       edges: "<<mesh_.n_logical_edges()<<std::endl;
    std::cout<<"       faces: "<<mesh_.n_logical_faces()<<std::endl;
    std::cout<<"       cells: "<<mesh_.n_logical_cells()<<std::endl;


}


}
