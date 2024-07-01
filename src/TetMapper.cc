#include "TetMapper.hh"

#include "VertexPriorityQueue.hh"

#include <queue>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

#include "BadTetFinder.hh"

namespace OpenVolumeMesh{

bool TetMapper::map_boundary_to_tet(TetrahedralMesh& mesh,
                                    bool restrict_single_triangle_tet_face,
                                    std::vector<Vec3d> corner_positions){



    std::cout<<" ======================================================================"<<std::endl;
    std::cout<<" MAPPING BOUNDARY TO "<<(restrict_single_triangle_tet_face ? "STIFF " : "")<<"TETRAHEDRON"<<std::endl;

    if(corner_positions.size() != 4){
        std::cout<<" ERROR - can't map boundary to tet without the four tet vertices"<<std::endl;
        return false;
    }

    const int max_attempts_count(100);
    bool success(false);
    bool is_flipped(true);
    int i(0);

    auto mesh_copy = mesh;

    //std::vector<VertexHandle> last_tet_corners;
    //NOTE: the flipped tet condition is their to retry if the corners make up a flipped tet
    while((!success || BadTetFinder::meshContainsFlippedTets(mesh)) && i<max_attempts_count ){
        mesh = mesh_copy;
        TetMapper mapper(mesh, corner_positions, i);
        success =  mapper.map_boundary_to_tet(restrict_single_triangle_tet_face);
        i++;
        if(!(i%20)){
            std::cout<<" attempts count = "<<i<<std::endl;
        }
    }
    std::cout<<" done mapping boundary to tet after "<<i<<" attempts. Success: "<<success<<std::endl;

    /*if(BadTetFinder::meshContainsFlippedTets(mesh)){
        std::cout<<" --> looks like the tet boundary is inverted, swapping two corners and recomputing mapping"<<std::endl;
        TetMapper mapper(mesh, corner_positions);
        std::vector<VertexHandle> swapped_tet_corners = {last_tet_corners[0],
                                                        last_tet_corners[1],
                                                        last_tet_corners[3],
                                                        last_tet_corners[2]};
        success = mapper.map_boundary_to_tet(&swapped_tet_corners);
    }*/
    std::cout<<" ======================================================================"<<std::endl;

    return success;

}


TetMapper::TetMapper(TetrahedralMesh& mesh,
                     std::vector<Vec3d> corner_positions,
                     int rng_seed) :
    mesh_(mesh),
    tet_corners_positions_(corner_positions),
    rng_(rng_seed),
    tet_corner_prop_(mesh_.request_vertex_property<bool>()),
    tet_edge_prop_(mesh_.request_vertex_property<bool>()),
    tet_face_prop_(mesh_.request_vertex_property<bool>())
    //tet_edge_edge_prop_(codomain_mesh_.request_edge_property<bool>())
{}


bool TetMapper::map_boundary_to_tet(bool restrict_single_triangle_tet_face,
                                    std::vector<VertexHandle>* corner_vertices){

    if(tet_corners_positions_.size() != 4){
        std::cout<<" ERROR - can't map boundary to tet without the four tet vertices"<<std::endl;
        return false;
    }

    if(corner_vertices){
        std::cout<<" -> using pre-set corners "<<(*corner_vertices)<<std::endl;

        if(corner_vertices->size() != 4){
            std::cout<<" ERROR - only passed "<<corner_vertices->size()<<" corner vertices to tetMapper "<<std::endl;
            return false;
        }
    }


    if(!embed_and_map_to_tet_edges(restrict_single_triangle_tet_face, corner_vertices)){
        //std::cout<<" error while embedding and mapping tet edges"<<std::endl;
        return false;
    }

    if(!embed_and_map_to_tet_faces()){
        //std::cout<<" error while embedding and mapping tet faces"<<std::endl;
        return false;
    }

    if(!all_boundary_vertices_are_tet_embedded()){
        //std::cout<<" ERROR - some boundary vertices are not tet-embedded"<<std::endl;
        return false;
    }


    return true;
}


bool TetMapper::embed_and_map_to_tet_edges(bool restrict_single_triangle_tet_face,
                                           std::vector<VertexHandle>* corner_vertices){

    for(int i(0); i<4;i++){
        VertexHandle corner(-1);
        if(corner_vertices){
            corner = (*corner_vertices)[i];
        }else{
            corner = find_next_corner(restrict_single_triangle_tet_face);
        }

        if(!corner.is_valid()){
            std::cout<<" ERROR - couldn't find corner "<<i<<std::endl;
            return false;
        }
        tet_corner_prop_[corner] = true;
        tet_corners_.push_back(corner);
        mesh_.set_vertex(corner, tet_corners_positions_[i]);

        //std::cout<<" - corner "<<i<<" is "<<corner<<" at "<<codomain_mesh_.vertex(corner)<<", boundary: "<<codomain_mesh_.is_boundary(corner)<<std::endl;

        for(auto j(0); j<i; j++){
            if(!find_edge_between_corners(i, j)){
                //std::cout<<" error while embedding edge"<<std::endl;
                return false;
            }
        }
    }

    return true;
}



bool TetMapper::find_edge_between_corners(int from_idx, int to_idx){

    //std::cout<<" ---------------------------"<<std::endl;
    //std::cout<<" finding edge between corners "<<from_idx<<" and "<<to_idx<<std::endl;
    if(from_idx < 0 || from_idx > 3 || to_idx < 0 || to_idx > 3){
        std::cout<<" ERROR - invalid indices"<<std::endl;
        return false;
    }

    if(from_idx >= (int)tet_corners_.size() || to_idx >= (int)tet_corners_.size()){
        std::cout<<" ERROR - indices out of tet corners range"<<std::endl;
        return false;
    }

    auto from_v = tet_corners_[from_idx];
    auto   to_v = tet_corners_[to_idx];

    std::vector<VertexHandle> edge_vertices;

    shortest_geodesic_path(from_v, to_v, edge_vertices);


    if(!edge_vertices.empty()){
        auto direction = (mesh_.vertex(to_v) - mesh_.vertex(from_v));
        double t = 1. / (double)(edge_vertices.size() + 2);

        //std::cout<<" edge vertices ("<<edge_vertices.size()<<"): "<<edge_vertices<<std::endl;
        for(int i(0); i<(int)edge_vertices.size(); i++){
            mesh_.set_vertex(edge_vertices[i], mesh_.vertex(from_v) + (i+1) * t * direction);
            tet_edge_prop_[edge_vertices[i]] = true;
        }

    }else if(!mesh_.halfedge(from_v, to_v).is_valid()){
        //std::cout<<" ERROR - tet edge ("<<from_idx<<", "<<to_idx<<") is empty and mesh edge doesn't exist"<<std::endl;
        return false;
    }else{
        //std::cout<<" WARNING - tet edge ("<<from_idx<<", "<<to_idx<<") is empty"<<std::endl;
    }

    //std::cout<<" ---------------------------"<<std::endl;
    return true;
}



VertexHandle TetMapper::find_next_corner(bool restrict_single_triangle_tet_face){


    //std::cout<<" ------- "<<std::endl;

    if(restrict_single_triangle_tet_face){
        if(tet_corners_.size() == 1){

            for(auto vv_it = mesh_.vv_iter(tet_corners_[0]); vv_it.valid(); vv_it++){
                if(mesh_.is_boundary(*vv_it)){
                    return *vv_it;
                }
            }
            return VertexHandle(-1);
        }else if(tet_corners_.size() == 2){

            for(auto vv_it = mesh_.vv_iter(tet_corners_[0]); vv_it.valid(); vv_it++){
                if(mesh_.is_boundary(*vv_it)){
                    for(auto vv_it2 = mesh_.vv_iter(*vv_it); vv_it2.valid(); vv_it2++){
                        if(*vv_it2 == tet_corners_[1]){
                            return *vv_it;
                        }
                    }
                }
            }

            return VertexHandle(-1);
        }
    }

    VertexHandle next_corner(-1);
    int i(0);
    do{
        int index = rng_() % mesh_.n_vertices();
        //std::cout<<" - index = "<<index<<std::endl;
        next_corner = VertexHandle(index);
        //std::cout<<" - checking candidate "<<next_corner<<std::endl;

        i++;
    }while((!mesh_.is_boundary(next_corner) ||
           is_tet_embedded(next_corner)) &&
           i<(int)mesh_.n_vertices());

    return next_corner;

}



bool TetMapper::embed_and_map_to_tet_faces(){


    std::vector<VertexHandle> free_vertices;
    for(auto v: mesh_.vertices()){
        if(mesh_.is_boundary(v) && !is_tet_embedded(v)){
            free_vertices.push_back(v);
        }
    }

    if(!map_to_tet_faces(free_vertices)){
        std::cout<<" error while mapping face to tet"<<std::endl;
        return false;
    }

    for(auto v: mesh_.vertices()){
        if(mesh_.is_boundary(v) && !is_tet_embedded(v)){
            tet_face_prop_[v] = true;
        }
    }

#if 0
    for(int i(0); i<4; i++){
        std::vector<VertexHandle> face_vertices;
        find_next_face(face_vertices);

        if(!face_vertices.empty()){

            std::cout<<" - found "<<face_vertices.size()<<" vertices for face "<<i<<": ";
            for(auto v: face_vertices){
                std::cout<<" "<<v;
            }
            std::cout<<std::endl;

            if(!map_to_tet_faces(face_vertices)){
                std::cout<<" error while mapping face to tet"<<std::endl;
                return false;
            }

            for(auto v: face_vertices){
                tet_face_prop_[v] = i+1;
                //codomain_mesh_.set_vertex(v, {0.25,0.25,0.25});
            }

        }else{
            std::cout<<" WARNING - empty tet face"<<std::endl;
        }
    }
#endif

    return true;
}


#if 0
void TetMapper::find_next_face(std::vector<VertexHandle>& face_vertices){

    std::cout<<" ------------------------------------------------------------------------ "<<std::endl;
    std::cout<<" ------------------------------------------------------------------------ "<<std::endl;

    //first, find a boundary vertex that is not embedded on the tet yet
    VertexHandle seed(-1);
    for(auto v: codomain_mesh_.vertices()){
        if(codomain_mesh_.is_boundary(v) &&
                !is_tet_embedded(v)){
            std::cout<<" - found face seed "<<v<<std::endl;
            seed = v;
            break;
        }
    }

    if(!seed.is_valid()){
        std::cout<<" WARNING - couldn't find any seed for next face"<<std::endl;
        return;
    }


    std::queue<FaceHandle> to_visit;
    auto visited_prop = codomain_mesh_.request_face_property<bool>();
    for(auto vf_it = codomain_mesh_.vf_iter(seed); vf_it.valid(); vf_it++){
        if(codomain_mesh_.is_boundary(*vf_it)){
            visited_prop[*vf_it] = true;
            to_visit.push(*vf_it);
        }
    }

    auto added_prop = codomain_mesh_.request_vertex_property<bool>();

    int i = 0;
    while(!to_visit.empty() && i < codomain_mesh_.n_faces()){

        std::cout<<" ---------------------"<<std::endl;
        auto f = to_visit.front();
        to_visit.pop();

        std::cout<<" - visiting face "<<f<<": "<<codomain_mesh_.get_halfface_vertices(codomain_mesh_.halfface_handle(f, 0))<<std::endl;

        //adding the current face's vertices
        for(auto fv_it = codomain_mesh_.fv_iter(f); fv_it.valid(); fv_it++){
            std::cout<<" -- checking vertex "<<(*fv_it)<<std::endl;
            if(codomain_mesh_.is_boundary(*fv_it) &&
                    !added_prop[*fv_it] &&
                    !is_tet_embedded(*fv_it)){
                face_vertices.push_back(*fv_it);
                added_prop[*fv_it] = true;
                std::cout<<" ---> added vertex "<<std::endl;
            }
        }

        //std::cout<<" ------------------------------------------------------------------------ "<<std::endl;
        //std::cout<<" - incident face to vertex "<<v<<std::endl;

        //and setting up the neighboring faces for visitation
        for(auto fe_it = codomain_mesh_.fe_iter(f); fe_it.valid(); fe_it++){
            auto from_v = codomain_mesh_.edge(*fe_it).from_vertex();
            auto   to_v = codomain_mesh_.edge(*fe_it).to_vertex();

            std::cout<<" -- checking edge "<<codomain_mesh_.edge(*fe_it)<<", tet-embedded: "<<is_tet_embedded(from_v)<<"/"<<is_tet_embedded(to_v)<<std::endl;

            /*if(!tet_edge_edge_prop_[*fe_it]){
                for(auto ef_it = codomain_mesh_.ef_iter(*fe_it); ef_it.valid(); ef_it++){
                    std::cout<<" --- checking neighbor face "<<codomain_mesh_.get_halfface_vertices(codomain_mesh_.halfface_handle(*ef_it, 0))<<std::endl;
                    if(codomain_mesh_.is_boundary(*ef_it) &&
                            *ef_it != f &&
                            !visited_prop[*ef_it]){

                        std::cout<<" ----> is boundary and not the visited face, adding it for visitation"<<std::endl;
                        visited_prop[*ef_it] = true;
                        to_visit.push(*ef_it);
                    }
                }
            }*/
        }

        i++;
    }

    std::cout<<" face vertices from faces: ";
    for(auto v: face_vertices){
        std::cout<<" "<<v;
    }
    std::cout<<std::endl;

}
#endif


bool TetMapper::map_to_tet_faces(const std::vector<VertexHandle>& free_vertices){


    using namespace Eigen;

    //stores the free vertices' index in the vector above
    auto variable_index_property = mesh_. template request_vertex_property<int>("variable index");

    int v_idx(0);
    for(auto v: free_vertices){
            variable_index_property[v] = v_idx;
            v_idx++;
    }

    if(free_vertices.empty()){
        std::cout<<" CANNOT MAP FACE WITH NO FREE VERTICES"<<std::endl;
        return false;
    }


    //setting up K and b
    const size_t variable_count(free_vertices.size());
    SparseMatrix<double> K(variable_count, variable_count);
    VectorXd bx(variable_count);
    VectorXd by(variable_count);
    VectorXd bz(variable_count);

    //std::cout<<"free vertices count : "<<variable_count<<std::endl;

    for(size_t i(0); i<variable_count; i++){

        //std::cout<<"checking vertex "<<free_vertices[i]<<std::endl;

        int degree(0);
        double coeff_x(0);
        double coeff_y(0);
        double coeff_z(0);

        auto ve_it = mesh_.ve_iter(free_vertices[i]);
        while(ve_it.valid()){
            if(mesh_.is_boundary(*ve_it)){
                auto edge_vertices = mesh_.edge_vertices(*ve_it);
                VertexHandle neighbor = edge_vertices[0] == free_vertices[i] ? edge_vertices[1] : edge_vertices[0];

                //std::cout<<" ---- checking connected vertex "<<neighbor<<std::endl;

                if(is_tet_embedded(neighbor)/*cube_edge_property_[neighbor] == surrounding_edge_indices[0] ||
                        cube_edge_property_[neighbor] == surrounding_edge_indices[1]  ||
                        cube_edge_property_[neighbor] == surrounding_edge_indices[2]  ||
                        cube_edge_property_[neighbor] == surrounding_edge_indices[3]  ||
                        cube_corner_property_[neighbor] >= 0*/){

                    //std::cout<<" added to coeff sum : "<<codomain_mesh_.vertex(neighbor)[free_coordinate_x]<<", "<<codomain_mesh_.vertex(neighbor)[free_coordinate_y]<<std::endl;

                    coeff_x += mesh_.vertex(neighbor)[0];
                    coeff_y += mesh_.vertex(neighbor)[1];
                    coeff_z += mesh_.vertex(neighbor)[2];
                    degree++;

                }else /*if(cube_face_property_[neighbor] == face_index)*/{
                    //std::cout<<"  added it at index "<<variable_index_property[neighbor]<<std::endl;
                    K.insert(i, variable_index_property[neighbor]) = -1;
                    degree++;
                }
            }

            ve_it++;
        }
        //std::cout<<"--------"<<std::endl;

        K.insert(i, i) = degree;
        bx[i] = coeff_x;
        by[i] = coeff_y;
        bz[i] = coeff_z;
    }


    /*std::cout<<"K : "<<K<<std::endl;
    std::cout<<"b_x : "<<bx<<std::endl;
    std::cout<<"b_y : "<<by<<std::endl;*/


    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(K);
    if(solver.info()!=Success) {
        std::cout<<"decomposition failed !"<<std::endl;
        // decompositionskype failed
        return false;
    }

    auto x = solver.solve(bx);
    if(solver.info()!=Success) {
        std::cout<<"Solving x failed ! "<<std::endl;
        return false;
    }

    auto y = solver.solve(by);
    if(solver.info()!=Success) {
        std::cout<<"Solving y failed ! "<<std::endl;
        return false;
    }

    auto z = solver.solve(bz);
    if(solver.info()!=Success) {
        std::cout<<"Solving z failed ! "<<std::endl;
        return false;
    }

    /*std::cout<<"x : "<<x<<std::endl;
    std::cout<<"y : "<<y<<std::endl;

    std::cout<<"x size : "<<x.size()<<std::endl;
    std::cout<<"y size : "<<y.size()<<std::endl;*/

    //taking the fixed coordinate of one vertex of one of the surrounding edges
    //double fixed_coordinate_value = codomain_mesh_.vertex(cube_edges_vertices_[surrounding_edge_indices[0]][0])[fixed_coordinate];

    //std::cout<<"fixed coordinate value : "<<fixed_coordinate_value<<std::endl;

    //move face vertices at the right locations
    for(size_t i(0); i<variable_count; i++){
        Vec3d new_pos;
        //std::cout<<" pos = "<<x[i]<<", "<<y[i]<<std::endl;
        new_pos[0] = x[i];
        new_pos[1] = y[i];
        new_pos[2] = z[i];
        mesh_.set_vertex(free_vertices[i], new_pos);
    }

    //std::cout<<" done !"<<std::endl;
    return true;
}


bool TetMapper::shortest_geodesic_path(VertexHandle source,
                                       VertexHandle target,
                                       std::vector<VertexHandle>& path) const{


    using Vertex = VertexHandle;

    path.clear();

    //std::cout<<"looking for shortest path from "<<source<<" to "<<target<<std::endl;

    if(source.idx() < 0 || target.idx() < 0){
        return false;
    }

    //inverting source and target to get the right path when backtracking at the end
    VertexHandle temp = target;
    target = source;
    source = temp;

    auto distance = mesh_.request_vertex_property<double>("distance");
    auto previous = mesh_.request_vertex_property<VertexHandle>("previous vertex");

    VertexPriorityQueue<VertexHandle, double> priority_queue;

    for(auto v: mesh_.vertices()){
        if(mesh_.is_boundary(v)){
            distance[v] = std::numeric_limits<float>::max();
            previous[v] = Vertex(-1);
            priority_queue.push(v, distance[v]);
        }
    }

    distance[source] = 0;
    priority_queue.update_priority(source, 0);

    Vertex current_vertex = source;

    bool found_target = false;

    size_t iteration_count(0);

    while(!priority_queue.empty() &&
          !found_target &&
          iteration_count < std::numeric_limits<int>::max()){

        current_vertex = priority_queue.front();
        priority_queue.pop();

        //std::cout<<"visiting vertex "<<current_vertex<<" at distance "<<distance[current_vertex]<<" and previous vertex "<<previous[current_vertex]<<std::endl;


        //iterate through neighbors by using an edge iterator.
        auto ve_it = mesh_.ve_iter(current_vertex);

        while(ve_it.valid()){

            //this ensures that ALL edges along the path are indeed on the boundary.
            //indeed, it might happen that two vertices are boundary but the
            //edge connecting them is an interior one
            if(mesh_.is_boundary(*ve_it)) {
                auto vertices = mesh_.edge_vertices(*ve_it);


                Vertex neighbor = vertices[0] == current_vertex ? vertices[1] : vertices[0];
                //std::cout<<"checking neighbor "<<neighbor<<std::endl;
                if(!is_tet_embedded(neighbor) ||
                        neighbor == target){

                    auto neighbor_to_current_distance = (mesh_.vertex(current_vertex) -
                                                         mesh_.vertex(neighbor)).norm();
                    //update the neighbor's distance
                    if(distance[neighbor] > distance[current_vertex] + neighbor_to_current_distance){
                        distance[neighbor] = distance[current_vertex] + neighbor_to_current_distance;
                        previous[neighbor] = current_vertex;
                        priority_queue.update_priority(neighbor, distance[neighbor]);

                        //std::cout<<"    vertex previous to "<<neighbor<<" is now "<<previous[neighbor]<<std::endl;
                        //std::cout<<"    updated distance to "<<distance[neighbor]<<std::endl;
                    }
                }
            }
            ve_it++;
        }

        if(current_vertex == target){
            //std::cout<<"found target vertex, stopping Shortest Path"<<std::endl;
            found_target = true;
        }


        //std::cout<<"   vertex queue size after "<<iteration_count<<" iterations : "<<priority_queue.size()<<std::endl;

        iteration_count++;
    }

    if(!found_target){
        std::cout<<" COULD NOT FIND PATH AFTER "<<iteration_count<<" ITERATIONS"<<std::endl;
        return false;
    }

    Vertex step = previous[target];

    //tet_edge_edge_prop_[codomain_mesh_.edge_handle(codomain_mesh_.halfedge(target, step))];

#warning that's where it fails for the furch ball sometimes
    while(step != source){
        if(step.idx() == -1){
            //std::cout<<" BAD STEP IN PATH FROM "<<target<<" TO "<<source<<std::endl;
            return false;
        }

        //tet_edge_edge_prop_[codomain_mesh_.edge_handle(codomain_mesh_.halfedge(step, previous[step]))];

        path.push_back(step);
        step = previous[step];

    }

    //tet_edge_edge_prop_[codomain_mesh_.edge_handle(codomain_mesh_.halfedge(source, path.back()))];

    //std::cout<<" SET PATH."<<std::endl;

    return true;

}

bool TetMapper::is_tet_embedded(const VertexHandle& vh) const{
    return tet_corner_prop_[vh] || tet_edge_prop_[vh] || tet_face_prop_[vh];
}

bool TetMapper::all_boundary_vertices_are_tet_embedded() const{

    for(auto v: mesh_.vertices()){
        if(mesh_.is_boundary(v) && !is_tet_embedded(v)){
            std::cout<<" - vertex "<<v<<" is boundary but not tet-embedded"<<std::endl;
            std::cout<<" - incident faces: "<<std::endl;
            for(auto vf_it = mesh_.vf_iter(v); vf_it.valid(); vf_it++){
                std::cout<<" -- "<<mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0))<<", is boundary: "<<mesh_.is_boundary(*vf_it)<<std::endl;
            }
            return false;
        }
    }
    return true;
}

}
