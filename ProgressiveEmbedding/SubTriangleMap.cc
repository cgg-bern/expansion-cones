#include "SubTriangleMap.hh"

namespace OpenVolumeMesh{


//############################# SUB-TRIANGLE FACE ########################
SubTriangleFace::SubTriangleFace(const FaceVertices& vertices){


    if(vertices.size() != 3){
        vertices_ = {VertexHandle(-1), VertexHandle(-1), VertexHandle(-1)};
        std::cerr<<" ERROR - created SubTriangleFace with "<<vertices.size()<<" instead of 3"<<std::endl;
    }else{
        int min_index(0);
        int min_vertex_idx(vertices[0].idx());

        for(size_t i(1); i<3; i++){
            if(vertices[i].idx() < min_vertex_idx){
                min_vertex_idx = vertices[i].idx();
                min_index = i;
            }
        }

        for(size_t i(0); i<3; i++){
            vertices_.push_back(vertices[(min_index + i) % 3]);
        }
    }


    if(vertices[0].idx() == -1 ||
            vertices[0].idx() == -1 ||
            vertices[1].idx() == -1){
        std::cerr<<" ERROR - created SubTriangleFace with an invalid vertex: "<<(*this)<<std::endl;
        vertices_ = {VertexHandle(-1), VertexHandle(-1), VertexHandle(-1)};
        //exit(EXIT_FAILURE);
    }

    if(vertices[0] == vertices[1] ||
            vertices[0] == vertices[2] ||
            vertices[1] == vertices[2]){
        std::cerr<<" ERROR - created SubTriangleFace with twice the same vertex: "<<(*this)<<std::endl;
        vertices_ = {VertexHandle(-1), VertexHandle(-1), VertexHandle(-1)};
        //exit(EXIT_FAILURE);
    }
}



SubTriangleFace::SubTriangleFace(const std::initializer_list<VertexHandle> vertices)
    : SubTriangleFace(FaceVertices(vertices)){}



bool SubTriangleFace::operator<(const SubTriangleFace& other) const{
    if(vertices_[0] < other.vertices_[0]){
        return true;
    }else if(vertices_[0] == other.vertices_[0]){
        if(vertices_[1] < other.vertices_[1]){
            return true;
        }else if(vertices_[1] == other.vertices_[1] &&
                 vertices_[2] < other.vertices_[2]){
            return true;
        }
    }
    return false;
}

bool SubTriangleFace::operator==(const SubTriangleFace& other) const{
    return vertices_[0] == other.vertices_[0] &&
            vertices_[1] == other.vertices_[1] &&
            vertices_[2] == other.vertices_[2];

}

bool SubTriangleFace::operator!=(const SubTriangleFace& other) const{
    return !(*this == other);

}

/* \return VertexHandle(-1) if the index is out of bounds */
VertexHandle SubTriangleFace::operator[](const int& index) const{

    return (index >= 0 && index <3) ? vertices_[index] : VertexHandle(-1);
}

const FaceVertices& SubTriangleFace::get_vertices() const{
    return vertices_;
}



SubTriangleFace SubTriangleFace::opposite_face() const{
    return {vertices_[2], vertices_[1], vertices_[0]};
}


bool SubTriangleFace::replace_vertex(const VertexHandle& source,
                                     const VertexHandle& target){
    bool replaced(false);
    for(auto& v: vertices_){
        if(v == source){
            v = target;
            replaced = true;
        }
    }

    //recreate the subface to make sure the vertices are in the right order
    (*this) = SubTriangleFace(vertices_);

    return replaced;
}


bool SubTriangleFace::contains_edge(const VertexHandle& from_vertex,
                                    const VertexHandle& to_vertex) const{

    int found_vertices(0);
    found_vertices += vertices_[0] == from_vertex || vertices_[0] == to_vertex;
    found_vertices += vertices_[1] == from_vertex || vertices_[1] == to_vertex;
    found_vertices += vertices_[2] == from_vertex || vertices_[2] == to_vertex;
    return found_vertices == 2;
}

bool SubTriangleFace::contains_vertex(const VertexHandle& vertex) const{
    return vertices_[0] == vertex ||
            vertices_[1] == vertex ||
            vertices_[2] == vertex;
}


//############################# SUB-TRIANGLE MAP ########################

bool SubTriangleMap::contains(const SubTriangleFace& face) const{
    return find(face) != end();
}




void SubTriangleMap::insert(const SubTriangleFace& face,
                            const std::vector<SubTriangleFace>& subfaces){


    for(auto subface: subfaces){
        this->insert(face, subface);
    }
}


void SubTriangleMap::insert(const SubTriangleFace& face,
                            const SubTriangleFace& subface){

    if(face == subface){
        std::cerr<<" ERROR - inserting face "<<subface<<" has a subface of itself"<<std::endl;
        return;
    }

    if(!count(face)){
        BaseMap::insert({face, {}});
        BaseMap::insert({face.opposite_face(), {}});
    }

    (*this)[face].insert(subface);
    (*this)[face.opposite_face()].insert(subface.opposite_face());
}

void SubTriangleMap::erase(const SubTriangleFace& face){
    BaseMap::erase(face);
    BaseMap::erase(face.opposite_face());
}


#if 0
bool SubTriangleMap::merge_node_with_super_node(const SubTriangleFace& face,
                                                const SubTriangleFace& super_face){


}
#endif


void SubTriangleMap::extract_all_subvertices(const SubTriangleFace& face,
                                            std::set<VertexHandle>& equatorial_disc_vertices) const{
    //std::cout<<" ------"<<std::endl;
    //std::cout<<" gathering face "<<face<<" subvertices..."<<std::endl;

    std::queue<SubTriangleFace> queue;
    queue.push(face);

    while(!queue.empty()){

        auto current_triangle = queue.front();
        queue.pop();

        //std::cout<<" - exploring triangle "<<current_triangle<<std::endl;
        //insert the current triangle's vertices
        for(auto v: current_triangle.get_vertices()){
            equatorial_disc_vertices.insert(v);
        }

        const auto& current_triangle_subfaces_value = find(current_triangle);
        if(current_triangle_subfaces_value != end()){
            //std::cout<<" - adding its subtriangles to queue: ";
            //and add the subtriangles to the queue
            for(auto subtri: current_triangle_subfaces_value->second){
                queue.push(subtri);
                //std::cout<<subtri<<", ";
            }
            //std::cout<<std::endl;
        }
    }



    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------"<<std::endl;

}




void SubTriangleMap::extract_all_subfaces(const SubTriangleFace& face,
                                          EquatorialDisc& equatorial_disc) const{
    //std::cout<<" ------"<<std::endl;
    //std::cout<<" gathering face "<<face<<" subfaces..."<<std::endl;

    std::queue<SubTriangleFace> queue;
    queue.push(face);

    while(!queue.empty()){

        auto current_face = queue.front();
        queue.pop();

        //std::cout<<" - exploring face "<<current_face<<std::endl;

        const auto& subfaces = find(current_face);
        if(subfaces != end()){

            //std::cout<<" -> node, adding its subfaces to queue: ";
            //and add the subtriangles to the queue
            for(auto subface: subfaces->second){
                if(subface != current_face){
                    queue.push(subface);
                    //std::cout<<subface<<", ";
                }else{
                    std::cerr<<" ERROR - face "<<current_face<<" contains itself as a subface"<<std::endl;
                    equatorial_disc.clear();
                    return;
                }
            }
            //std::cout<<std::endl;
        }else{
            //std::cout<<" --> no sub-face, inserting leaf"<<std::endl;
            equatorial_disc.insert(current_face);
        }
    }

    //std::cout<<" ...done"<<std::endl;
    //std::cout<<" ------"<<std::endl;
}




void SubTriangleMap::extract_subvertices_intersection(const SubTriangleFace& first_face,
                                                      const SubTriangleFace& second_face,
                                                      std::set<VertexHandle>& subvertices_intersection) const{

    //std::cout<<" ------ "<<std::endl;
    //std::cout<<" extracting intersection of subvertices of face "<<first_face<<" and face "<<second_face<<std::endl;

    std::set<VertexHandle> first_face_subvertices;
    std::set<VertexHandle> second_face_subvertices;

    extract_all_subvertices(first_face,
                            first_face_subvertices);

    extract_all_subvertices(second_face,
                            second_face_subvertices);


    for(auto v_first: first_face_subvertices){
        if(second_face_subvertices.find(v_first) != second_face_subvertices.end()){
            subvertices_intersection.insert(v_first);
        }
    }

    //std::cout<<" done"<<std::endl;
    //std::cout<<" ------ "<<std::endl;
}



void SubTriangleMap::replace_vertex_in_all_faces_containing_vertex(const VertexHandle& source,
                                                                   const VertexHandle& target){

    //std::cout<<" ------ "<<std::endl;
    //std::cout<<" replacing vertex "<<source<<" with vertex "<<target<<" in ALL entries containing vertex "<<source<<std::endl;

    std::set<SubTriangleFace> to_delete;

    std::vector<value_type> to_insert;

    for(auto value_pair: *this){

        auto face(value_pair.first);
        //std::cout<<" - checking face "<<face<<std::endl;

        auto& subfaces = value_pair.second;
        std::set<SubTriangleFace> updated_subfaces;

        bool replaced_at_least_one_vertex(false);
        for(auto subface: subfaces){
            if(subface.contains_vertex(source)){
               // std::cout<<" ---> found subface "<<subface<<std::endl;
            }

            replaced_at_least_one_vertex |= subface.replace_vertex(source, target);
            updated_subfaces.insert(subface);
        }

        auto updated_face = face;
        replaced_at_least_one_vertex |= updated_face.replace_vertex(source, target);

        if(replaced_at_least_one_vertex){
            //std::cout<<" --> at least one subface contained vertex "<<source<<
            //           ", preparing updated entry"<<std::endl;

            //prepare old entry deletion
            to_delete.insert(face);

            //and prepare the new entry
            to_insert.push_back({updated_face, updated_subfaces});
        }
    }

    /*std::cout<<" done finding faces. "<<
                   to_delete.size()<<" faces to delete and "<<
                   to_insert.size()<<" faces to insert."<<std::endl<<
                   " Replacing old entries with new ones..."<<std::endl;*/

    //erase the current entry for the updated faces
       for(auto face_to_insert: to_insert){
           erase(face_to_insert.first);
       }



    //delete the old entries
    for(auto face_to_delete: to_delete){
        erase(face_to_delete);
        //std::cout<<" - deleted "<<face_to_delete<<std::endl;
    }

    //and insert the new ones
    for(auto face_to_insert: to_insert){

        insert(face_to_insert.first, {face_to_insert.second.begin(),
                                      face_to_insert.second.end()});

        /*std::cout<<" - inserted "<<face_to_insert.first<<" -> ";
        for(auto subface: face_to_insert.second){
            std::cout<<subface<<" ";
        }
        std::cout<<std::endl;*/
    }

    //std::cout<<" done"<<std::endl;
    //std::cout<<" ------ "<<std::endl;
}



#if 0
SubTriangleMap SubTriangleMap::extract_sub_map(const SubTriangleFace& face) const{

    std::cout<<" ------ "<<std::endl;

    std::cout<<"extracting submap starting at face "<<face<<std::endl;

    SubTriangleMap submap;


    std::queue<SubTriangleFace> queue;
    queue.push(face);

    while(!queue.empty()){

        auto current_face = queue.front();
        queue.pop();

        std::cout<<" - exploring face "<<current_face<<std::endl;

        const auto& subfaces = find(current_face);
        if(subfaces != end()){

            std::cout<<" -> node, adding a map entry"<<std::endl;;
            //and add the subtriangles to the queue
            for(auto subface: subfaces->second){
                if(subface != current_face){
                    queue.push(subface);
                    //std::cout<<subface<<", ";
                }else{
                    std::cerr<<" ERROR - face "<<current_face<<" contains itself as a subface"<<std::endl;
                    submap.clear();
                    return submap;
                }
            }

            submap.insert(subfaces->first, {subfaces->second.begin(), subfaces->second.end()});

            //std::cout<<std::endl;
        }
    }

    std::cout<<" done"<<std::endl;
    std::cout<<" ------ "<<std::endl;

    return submap;
}


bool SubTriangleMap::replace_vertex_in_map_and_update(const VertexHandle& source,
                                                      const VertexHandle& target) const{




    return false;
}
#endif

//############################# EXTERNAL OPERATORS ########################


std::ostream& operator<<(std::ostream& os, const SubTriangleFace& face){
    os<<"("<<face[0]<<", "<<face[1]<<", "<<face[2]<<")";
    return os;
}

}



