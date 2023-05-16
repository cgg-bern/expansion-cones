#include "ExpansionCone.hh"
#include "ProgEmbeddingHelpers.hh"
#include <queue>

#include <OpenVolumeMesh/FileManager/FileManager.hh>


namespace OpenVolumeMesh{


    ExpansionCone::ExpansionCone()
            : tip_vertices_prop_(request_vertex_property<bool>()),
              cone_to_mesh_v_handle_prop_(request_vertex_property<VertexHandle>()),
              cone_to_mesh_c_handle_prop_(request_cell_property<CellHandle>()),
              vertex_position_prop_(request_vertex_property<VertexPosition>()){

    }

    ExpansionCone::ExpansionCone(const ExpansionCone& other_cone) : ExpansionCone(){
        (*this) = other_cone;
    }


    ExpansionCone::ExpansionCone(const TetrahedralMesh& mesh) : ExpansionCone(){

        for(int i(0); i<(int)mesh.n_vertices(); i++){
            VertexHandle v(i);
            if(!mesh.is_deleted(v)){
                add_vertex(v, vec2vec(mesh.vertex(v)));
            }
        }

        for(int i(0); i<(int)mesh.n_cells(); i++){
            CellHandle c(i);
            if(!mesh.is_deleted(c)){
                add_cell(c, mesh.get_cell_vertices(c));
            }
        }

        if(n_vertices() != mesh.n_logical_vertices()){
            std::cout<<" ERROR - couldn't convert tet mesh into ExpansionCone: wrong number of vertices"<<std::endl;
            this->clear();
            return;
        }

        if(n_edges() != mesh.n_logical_edges()){
            std::cout<<" ERROR - couldn't convert tet mesh into ExpansionCone: wrong number of edges"<<std::endl;
            this->clear();
            return;
        }

        if(n_faces() != mesh.n_logical_faces()){
            std::cout<<" ERROR - couldn't convert tet mesh into ExpansionCone: wrong number of faces"<<std::endl;
            this->clear();
            return;
        }

        if(n_cells() != mesh.n_logical_cells()){
            std::cout<<" ERROR - couldn't convert tet mesh into ExpansionCone: wrong number of cells"<<std::endl;
            this->clear();
            return;
        }
    }



    void ExpansionCone::operator=(const ExpansionCone& other_cone){

        /*std::cout<<" called copy operator"<<std::endl;
        std::cout<<" initial cone: "<<std::endl;
        std::cout<<"  - vertices: "<<n_vertices()<<std::endl;
        std::cout<<"  - edges: "<<n_edges()<<std::endl;
        std::cout<<"  - faces: "<<n_faces()<<std::endl;
        std::cout<<"  - cells: "<<n_cells()<<std::endl;


        std::cout<<" other cone: "<<std::endl;
        std::cout<<"  - vertices: "<<other_cone.n_vertices()<<std::endl;
        std::cout<<"  - edges: "<<other_cone.n_edges()<<std::endl;
        std::cout<<"  - faces: "<<other_cone.n_faces()<<std::endl;
        std::cout<<"  - cells: "<<other_cone.n_cells()<<std::endl;*/


        //copy the mesh itself
        TetrahedralMesh::operator=(other_cone);


        /*std::cout<<" cone after copy: "<<std::endl;
        std::cout<<"  - vertices: "<<n_vertices()<<std::endl;
        std::cout<<"  - edges: "<<n_edges()<<std::endl;
        std::cout<<"  - faces: "<<n_faces()<<std::endl;
        std::cout<<"  - cells: "<<n_cells()<<std::endl;*/

        cone_tip_vertices_ = other_cone.cone_tip_vertices_;

        //and then copy the properties
        for(auto v: other_cone.vertices()){
            /*if(other_cone.tip_vertices_prop_[v]){
                cone_tip_vertices_.insert(v);
            }*/
            tip_vertices_prop_[v]          = other_cone.tip_vertices_prop_[v];

            cone_to_mesh_v_handle_prop_[v] = other_cone.cone_to_mesh_v_handle_prop_[v];
            vertex_position_prop_[v]       = other_cone.vertex_position_prop_[v];
        }

        for(auto c: other_cone.cells()){
            cone_to_mesh_c_handle_prop_[c] = other_cone.cone_to_mesh_c_handle_prop_[c];
        }

        mesh_to_cone_v_handle_map_ = other_cone.mesh_to_cone_v_handle_map_;
        //std::cout<<" called EC copy ctor"<<std::endl;
    }


    EXPANSION_CHECK_RESULT ExpansionCone::is_expandable(VertexPosition& pos,
                                                        bool print_debug){

        auto topo_exp = is_topo_expandable(print_debug);

        if(topo_exp == IS_EXPANDABLE){
            auto geo_exp = is_geo_expandable(pos, print_debug);
            return geo_exp;
        }else if(topo_exp == EXPANDABILITY_CHECK_ERROR){
            return EXPANDABILITY_CHECK_ERROR;
        }else{
            return topo_exp;
        }
    }



    void ExpansionCone::remove_tips(){
        for(auto tip_v: cone_tip_vertices()){
            delete_vertex(tip_v);
        }
    }



    VertexHandle ExpansionCone::merge_tips(){

        std::vector<std::vector<VertexHandle>> new_cells_vertices;
        //add the new tip at the same location as the tips'
        auto new_tip = add_vertex(VertexHandle(-1),
                                  this->vertex(*cone_tip_vertices().begin()));
        std::cout<<" - added new tip at "<<vec2vec(vertex(new_tip))<<std::endl;

        for(auto hf: halffaces()){
            auto op_vertex = halfface_opposite_vertex(hf);
            if(op_vertex.idx() != -1 && tip_vertices_prop_[op_vertex]){
                auto hf_vertices = get_halfface_vertices(hf);
                if(!tip_vertices_prop_[hf_vertices[0]] &&
                        !tip_vertices_prop_[hf_vertices[1]] &&
                        !tip_vertices_prop_[hf_vertices[2]]){
                    hf_vertices.push_back(new_tip);
                    new_cells_vertices.push_back(hf_vertices);
                    //std::cout<<" -- added cell "<<hf_vertices<<" to list"<<std::endl;
                }
            }
        }

        //remove everything
        remove_tips();
        for(auto e: edges()){
            //std::cout<<" -- deleted edge "<<edge(e)<<std::endl;
            delete_edge(e);
        }

        //std::cout<<" - cone after edge deletion : "; print_details();

        //and add the new cells
        for(const auto& new_cell_vertices: new_cells_vertices){
            auto ch = TetrahedralMesh::add_cell(new_cell_vertices);
            //not sure it's really useful but just in case
            cone_to_mesh_c_handle_prop_[ch] = CellHandle(-1);
        }

        //clear the current tips
        cone_tip_vertices_.clear();

        //and replace them with the new one
        set_as_tip_internal(new_tip);

        return new_tip;
    }


    VertexHandle ExpansionCone::find_removable_vertex_with_cell_valence_lower_than_3(const std::vector<VertexHandle>& candidate_vertices,
                                                                                     const VertexPropertyT<bool>& to_ignore_prop,
                                                                                     const int min_index){


        const bool print_debug(false);
        if(print_debug){
            std::cout<<" --------------------------------------------------"<<std::endl;
            std::cout<<" ----- looking for best candidate to remove..."<<std::endl;
            std::cout<<" ----- min index = "<<min_index<<std::endl;
        }

        /*ExpansionCone uncollapsed_cone = *this;
        uncollapsed_cone.enable_fast_deletion(false);
        for(auto candidate_v: candidate_vertices){
            uncollapsed_cone.delete_vertex(collapsed_v);
            std::cout<<" -- deleted collapsed vertex "<<collapsed_v<<std::endl;
        }*/

        //std::cout<<" - current cone: "; print_details();
        //std::cout<<" - current cone: "<<(*this)<<std::endl;

        //std::cout<<" - base boundary egdes: "<<std::endl;
        //first, define which are on the boundary (because only those can be potentially removed)
        auto base_boundary_prop = request_vertex_property<bool>();
        for(auto e: edges()){
            if(is_base_boundary_edge(e)){
                base_boundary_prop[edge(e).from_vertex()] = true;
                base_boundary_prop[edge(e).to_vertex()]   = true;
                //std::cout<<" -- "<<edge(e)<<std::endl;
            }
        }

        if(print_debug){
            std::cout<<" -------------"<<std::endl;
        }

        VertexHandle best_removable_vertex(-1);
        int best_removable_vertex_valence(std::numeric_limits<int>::max());

        for(auto candidate: candidate_vertices){

            if(print_debug){
                std::cout<<" ----------------------------"<<std::endl;
                std::cout<<" -- checking candidate "<<candidate<<std::endl;
                std::cout<<" -- with neighbors: "<<std::endl;
                for(auto vv_it = vv_iter(candidate); vv_it.valid(); vv_it++){
                    if(!is_cone_tip(*vv_it)){
                        std::cout<<"   --- "<<(*vv_it)<<" with neighbors:"<<std::endl;
                        for(auto vv_it2 = vv_iter(*vv_it); vv_it2.valid(); vv_it2++){
                            std::cout<<"           --- "<<(*vv_it2)<<std::endl;

                        }
                    }
                }
            }

            if(candidate.idx() < min_index){
                if(print_debug){
                    std::cout<<" ---> index greater than min index "<<min_index<<std::endl;
                }
                continue;
            }

            if(is_deleted(candidate)){

                if(print_debug){
                    std::cout<<" ---> already deleted"<<std::endl;
                }
                continue;
            }

            if(to_ignore_prop[candidate]){
                if(print_debug){
                    std::cout<<" ---> already handled"<<std::endl;
                }
                continue;
            }

            if(!base_boundary_prop[candidate]){
                if(print_debug){
                    std::cout<<" - candidate "<<candidate<<" is not on the base's boundary, skipping"<<std::endl;
                }
                continue;
            }

            int candidate_valence(cell_valence(candidate));
            if(!candidate_valence){
                std::cout<<" ERROR - vertex "<<candidate<<" has cell valence "<<candidate_valence<<std::endl;
                std::cout<<" cone details: "; print_details();
                return VertexHandle(-1);
            }

            //if the best removable candidate already has a lower cell valence that
            //this new candidate there's no need to see if this one's better
            if(best_removable_vertex.idx() != -1 &&
                    candidate_valence >= best_removable_vertex_valence){

                if(print_debug){
                    std::cout<<" -- save valence as best one, skipping"<<std::endl;
                }
                continue;
            }

            if(print_debug){
                std::cout<<" ---> valence is "<<candidate_valence<<", best is "<<best_removable_vertex_valence<<std::endl;
            }

            ExpansionCone cone_copy = *this;
            cone_copy.enable_fast_deletion(false);
            cone_copy.delete_vertex(candidate);
            //cone_copy.collect_garbage();
            //std::cout<<" -- cone copy after deleting vertex "<<candidate<<": "; cone_copy.print_details();
            //std::cout<<" -- cone copy: "<<cone_copy<<std::endl;

            auto topo_result = cone_copy.is_topo_expandable();
            if(topo_result == -1){
                std::cout<<" error while checking expandability of cone copy with vertex "<<candidate<<" removed"<<std::endl;
                return VertexHandle(-1);
            }else if(!topo_result){
                //if it's still topo-expandable we can update the best candidate
                //since we already checked that its cell valence was lower
                best_removable_vertex = candidate;
                best_removable_vertex_valence = candidate_valence;

                if(print_debug){
                    std::cout<<" ---> updated best vertex to "<<candidate<<" with valence "<<candidate_valence<<std::endl;
                }
                if(candidate_valence <= 2){
                    break;
                }
            }else{

                if(print_debug){
                    std::cout<<" ---> not topo-expandable, result: "<<topo_result<<std::endl;
                }
            }
        }

        if(print_debug){
            std::cout<<" -----> best vertex: "<<best_removable_vertex<<" with valence "<<best_removable_vertex_valence<<std::endl;
            std::cout<<" --------------------------------------------------"<<std::endl;
        }

        if(best_removable_vertex.idx() == -1){
            std::cout<<" --> ERROR - couldn't find a boundary vertex to remove! "<<std::endl;
        }


        return best_removable_vertex;

    }


    int ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(TetrahedralMesh& mesh,
                                                                    const VertexPropertyT<VertexPosition>& exact_vertex_position_prop,
                                                                    const VertexHandle& vh,
                                                                    ExpansionCone& one_ring_EC){
        return set_up_1_ring_neighborhood_as_expansion_cone(mesh,
                                                            exact_vertex_position_prop,
                                                            std::vector<VertexHandle>({vh}),
                                                            one_ring_EC);
    }


    int ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(TetrahedralMesh& mesh,
                                                                    const VertexPropertyT<VertexPosition>& exact_vertex_position_prop,
                                                                    const std::vector<VertexHandle>& vs,
                                                                    ExpansionCone& one_ring_EC){

        one_ring_EC.clear();

        auto added_vertex_prop = mesh.request_vertex_property<bool>();
        auto added_cell_prop = mesh.request_cell_property<bool>();

        for(auto v: vs){
            one_ring_EC.add_vertex(v, exact_vertex_position_prop[v], true);
            added_vertex_prop[v] = true;
        }

        for(auto v: vs){
            //std::cout<<" --------------------------- checking cells around vertex "<<v<<std::endl;
            for(auto vc_it = mesh.vc_iter(v); vc_it.valid(); vc_it++){
                //std::cout<<"  -- checking cell "<<mesh.get_cell_vertices(*vc_it)<<std::endl;
                if(added_cell_prop[*vc_it]){
                    //std::cout<<"   ---> already added "<<std::endl;
                    continue;
                }
                //std::cout<<"   --- adding its vertices"<<std::endl;
                auto c_vertices = mesh.get_cell_vertices(*vc_it);
                for(auto cv: c_vertices){
                    //std::cout<<"      ---- adding vertex "<<cv<<std::endl;
                    if(!added_vertex_prop[cv]){
                        one_ring_EC.add_vertex(cv, exact_vertex_position_prop[cv]);
                        added_vertex_prop[cv] = true;
                        //std::cout<<"        --> added vertex "<<cv<<std::endl;
                    }else{
                        //std::cout<<"        --> already added "<<std::endl;
                    }
                }

                auto cone_c = one_ring_EC.add_cell(*vc_it, c_vertices);
                if(cone_c.idx() == -1){
                    std::cout<<" error while adding cell "<<c_vertices<<" incident to vertex "<<v<<std::endl;
                    return -1;
                }
                added_cell_prop[*vc_it] = true;
                //std::cout<<"  --> updated cell count = "<<one_ring_EC.n_cells()<<std::endl;
            }
        }

        return 0;
    }



    int ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(TetrahedralMesh& mesh,
                                                                    const std::vector<VertexHandle>& vs,
                                                                    ExpansionCone& one_ring_EC){

        one_ring_EC.clear();

        auto added_vertex_prop = mesh.request_vertex_property<bool>();
        auto added_cell_prop = mesh.request_cell_property<bool>();

        for(auto v: vs){
            one_ring_EC.add_vertex(v, vec2vec(mesh.vertex(v)), true);
            added_vertex_prop[v] = true;
        }


        for(auto v: vs){
            for(auto vc_it = mesh.vc_iter(v); vc_it.valid(); vc_it++){
                if(added_cell_prop[*vc_it]){
                    continue;
                }
                auto c_vertices = mesh.get_cell_vertices(*vc_it);
                for(auto cv: c_vertices){
                    if(!added_vertex_prop[cv]){
                        one_ring_EC.add_vertex(cv, vec2vec(mesh.vertex(cv)));
                        added_vertex_prop[cv] = true;
                    }
                }

                auto cone_c = one_ring_EC.add_cell(*vc_it, c_vertices);
                if(cone_c.idx() == -1){
                    std::cout<<" error while adding cell "<<c_vertices<<" incident to vertex "<<v<<std::endl;
                    return -1;
                }
                added_cell_prop[*vc_it] = true;
            }
        }

        return 0;
    }


    EXPANSION_CHECK_RESULT ExpansionCone::is_topo_expandable(bool print_debug) {

        //print_debug = true;

        /*auto cone_center_vertex = mesh_to_cone_handle(center_vertex);
        if(cone_center_vertex.idx() == -1){
            std::cerr<<" ERROR - couldn't find mesh center vertex "<<center_vertex<<" in cone"<<std::endl;
            return EXPANDABILITY_CHECK_ERROR;
        }*/

        if(print_debug){
            std::cout<<" ------ checking for topo-expandability of cone "<<std::endl;
            print_details();
        }

        //if the cone contains no cell it's not expandable
        if(!n_logical_cells()){
            return NO_CELL_IN_EC;
        }


        //extracting the base from the cone
        TetrahedralMesh cone_base = static_cast<TetrahedralMesh>(*this);


        for(auto cone_tip: cone_tip_vertices()) {
            //std::cout<<" deleted vertex "<<cone_tip<<std::endl;
            cone_base.delete_vertex(cone_tip);
        }


        if(print_debug) {
            std::cout<<" --------------------- "<<std::endl;
            std::cout << " vertices: ("<<n_vertices()<<")"<< std::endl;
            for (auto v: vertices()) {
                std::cout << " - " << v << " -> " << cone_to_mesh_v_handle_prop_[v] << std::endl;
            }
            std::cout<<" --------------------- "<<std::endl;

            std::cout<<" cone base: ("<<cone_base.n_logical_vertices()<<
                       ","<<cone_base.n_logical_edges()<<
                       ","<<cone_base.n_logical_faces()<<
                       ","<<cone_base.n_logical_cells()<<")"<<std::endl;
            for(auto f: cone_base.faces()){
                std::cout<<" - face "<<f<<" : "<<cone_base.get_halfface_vertices(cone_base.halfface_handle(f, 0))<<std::endl;
            }
            for(auto e: cone_base.edges()){
                std::cout<<" - edge "<<e<<
                           " : "<<cone_base.from_vertex_handle(cone_base.halfedge_handle(e, 0))<<
                           "-"<<cone_base.to_vertex_handle(cone_base.halfedge_handle(e, 0))<<std::endl;
            }
            //print_debug = false;
        }

        //cone_base.collect_garbage();

        if(cone_base.n_logical_cells()){


            std::cout<<" ERROR - cone base ("<<cone_base.n_logical_vertices()<<
                     ","<<cone_base.n_logical_edges()<<
                     ","<<cone_base.n_logical_faces()<<
                     ","<<cone_base.n_logical_cells()<<
                     ") of cone "<<(*this)<<" contains cells"<<std::endl;
            std::cout<<" cone base vertices: "<<std::endl;
            for(auto v: cone_base.vertices()){
                std::cout<<v<<"("<<cone_base.is_deleted(v)<<") -> "<<cone_to_mesh_handle(v)<<std::endl;
            }
            std::cout<<std::endl;

            //cone_base.collect_garbage();

            std::cout<<" cone base cells: "<<std::endl;
            for(int i(0); i<(int)cone_base.n_cells(); i++){
                CellHandle c(i);
                if(!cone_base.is_deleted(c)) {

                    std::cout<<" checking cone base cell "<<c<<std::endl;

                    auto cell_vertices = cone_base.get_cell_vertices(c);
                    if(cell_vertices.empty()){
                        std::cout<<" ERROR - empty cell "<<c<<std::endl;
                        return EXPANDABILITY_CHECK_ERROR;
                    }
                    std::cout << " - " << c << "(" << cone_base.is_deleted(c) << "): " << cell_vertices << " -> ";
                    for(auto cone_v: cell_vertices){
                        std::cout<<cone_to_mesh_handle(cone_v)<<" ";
                    }
                    std::cout << std::endl;
                }
            }

            std::cout<<" full cone cells: "<<std::endl;
            for(auto c: this->cells()){

                auto cell_vertices = this->get_cell_vertices(c);
                std::cout << " - " << c << "(" << this->is_deleted(c) << "): " << cell_vertices << " -> ";
                for(auto cone_v: cell_vertices){
                    std::cout<<cone_to_mesh_handle(cone_v)<<" ";
                }
                std::cout << std::endl;

            }


            return EXPANDABILITY_CHECK_ERROR;
        }


        if(print_debug) {
            std::cout<<" ----------- cone:"<<std::endl;
            std::cout << " vertices: " << n_logical_vertices() << std::endl;
            std::cout << "    edges: " << n_logical_edges() << std::endl;
            std::cout << "    faces: " << n_logical_faces() << std::endl;
            std::cout << "    cells: " << n_logical_cells() << std::endl;
            std::cout<<" ----------- base:"<<std::endl;
            std::cout << " vertices: " << cone_base.n_logical_vertices() << std::endl;
            std::cout << "    edges: " << cone_base.n_logical_edges() << std::endl;
            std::cout << "    faces: " << cone_base.n_logical_faces() << std::endl;
        }



        //std::cout<<" euler characteristic = "<<euler_characteristic<<std::endl;
        //if the cone base is not disc topology then it's not expandable
        if(euler_characteristic() != 1){
            //std::cout<<" --> euler characteristic = "<<cone_base.euler_characteristic()<<
            //           " != 1 -> cone base is not disc topology and thus cone is not expandable"<<std::endl;
            return EULER_CHARAC_IS_NOT_ONE;
        }

        int base_euler_characteristic = (int)cone_base.n_logical_vertices()
                                      - (int)cone_base.n_logical_edges()
                                      + (int)cone_base.n_logical_faces();

        if(base_euler_characteristic != 1){
            //std::cout<<" --> euler characteristic = "<<cone_base.euler_characteristic()<<
            //           " != 1 -> cone base is not disc topology and thus cone is not expandable"<<std::endl;
            return CONE_BASE_EULER_CHARAC_IS_NOT_ONE;
        }

        //std::cout<<" --> euler characteristic = 1"<<std::endl;


        auto visited_prop = cone_base.request_vertex_property<bool>();
        std::queue<VertexHandle> to_visit;
        //use the first non-tip vertex as starting point
        for(auto base_v: cone_base.vertices()){
            if(!tip_vertices_prop_[base_v]){
                to_visit.push(base_v);
                break;
            }
        }
        int visited_count(0);

        while(!to_visit.empty()){
            auto current_v = to_visit.front();
            to_visit.pop();

            if(!visited_prop[current_v]) {
                visited_prop[current_v] = true;
                visited_count++;

                //std::cout << " - current v: " << current_v << std::endl;

                //quick check for spikes
                if (cone_base.valence(current_v) == 1) {
                    //std::cout << " - vertex " << current_v << "->" << cone_to_mesh_handle(current_v) <<
                    //             " has valence 1 and the cone thus contains an antenna" << std::endl;
                    return ANTENNA;

                } else if (!cone_base.valence(current_v)) {
                    //std::cout << " - vertex " << current_v << "->" << cone_base.cone_to_mesh_handle(current_v) <<
                    //          " has valence 0 and the cone thus contains a spike" << std::endl;
                    return MULTIPLE_CONNECTED_COMPONENTS;
                }

                for (auto out_he: cone_base.outgoing_halfedges(current_v)) {
                    auto neighbor = cone_base.to_vertex_handle(out_he);

                    if (!visited_prop[neighbor]) {
                        to_visit.push(neighbor);
                        //std::cout << " -- neighbor " << neighbor << " is not visited yet, adding it to the list"<< std::endl;
                    }
                }
            }
        }

        if(visited_count != (int)cone_base.n_logical_vertices()){
            //std::cout<<" --> cone base contains "<<cone_base.n_logical_vertices()<<
            //           " vertices but only "<<visited_count<<" were visited"<<std::endl;
            return MULTIPLE_CONNECTED_COMPONENTS;
        }


        //then check for antennas (zero-valence edges)
        for(auto e: cone_base.edges()){
            if(!cone_base.valence(e)){
                //std::cout<<" --> edge "<<edge(e)<<" has 0-valence"<<std::endl;
                return ANTENNA;
            }
        }





        //std::cout<<" --> all edges visited"<<std::endl;


        // and finally that all base vertices have no more than 1 more incident edge
        // that they have incident faces.
        // This checks that all tets surrounding the spokes are contiguously expanded
        // i.e. that the base doesn't look like
        // ___
        // \./
        // /_|

        for(auto v: cone_base.vertices()) {

            if(print_debug){
                std::cout<<" - checking vertex "<<v<<std::endl;
            }
            auto edge_valence = cone_base.valence(v);
            int face_valence(0);
            for (auto vf_it = cone_base.vf_iter(v); vf_it.valid(); vf_it++) {
                face_valence++;
            }

            if(print_debug){
                std::cout<<" edge valence = "<<edge_valence<<", face valence = "<<face_valence<<std::endl;
            }

            if ((edge_valence - face_valence) > 1) {
                /*std::cout<<" --> spoke edge "<<halfedge(out_he_it)<<
                               " has face valence "<<face_valence<<
                               " and cell valence "<<cell_valence<<std::endl;
                    */
                return SPOKES_NOT_CONTINUOUSLY_EXPANDED;
            }
        }


        //finally, a cluster-specific test is to check that the set of tip vertices are also disc-topology
        if(cone_tip_vertices_.size() > 1) {

            if(print_debug){
                std::cout<<" RUNNING CLUSTER CHECKS"<<std::endl;
            }


            TetrahedralMesh cone_tips = static_cast<TetrahedralMesh>(*this);

            for(auto v: vertices()){
                if(!tip_vertices_prop_[v]){
                    cone_tips.delete_vertex(v);
                 }
            }

            if(print_debug){
                std::cout<<" - cone tips mesh: "<<std::endl;
                std::cout<<" VERTICES:"<<std::endl;
                for(auto v: cone_tips.vertices()){
                    std::cout<<" -- "<<v<<"->"<<cone_to_mesh_handle(v)<<std::endl;
                }
                std::cout<<" FACES:"<<std::endl;
                for(auto f: cone_tips.faces()){
                    std::cout<<" -- "<<f<<": "<<cone_tips.get_halfface_vertices(cone_tips.halfface_handle(f, 0))<<std::endl;
                }
            }

            if(cone_tips.n_logical_vertices() != cone_tip_vertices().size()){
                std::cout<<" ERROR - cone tips mesh contains "<<cone_tips.n_logical_vertices()<<
                           " vertices but there are "<<cone_tip_vertices().size()<<" cone tip vertices"<<std::endl;
                std::cout<<" cone tips : "<<cone_tip_vertices()<<std::endl;
                return EXPANDABILITY_CHECK_ERROR;
            }


            int tip_e_count(cone_tips.n_logical_edges());
            int tip_f_count(cone_tips.n_logical_faces());
            int tip_c_count(cone_tips.n_logical_cells());


            if(((int)cone_tip_vertices_.size() - tip_e_count + tip_f_count - tip_c_count) != 1){
                if(print_debug) {
                    std::cout << " --> Euler characteristic = " << cone_tip_vertices_.size()<<
                                 " - "<<tip_e_count<<" + "<<tip_f_count<<" - "<<tip_c_count<<
                                 " = "<<((int)cone_tip_vertices_.size() - tip_e_count + tip_f_count - tip_c_count) << std::endl;
                }
                return CLUSTER_TIPS_EULER_CHARACTERISTIC_IS_NOT_ONE;
            }


            auto visited_tip_prop = cone_tips.request_vertex_property<bool>();
            std::queue<VertexHandle> tips_to_visit;
            tips_to_visit.push(*cone_tip_vertices_.begin());
            while(!tips_to_visit.empty()){
                auto current_tip = tips_to_visit.front();
                tips_to_visit.pop();

                if(!visited_tip_prop[current_tip]) {
                    for (auto out_he: cone_tips.outgoing_halfedges(current_tip)){
                        auto neighbor = cone_tips.to_vertex_handle(out_he);
                        if(!visited_tip_prop[neighbor]){
                            tips_to_visit.push(neighbor);
                        }
                    }
                }

                visited_tip_prop[current_tip] = true;
            }

            for(auto tip: cone_tip_vertices_){
                if(!visited_tip_prop[tip]){
                    return CLUSTER_TIPS_ARE_NOT_A_SINGLE_CONNECTED_COMPONENT;
                }
            }

            for(auto f: faces()){
                if(is_boundary(halfface_handle(f,0)) && is_boundary(halfface_handle(f,1))){
                    return CLUSTER_NOT_ALL_FACES_ARE_INCIDENT_TO_AT_LEAST_ONE_TET;
                }
            }
        }



        //std::cout<<" ... done -> topo-expandable"<<std::endl;
        return IS_EXPANDABLE;
    }



    EXPANSION_CHECK_RESULT ExpansionCone::is_geo_expandable(VertexPosition& new_position,
                                                            bool print_debug) const {



        int result = find_chebyshev_center(new_position, print_debug);

        if(!result){
            return IS_EXPANDABLE;
        }else if(result == -1){
            return EXPANDABILITY_CHECK_ERROR;
        }else{
            return IS_NOT_GEO_EXPANDABLE;
        }
    }



    bool ExpansionCone::same_entities_count(const ExpansionCone& other) const{

        return n_vertices() == other.n_vertices()  &&
               n_edges() == other.n_edges() &&
               n_faces() == other.n_faces() &&
               n_cells() == other.n_cells();
    }








    VertexPosition ExpansionCone::vertex(const VertexHandle& v) const{

        return vertex_position_prop_[v];
    }


    void ExpansionCone::set_vertex(const VertexHandle& v,
                                   const VertexPosition& pos){
        vertex_position_prop_[v] = pos;
        TetrahedralMesh::set_vertex(v, vec2vec(pos));
    }

    const std::set<VertexHandle>& ExpansionCone::cone_tip_vertices() const{
        return cone_tip_vertices_;
    }


    VertexHandle ExpansionCone::add_vertex(const VertexHandle& mesh_vertex,
                                           const VertexPosition& pos,
                                           bool tip_vertex){

        auto cone_vertex = TetrahedralMesh::add_vertex(vec2vec(pos));
        set_vertex(cone_vertex, pos);
        cone_to_mesh_v_handle_prop_[cone_vertex] = mesh_vertex;
        mesh_to_cone_v_handle_map_.insert({mesh_vertex, cone_vertex});

        if (n_vertices() == 1 || tip_vertex) {
            tip_vertices_prop_[cone_vertex] = true;
            cone_tip_vertices_.insert(cone_vertex);
        }

        return cone_vertex;
    }


    void ExpansionCone::set_as_tip(const VertexHandle& mesh_vertex){
        auto cone_vertex = mesh_to_cone_handle(mesh_vertex);
        if(cone_vertex.is_valid()){
            set_as_tip_internal(cone_vertex);
        }
    }


    void ExpansionCone::set_as_tip_internal(const VertexHandle& cone_vertex){
        tip_vertices_prop_[cone_vertex] = true;
        cone_tip_vertices_.insert(cone_vertex);
    }


    CellHandle ExpansionCone::add_cell(const CellHandle& mesh_cell,
                                       const std::vector<VertexHandle>& mesh_vertices){

        if(mesh_vertices.size() != 4){
            std::cerr<<" ERROR - added cell with "<<mesh_vertices.size()<<" vertices"<<std::endl;
            return CellHandle(-1);
        }

        /*std::cout<<" cone vertices: "<<std::endl;
        for(auto v: vertices()){
            std::cout<<" -- "<<v<<" -> "<<cone_to_mesh_v_handle_prop_[v]<<std::endl;
        }*/

        //std::cout<<" mesh vertices: "<<std::endl;
        std::vector<VertexHandle> cone_vertices;
        for(auto v: mesh_vertices){
            //std::cout<<" -- "<<v<<"/"<<mesh_vertices[i]<<
            //           " -> "<<mesh_to_cone_v_handle_map_.find(v)->second<<std::endl;

            auto it = mesh_to_cone_v_handle_map_.find(v);
            if(it == mesh_to_cone_v_handle_map_.end()){
                std::cerr<<" ERROR - couldn't find mesh vertex "<<v<<" in mesh to cone map"<<std::endl;
                return CellHandle(-1);
            }

            cone_vertices.push_back(it->second);
        }

        //std::cout<<" adding cell with vertices "<<mesh_vertices<<std::endl;
        //std::cout<<" --> corresponding to cone vertices "<<cone_vertices<<std::endl;

        auto cone_cell = TetrahedralMesh::add_cell(cone_vertices);

        cone_to_mesh_c_handle_prop_[cone_cell] = mesh_cell;

        return cone_cell;
    }


    FaceHandle ExpansionCone::add_face(const std::vector<VertexHandle>& mesh_vertices){
        if(mesh_vertices.size() != 3){
            std::cerr<<" ERROR - added face with "<<mesh_vertices.size()<<" vertices"<<std::endl;
            return FaceHandle(-1);
        }

        /*std::cout<<" cone vertices: "<<std::endl;
        for(auto v: vertices()){
            std::cout<<" -- "<<v<<" -> "<<cone_to_mesh_v_handle_prop_[v]<<std::endl;
        }*/

        //std::cout<<" mesh vertices: "<<std::endl;
        std::vector<VertexHandle> cone_vertices;
        for(auto v: mesh_vertices){
            //std::cout<<" -- "<<v<<"/"<<mesh_vertices[i]<<
            //           " -> "<<mesh_to_cone_v_handle_map_.find(v)->second<<std::endl;

            auto it = mesh_to_cone_v_handle_map_.find(v);
            if(it == mesh_to_cone_v_handle_map_.end()){
                std::cerr<<" ERROR - couldn't find mesh vertex "<<v<<" in mesh to cone map"<<std::endl;
                return FaceHandle(-1);
            }

            cone_vertices.push_back(it->second);
        }

        auto hf = TetrahedralMesh::halfface(cone_vertices);
        if(hf.idx() != -1){
            return face_handle(hf);
        }
        //std::cout<<" adding face with vertices "<<mesh_vertices<<std::endl;
        //std::cout<<" --> corresponding to cone vertices "<<cone_vertices<<std::endl;

        return TetrahedralMesh::add_face(cone_vertices);
    }


    EdgeHandle ExpansionCone::add_edge(const VertexHandle& mesh_from_vertex,
                                       const VertexHandle& mesh_to_vertex){

        auto it = mesh_to_cone_v_handle_map_.find(mesh_from_vertex);
        if(it == mesh_to_cone_v_handle_map_.end()){
            std::cerr<<" ERROR - couldn't find mesh from-vertex "<<mesh_from_vertex<<" in mesh to cone map"<<std::endl;
            return EdgeHandle(-1);
        }
        auto cone_from_vertex = it->second;

        it = mesh_to_cone_v_handle_map_.find(mesh_to_vertex);
        if(it == mesh_to_cone_v_handle_map_.end()){
            std::cerr<<" ERROR - couldn't find mesh to-vertex "<<mesh_to_vertex<<" in mesh to cone map"<<std::endl;
            return EdgeHandle(-1);
        }
        auto cone_to_vertex = it->second;

        return TetrahedralMesh::add_edge(cone_from_vertex,
                                         cone_to_vertex);
    }

    VertexHandle ExpansionCone::split_edge(const HalfEdgeHandle& cone_e){
        return TetrahedralMesh::split_edge(cone_e);
    }

    VertexHandle ExpansionCone::split_edge(const EdgeHandle& cone_e){
        return TetrahedralMesh::split_edge(cone_e);
    }

    VertexHandle ExpansionCone::split_edge(const EdgeHandle& cone_e,
                                           const VertexHandle& mesh_mid_v){

        auto cone_mid_v = split_edge(cone_e);
        auto mid_v_pos = (vertex(edge(cone_e).from_vertex()) + vertex(edge(cone_e).to_vertex()))/2;
        this->set_vertex(cone_mid_v, mid_v_pos);

        cone_to_mesh_v_handle_prop_[cone_mid_v] = mesh_mid_v;
        mesh_to_cone_v_handle_map_.insert({mesh_mid_v, cone_mid_v});

        return cone_mid_v;
    }

#if 0
    void ExpansionCone::merge(const ExpansionCone& other){
        if(!other.n_vertices()){
            return;
        }


    }

#endif


    VertexHandle ExpansionCone::cone_to_mesh_handle(const VertexHandle& cone_handle) const{
        if(cone_handle.idx() >= (int)n_vertices()){
            return VertexHandle(-1);
        }
        return cone_to_mesh_v_handle_prop_[cone_handle];
    }

    VertexHandle ExpansionCone::mesh_to_cone_handle(const VertexHandle& mesh_handle) const{
        auto it = mesh_to_cone_v_handle_map_.find(mesh_handle);
        if(it == mesh_to_cone_v_handle_map_.end()){
            //std::cerr<<" ERROR - couldn't find mesh vertex "<<mesh_handle<<" in mesh to cone map"<<std::endl;
            return VertexHandle(-1);
        }
        return it->second;
    }


    CellHandle ExpansionCone::cone_to_mesh_handle(const CellHandle& cone_handle) const{
        if(cone_handle.idx() >= (int)n_cells()){
            return CellHandle(-1);
        }
        return cone_to_mesh_c_handle_prop_[cone_handle];
    }


    void ExpansionCone::clear(){
        TetrahedralMesh::clear_mesh_props();
        (*this) = ExpansionCone();
    }


    bool ExpansionCone::contains_flipped_tets() const{
        return ExactBadTetFinder::meshContainsFlippedTets(*this, vertex_position_prop_);
    }

    bool ExpansionCone::contains_degenerate_tets() const{
        return ExactBadTetFinder::meshContainsDegenerateTets(*this, vertex_position_prop_);
    }

    bool ExpansionCone::is_base_boundary_vertex(const VertexHandle& vh) const{
        for(auto out_he: outgoing_halfedges(vh)){
            if(is_base_boundary_edge(edge_handle(out_he))){
                return true;
            }
        }
        return false;
    }


    bool ExpansionCone::is_base_boundary_edge(const EdgeHandle& eh) const{
        return valence(eh) == 2 && is_base_edge(eh);
    }

    bool ExpansionCone::is_base_edge(const EdgeHandle& eh) const{
        return !tip_vertices_prop_[edge(eh).from_vertex()] &&
                !tip_vertices_prop_[edge(eh).to_vertex()];
    }


   /* int ExpansionCone::find_max_min_volume_center(VertexPosition& new_position,
                                                  int library_to_use,
                                                  double eps,
                                                  bool print_debug){

        return find_max_min_volume_center_with_CGAL(new_position, eps, print_debug);

    }*/



    int ExpansionCone::find_max_min_volume_center_with_CGAL(VertexPosition& new_position,
                                                            double eps,
                                                            bool print_debug) const{

        if(print_debug){
            std::cout<<" printing CGAL debug information"<<std::endl;
        }

        new_position = {0,0,0};

        using Program  = CGAL::Quadratic_program<double>;
        using Solution = CGAL::Quadratic_program_solution<ExactType>;

        if(false && print_debug){
            std::cout<<" gathering cone "<<*this<<" boundary normals..."<<std::endl;
            print_details();

            std::cout<<" vertices: "<<std::endl;
            for(auto v: vertices()){
                std::cout<<" - "<<v<<" -> "<<cone_to_mesh_handle(v)<<" : "<<std::setprecision(20)<<vec2vec(vertex(v))<<std::endl;
            }

            std::cout<<" cells: "<<std::endl;
            for(auto c: cells()){
                std::cout<<" - "<<c<<": "<<get_cell_vertices(c)<<std::endl;
            }
        }

        int max_precision(0);


        //create solver
        // by default, we have an LP with Ax <= b
        Program lp(CGAL::SMALLER, false, 0, false, 0);

        double lower_bound(0);

        //slack variable weight
        double w(1);

        // now set the non-default entries
        const int X(0), Y(1), Z(2);//, R(3);
        int i(0);
        for(auto hf: halffaces()){

            if(is_boundary(hf)){
                VertexPosition normal, centroid;

                if(triangle_normal_and_centroid(hf,
                                                normal,
                                                centroid) == -1){
                    std::cerr<<" error while computing triangle normal and centroid"<<std::endl;
                    return -1;
                }

                //skipping the zero-norm for clusters
                //if(norm(normal) == ExactType(0)){
                if(vec2vec(normal).norm() < 1e-10){
                    //std::cout<<" - 0-norm for face "<<get_halfface_vertices(hf)<<", skipping"<<std::endl;
                    continue;
                }

                double b = normal.dot(centroid).to_double() - eps;

                //-Ai^T . x - Si <= -ni^T.pi - epsilon
                lp.set_a(X, i, normal[X].to_double());
                lp.set_a(Y, i, normal[Y].to_double());
                lp.set_a(Z, i, normal[Z].to_double());

                lp.set_a(i+3, i, -1);


                lp.set_b(i, b );

                //constraint on slack variable
                lp.set_l(i+3, true, lower_bound);
                lp.set_c(i+3, w);

                i++;

                if(print_debug){
                    std::cout<<" ------------"<<std::endl;
                    std::cout<<" - halfface "<<get_halfface_vertices(hf)<<" -> ";
                    for(auto v: get_halfface_vertices(hf)){
                        std::cout<<cone_to_mesh_v_handle_prop_[v]<<" ";
                    }
                    std::cout<<": [";
                    for(auto v: get_halfface_vertices(hf)){
                        std::cout<<"("<<vec2vec(vertex(v))<<") ";
                    }
                    std::cout<<" ]"<<std::endl;

                    std::cout<<" - normal = "<<vec2vec(normal)<<
                             " (norm = "<<norm(normal).to_double()<<
                             "), centroid = "<<vec2vec(centroid)<<
                             ", rhs = "<<b<<std::endl;
                    /*std::cout<<" --> eigen norm = "<<CGAL::to_double(eigen_norm(normal))<<std::endl;
                    std::cout<<" --> quake norm = "<<CGAL::to_double(quake_norm(normal))<<std::endl;
                    std::cout<<" --> l-inf norm = "<<CGAL::to_double(linf_norm(normal))<<std::endl;*/
                }

            }else{
                /*std::cout<<" - hf "<<get_halfface_vertices(hf)<<" -> ";
                for(auto v: get_halfface_vertices(hf)){
                    std::cout<<cone_to_mesh_v_handle_prop_[v]<<" ";
                }
                std::cout<<" is not on the boundary"<<std::endl;*/
            }
        }

        //solve LP
        CGAL::Quadratic_program_options options;
        options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule
        //options.set_verbosity(1);

        //std::cout<<" - about to solve system..."<<std::endl;
        Solution s = CGAL::solve_linear_program(lp, ExactType(0), options);

        if(print_debug){
            std::cout<<" - solved system"<<std::endl;
            std::cout<<" -- is infeasible: "<<s.is_infeasible() <<std::endl;
            std::cout<<" -- is unbounded: "<<s.is_unbounded()<<std::endl;
            std::cout<<" -- obj function = "<<(s.objective_value_numerator()/s.objective_value_denominator()).to_double()<<std::endl;
        }

        std::vector<ExactType> slack_variables;
        i = 0;
        int violated_constraints_count(0);
        for(auto val_it = s.variable_values_begin(); val_it != s.variable_values_end(); val_it++){
            ExactType sol_value = (val_it->numerator()/val_it->denominator());

            if(i<3){
                new_position[i] = sol_value;
            }else{
                slack_variables.push_back(sol_value);
                //std::cout<<" s"<<i<<" = "<<sol_value<<std::endl;
                violated_constraints_count += sol_value != 0;
            }
            i++;
        }

        if(!s.solves_linear_program(lp)){

            bool found_too_large(true);
            /*for(size_t i(0); i<A.size(); i++){
                found_too_large |= (std::abs((A[i].dot(cgal_new_pos) + this->norm(A[i]) * cgal_radius) - B[i]) > 1e-7);
            }*/

            if(found_too_large){
                std::cerr<<" WARNING - couldn't solve LP with CGAL solver"<<std::endl;
                std::cout<<" error: "<<s.get_error()<<std::endl;

                return -1;
            }
        }

        return violated_constraints_count;

    }





    int ExpansionCone::find_chebyshev_center(VertexPosition& new_position,
                                             bool print_debug) const{

        new_position = {0,0,0};

        using Program  = CGAL::Quadratic_program<ExactType>;
        using Solution = CGAL::Quadratic_program_solution<ExactType>;

        //first, gather all triangles' normals and centroids
        //std::vector<VertexPosition> A;
        //std::vector<ExactType> B;

        if(false && print_debug){
            std::cout<<" gathering cone "<<*this<<" boundary normals..."<<std::endl;
            print_details();

            std::cout<<" vertices: "<<std::endl;
            for(auto v: vertices()){
                std::cout<<" - "<<v<<" -> "<<cone_to_mesh_handle(v)<<" : "<<std::setprecision(20)<<vec2vec(vertex(v))<<std::endl;
            }

            std::cout<<" cells: "<<std::endl;
            for(auto c: cells()){
                std::cout<<" - "<<c<<": "<<get_cell_vertices(c)<<std::endl;
            }
        }

        int max_precision(0);
        int n_constraints(0);

        CGAL::Quadratic_program_options options;
        options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule

        //create solver

        // by default, we have an LP with Ax <= b
        Program lp(CGAL::SMALLER, false, 0, false, 0);

        // now set the non-default entries
        const int X(0), Y(1), Z(2), R(3);

        for(auto hf: halffaces()){

            if(is_boundary(hf)){
                VertexPosition normal, centroid;

                /*std::cout<<" - computing normal for halfface "<<hf<<std::endl;
            for(auto v: get_vertices(hf)){
                std::cout<<" -- "<<v<<
                           " at "<<vec2vec(this->vertex(v))<<
                           " boundary: "<<mesh_.is_boundary(v)<<
                           " expanded: "<<expanded_prop_[v]<<std::endl;
            }*/

                if(triangle_normal_and_centroid(hf,
                                                normal,
                                                centroid) == -1){
                    std::cerr<<" error while computing triangle normal and centroid"<<std::endl;
                    return -1;
                }

                //skipping the zero-norm for clusters
                if(norm(normal) == ExactType(0)){
                    //std::cout<<" - 0-norm for face "<<get_halfface_vertices(hf)<<", skipping"<<std::endl;
                    continue;
                }

                auto b = normal.dot(centroid);


                //Ai^T . x + |Ai| . r <= Bi
                lp.set_a(X, n_constraints, normal[0]);
                lp.set_a(Y, n_constraints, normal[1]);
                lp.set_a(Z, n_constraints, normal[2]);

                lp.set_a(R, n_constraints, norm(normal));

                lp.set_b(n_constraints, b);

                n_constraints++;


                if(false && print_debug){
                    std::cout<<" ------------"<<std::endl;
                    std::cout<<" - halfface "<<get_halfface_vertices(hf)<<" -> ";
                    for(auto v: get_halfface_vertices(hf)){
                        std::cout<<cone_to_mesh_v_handle_prop_[v]<<" ";
                    }
                    std::cout<<": [";
                    for(auto v: get_halfface_vertices(hf)){
                        std::cout<<"("<<vec2vec(vertex(v))<<") ";
                    }
                    std::cout<<" ]"<<std::endl;

                    std::cout<<" - normal = "<<vec2vec(normal)<<
                             " (norm = "<<norm(normal).to_double()<<
                             "), centroid = "<<vec2vec(centroid)<<
                             ", rhs = "<<b.to_double()<<std::endl;
                    /*std::cout<<" --> eigen norm = "<<CGAL::to_double(eigen_norm(normal))<<std::endl;
                    std::cout<<" --> quake norm = "<<CGAL::to_double(quake_norm(normal))<<std::endl;
                    std::cout<<" --> l-inf norm = "<<CGAL::to_double(linf_norm(normal))<<std::endl;*/
                }

            }else{
                /*std::cout<<" - hf "<<get_halfface_vertices(hf)<<" -> ";
                for(auto v: get_halfface_vertices(hf)){
                    std::cout<<cone_to_mesh_v_handle_prop_[v]<<" ";
                }
                std::cout<<" is not on the boundary"<<std::endl;*/
            }
        }


        if(print_debug){
            std::cout<<" - added "<<n_constraints<<" constraints "<<std::endl;
        }



        //r >= 0
        ExactType lower_bound(CGAL::Gmpz(1),
                              CGAL::Gmpz(std::numeric_limits<double>::max()));
        lower_bound = lower_bound * lower_bound;
        lower_bound = 0;


        lp.set_l(R, true, lower_bound);

        if(print_debug){
            std::cout<<"   Gmpq lower bound = "<<lower_bound<<std::endl;
            std::cout<<" double lower bound = "<<CGAL::to_double(lower_bound)<<std::endl;
            std::cout<<" actual lower bound = "<<lp.get_l()[R]<<std::endl;
        }

        // Objective function: -r = (0x + 0y + 0z - 1r)
        lp.set_c(X,  0);
        lp.set_c(Y,  0);
        lp.set_c(Z,  0);
        lp.set_c(R, -1);



        //solve LP
        //options.set_verbosity(1);

        //LightWeightStopWatch sw;

        //std::cout<<" - about to solve system..."<<std::endl;
        Solution s = CGAL::solve_linear_program(lp, ExactType(), options);
        //std::cout<<" - LP solving duration: "<<sw.total_duration()<<std::endl;

        if(print_debug){
            std::cout<<" - solved system"<<std::endl;
            std::cout<<" -- is infeasible: "<<s.is_infeasible() <<std::endl;
            std::cout<<" -- is unbounded: "<<s.is_unbounded()<<std::endl;
        }

        /*ExactType first_sol =
                s.variable_values_begin()->numerator()/
                s.variable_numerators_begin()->denominator();*/


        ExactType radius;
        int i(0);
        for(auto val_it = s.variable_values_begin(); val_it != s.variable_values_end(); val_it++){
            ExactType sol_value = (val_it->numerator()/ val_it->denominator());

            //ExactType sol_value = CGAL::to_double(*val_it);

            if(false && print_debug){

                std::cout<<" -   raw sol = "<<*val_it<<std::endl;
                std::cout<<" - sol value = "<<sol_value<<std::endl;
                std::cout<<" -   raw sol double = "<<CGAL::to_double(*val_it)<<std::endl;
                std::cout<<" - sol value double = "<<sol_value.to_double()<<std::endl;
                std::cout<<" ---------"<<std::endl;
            }

            if(i<3){
                new_position[i] = sol_value;
            }else{
                radius = sol_value;
                if(false && print_debug){

                    std::cout<<" radius = "<<radius<<std::endl;
                    std::cout<<" radius double = "<<radius.to_double()<<std::endl;
                }
            }
            i++;
        }
        //std::cout<<" -------------------------"<<std::endl;



        if(print_debug){
            std::cout<<" - updated new_position"<<std::endl;
        }

        /*if(s.is_infeasible() && radius == chebyshev_radius_lower_bound){
            std::cerr<<" ERROR - reached lower chebyshev radius lower bound"<<std::endl;
            return -1;
        }*/




#warning fix this for double representation
        if(!s.solves_linear_program(lp)){

            std::cerr<<" WARNING - couldn't solve LP with CGAL solver"<<std::endl;
            std::cout<<" error: "<<s.get_error()<<std::endl;

            return 1;
        }



        /*std::cout<<" lower bound = "<<lower_bound.to_double()<<std::endl;
        std::cout<<" Chebyshev center at = "<<vec2vec(new_position)<<std::endl;
        std::cout<<"              radius = "<<radius.to_double()<<std::endl;*/

        //std::cout<<" infeasible: "<<s.is_infeasible()<<std::endl;

        auto double_precision_position = vec2vec(new_position);
        //using the l-inf norm for safety (0.7 as an approximation of sqrt(2)/2)
        if(linf_norm(new_position - vec2vec(double_precision_position)) < radius*0.7){
            new_position = vec2vec(double_precision_position);
            //std::cout<<" -----> double precision Chebyshev center is enough"<<std::endl;
        }

        if(s.is_infeasible() || s.is_unbounded() || !s.is_valid()){
            //std::cout<<" --> infeasible: "<<s.is_infeasible()<<std::endl;
            //std::cout<<" -->  unbounded: "<<s.is_unbounded()<<std::endl;
            //std::cout<<" -->      valid: "<<s.is_valid()<<std::endl;
            return 1;
        }

#if ENABLE_EXACT_REPRESENTATION
        if(radius == 0){
            //std::cout<<" --> zero max radius"<<std::endl;
            return 1;
        }

        if(new_position == VertexPosition(0,0,0)){
            std::cout<<" --> feasible point at origin"<<std::endl;
            return 1;
        }

#else
        if(std::abs(cgal_radius) < radius_eps_){
        //std::cout<<" - radius "<<r->solution_value()<<" is smaller than epsilon "<<radius_eps_<<", unvalid result"<<std::endl;
        return 1;
    }

    if(this->norm(new_position) < (distance_eps_ * distance_eps_)){
        //std::cout<<" - Chebyshev center at "<<new_position<<" has norm = "<<new_position.norm()<<", which is too close to the center"<<std::endl;
        return 1;
    }
#endif


        //std::cout<<" ---->  problem max size = "<<max_precision<<std::endl;
        //std::cout<<" ----> solution max size = "<<max_size(new_position)<<std::endl;
        return 0;

    }



    int ExpansionCone::find_chebyshev_center_from_halffaces(const std::vector<std::vector<VertexHandle>>& halffaces,
                                                            VertexPosition& new_position,
                                                            bool print_debug) const{

        using Program  = CGAL::Quadratic_program<ExactType>;
        using Solution = CGAL::Quadratic_program_solution<ExactType>;

        //first, gather all triangles' normals and centroids
        std::vector<VertexPosition> A;
        std::vector<ExactType> B;

        int max_precision(0);

        for(auto& hf: halffaces){

            VertexPosition normal, centroid;

            if(triangle_normal_and_centroid(hf,
                                            normal,
                                            centroid) == -1){
                std::cerr<<" error while computing triangle normal and centroid"<<std::endl;
                return -1;
            }

            //skipping the zero-norm for clusters
            if(norm(normal) == ExactType(0)){
                //std::cout<<" - 0-norm for face "<<get_halfface_vertices(hf)<<", skipping"<<std::endl;
                continue;
            }

            A.push_back(normal);
            B.push_back(normal.dot(centroid));

            if(print_debug){
                std::cout<<" ------------"<<std::endl;
                std::cout<<" - halfface "<<hf<<" -> ";
                for(auto v: hf){
                    std::cout<<cone_to_mesh_v_handle_prop_[v]<<" ";
                }
                std::cout<<": [";
                for(auto v: hf){
                    std::cout<<"("<<vec2vec(vertex(v))<<") ";
                }
                std::cout<<" ]"<<std::endl;

                std::cout<<" - normal = "<<vec2vec(normal)<<
                         " (norm = "<<norm(normal).to_double()<<
                         "), centroid = "<<vec2vec(centroid)<<
                         ", rhs = "<<B.back().to_double()<<std::endl;
                /*std::cout<<" --> eigen norm = "<<CGAL::to_double(eigen_norm(normal))<<std::endl;
                std::cout<<" --> quake norm = "<<CGAL::to_double(quake_norm(normal))<<std::endl;
                std::cout<<" --> l-inf norm = "<<CGAL::to_double(linf_norm(normal))<<std::endl;*/
            }

            max_precision = std::max(max_precision, max_byte_size(A.back()));
            max_precision = std::max(max_precision, (int)B.back().size());
        }


        //create solver

        // by default, we have an LP with Ax <= b
        Program lp(CGAL::SMALLER, false, 0, false, 0);

        // now set the non-default entries
        const int X(0), Y(1), Z(2), R(3);

        //std::cout<<" - set-up lp, A size = "<<A.size()<<std::endl;

        //set-up constraints
        for(size_t i(0); i<A.size(); i++){
            //Ai^T . x + |Ai| . r <= Bi
            lp.set_a(X, i, A[i][0]);
            lp.set_a(Y, i, A[i][1]);
            lp.set_a(Z, i, A[i][2]);

            lp.set_a(R, i, norm(A[i]));

            lp.set_b(i, B[i]);
        }



        //r >= 0
        ExactType lower_bound(CGAL::Gmpz(1),
                              CGAL::Gmpz(std::numeric_limits<double>::max()));
        lower_bound = lower_bound * lower_bound;
        lower_bound = 0;


        lp.set_l(R, true, lower_bound);

        // Objective function: -r = (0x + 0y + 0z - 1r)
        lp.set_c(X,  0);
        lp.set_c(Y,  0);
        lp.set_c(Z,  0);
        lp.set_c(R, -1);

        //solve LP
        CGAL::Quadratic_program_options options;
        options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule

        //options.set_verbosity(1);

        //std::cout<<" - about to solve system..."<<std::endl;
        Solution s = CGAL::solve_linear_program(lp, ExactType(), options);

        if(print_debug){
            std::cout<<" - solved system"<<std::endl;
            std::cout<<" -- is infeasible: "<<s.is_infeasible() <<std::endl;
            std::cout<<" -- is unbounded: "<<s.is_unbounded()<<std::endl;
        }

        ExactType radius;
        int i(0);
        for(auto val_it = s.variable_values_begin(); val_it != s.variable_values_end(); val_it++){
            ExactType sol_value = (val_it->numerator()/val_it->denominator());
            if(i<3){
                new_position[i] = sol_value;
            }else{
                radius = sol_value;
                //std::cout<<" radius size = "<<radius.size()<<std::endl;
            }
            i++;
        }
        //std::cout<<" -------------------------"<<std::endl;

#warning fix this for double representation
        if(!s.solves_linear_program(lp)){

            bool found_too_large(true);
            /*for(size_t i(0); i<A.size(); i++){
                found_too_large |= (std::abs((A[i].dot(cgal_new_pos) + this->norm(A[i]) * cgal_radius) - B[i]) > 1e-7);
            }*/

            if(found_too_large){
                std::cerr<<" WARNING - couldn't solve LP with CGAL solver"<<std::endl;
                std::cout<<" error: "<<s.get_error()<<std::endl;

                for(size_t i(0); i<A.size(); i++){
                    std::cout<<" - A("<<i<<")t . x + |A("<<i<<")| . r = "<<(A[i].dot(new_position) + norm(A[i]) * radius) <<
                             ", b("<<i<<") = "<<B[i]<<std::endl;
                }
                return 1;
            }
        }



        /*std::cout<<" lower bound = "<<lower_bound.to_double()<<std::endl;
        std::cout<<" Chebyshev center at = "<<vec2vec(new_position)<<std::endl;
        std::cout<<"              radius = "<<radius.to_double()<<std::endl;*/

        //std::cout<<" infeasible: "<<s.is_infeasible()<<std::endl;

        auto double_precision_position = vec2vec(new_position);
        //using the l-inf norm for safety (0.7 as an approximation of sqrt(2)/2)
        if(linf_norm(new_position - vec2vec(double_precision_position)) < radius*0.7){
            new_position = vec2vec(double_precision_position);
            //std::cout<<" -----> double precision Chebyshev center is enough"<<std::endl;
        }

        if(s.is_infeasible() || s.is_unbounded() || !s.is_valid()){
            return 1;
        }

#if ENABLE_EXACT_REPRESENTATION
        if(radius == 0){
            std::cout<<" --> zero max radius"<<std::endl;
            return 1;
        }

        if(new_position == VertexPosition(0,0,0)){
            std::cout<<" --> feasible point at origin"<<std::endl;
            return 1;
        }

#else
        if(std::abs(cgal_radius) < radius_eps_){
        //std::cout<<" - radius "<<r->solution_value()<<" is smaller than epsilon "<<radius_eps_<<", unvalid result"<<std::endl;
        return 1;
    }

    if(this->norm(new_position) < (distance_eps_ * distance_eps_)){
        //std::cout<<" - Chebyshev center at "<<new_position<<" has norm = "<<new_position.norm()<<", which is too close to the center"<<std::endl;
        return 1;
    }
#endif


        //std::cout<<" ---->  problem max size = "<<max_precision<<std::endl;
        //std::cout<<" ----> solution max size = "<<max_size(new_position)<<std::endl;
        return 0;

    }


    int ExpansionCone::triangle_normal_and_centroid(const HalfFaceHandle& face,
                                                    VertexPosition& normal,
                                                    VertexPosition& centroid) const{

        return triangle_normal_and_centroid(get_halfface_vertices(face), normal, centroid);
    }

    int ExpansionCone::triangle_normal_and_centroid(const std::vector<VertexHandle>& face_vertices,
                                                    VertexPosition& normal,
                                                    VertexPosition& centroid) const{

        auto c_evi  = this->vertex(face_vertices[1]) - this->vertex(face_vertices[0]);
        auto c_evi1 = this->vertex(face_vertices[2]) - this->vertex(face_vertices[0]);


        normal = (c_evi.cross(c_evi1));
        /*if(norm(normal) == 0){
            std::cerr<<" ERROR - found zero-norm normal"<<std::endl;
            //return -1;
        }*/
        normalize(normal);

        centroid = (this->vertex(face_vertices[0]) +
                    this->vertex(face_vertices[1]) +
                    this->vertex(face_vertices[2]))/3.f;

        return 0;
    }


    bool ExpansionCone::is_visible_from(const VertexHandle& cone_from_vertex,
                                        const VertexHandle& cone_to_vertex) const{


        //std::cout<<" ================================================ intersections with segment "<<cone_from_vertex<<"-"<<cone_to_vertex<<std::endl;

        //first, check if the vertices are neighbors
        for(auto out_he: outgoing_halfedges(cone_from_vertex)){
            if(to_vertex_handle(out_he) == cone_to_vertex){
                //std::cout<<" --> neighbors -> visible"<<std::endl;
                return true;
            }
        }

        CGAL::Segment_3<CGAL_ExactKernel> segment(OVMvec3ToCGALPoint3(vertex(cone_from_vertex)),
                                                  OVMvec3ToCGALPoint3(vertex(cone_to_vertex)));

        const int scale_factor(std::numeric_limits<int>::max());

        bool intersected_with_interior_triangle(false);

        for(auto f: faces()){

            auto hf = halfface_handle(f, 0);
            auto hf_vertices = get_halfface_vertices(hf);

            //skipping triangles not containing any tip vertex
            if(!tip_vertices_prop_[hf_vertices[0]] &&
                    !tip_vertices_prop_[hf_vertices[1]] &&
                    !tip_vertices_prop_[hf_vertices[2]]){
                continue;
            }

            //skipping triangles containing either the start or end vertex
            //because the intersection with those can be misleading.
            //e.g. a segment completely outside of the cone, connecting a vertex on
            //a bunch of triangles will intersect all those triangles.
            //This might make it look like the segment passes through an interior
            //triangle, i.e. that it's an interior segment although it isn't
            if(TopoHelper::face_contains_vertex(*this, cone_from_vertex, f) ||
                    TopoHelper::face_contains_vertex(*this, cone_to_vertex, f)){
                continue;
            }

            CGAL::Triangle_3<CGAL_ExactKernel> tri(
                        OVMvec3ToCGALPoint3(vertex(hf_vertices[0]), scale_factor),
                    OVMvec3ToCGALPoint3(vertex(hf_vertices[1]), scale_factor),
                    OVMvec3ToCGALPoint3(vertex(hf_vertices[2]), scale_factor));


            /*std::cout<<" - checking intersection with tri "<<hf_vertices<<": ("
                     <<tri.vertex(0)<<"), ("
                     <<tri.vertex(1)<<"), ("
                     <<tri.vertex(2)<<")"<<std::endl;*/


            auto inter = CGAL::intersection(segment, tri);

            if(inter){
                //std::cout<<" - intersection with tri "<<hf_vertices<<": "<<inter<<std::endl;

                //this test makes the visibility compatible with zero-volume tets.
                //with such tets, you can have a vertex lying exactly on an interior face.
                //in that case, there's an intersection with this interior face, although the
                //segment is coming from the outside.
                //The problem is that then, this is indistinguishable from a segment reaching
                //this point from the *inside*
                //Tricky stuff, I know.
                /*auto* point_result = boost::get<CGAL_ExactPoint3>(&*inter);
                if(point_result){
                    std::cout<<" -- intersection is a point"<<std::endl;
                    if(*point_result == OVMvec3ToCGALPoint3(vertex(cone_to_vertex))){
                        std::cout<<" ---> which is the to-vertex"<<std::endl;
                        //continue;
                    }
                }*/

                if(!is_boundary(f)){
                    //std::cout<<" ----> intersected with interior triangle"<<std::endl;
                    intersected_with_interior_triangle = true;
                }else{
                    //std::cout<<" ----> intersection with a boundary face that doesn't contain either vertex -> not visible"<<std::endl;
                    return false;
                }
            }
        }

        //std::cout<<" ------> visible: "<<intersected_with_interior_triangle<<std::endl;


        return intersected_with_interior_triangle;

    }

    VertexPosition ExpansionCone::north_pole() const{
        VertexPosition N(0,0,0);
        for(auto vhf_it = vhf_iter(VertexHandle(0)); vhf_it.valid(); vhf_it++){
            if(is_boundary(*vhf_it)){
                VertexPosition normal, centroid;
                triangle_normal_and_centroid(*vhf_it, normal, centroid);
                N += normal;
            }
        }

        normalize(N);
        return N;
    }



    void ExpansionCone::project_base_to_sphere(){

        for(auto v: vertices()){
            if(!tip_vertices_prop_[v]){
                auto pos = this->vertex(v);
                auto pos_norm = norm(pos);
                if(pos_norm == CGAL::Gmpq(0)){
                    std::cout<<" --> can't project vertex "<<v<<" because it's at the origin"<<std::endl;
                    continue;
                }
                pos /= pos_norm;
                this->set_vertex(v, pos);
            }
        }
    }


    void ExpansionCone::project_base_stereographically(const VertexPosition& north_pole){
        project_base_to_sphere();
        const auto& N = north_pole;
        const auto& Nz = N[2];
        const VertexPosition target_N(0,0,-1);

        //std::cout<<" - north pole = "<<vec2vec(north_pole)<<std::endl;

        //to rotate projection plane to horizontal plane with N = (0,0,-1)
        //this makes the cone base face the Z axis without the tip vertex obstructing the field

        //angle = acos(N.(0,0,-1)) = acos(-N_z)
        auto theta = std::acos(-Nz.to_double());

        //axis computation
        VertexPosition A;
        //shouldn't happen, but just in case
        if(N[2] == CGAL::Gmpq(1)){
            A = {1,0,0};
            theta = 0;
        }else{
            A = N.cross(target_N);
        }
        normalize(A);
        //std::cout<<" - axis = "<<vec2vec(A)<<", angle = "<<theta<<std::endl;
        //std::cout<<" - A . N = "<<A.dot(N)<<std::endl;


        for(auto v: vertices()){
            if(!tip_vertices_prop_[v]){
                const auto& P = this->vertex(v);
                //std::cout<<" - moved vertex "<<v<<" from "<<vec2vec(P);

                //stereographic projection
                const auto P_p = N + 2*(P-N)/(1-N.dot(P));

                //axis-angle rotation
                const auto P_rot = -Nz*P_p + std::sin(theta)*(A.cross(P_p)) + (1+Nz)*A.dot(P_p)*A;

                this-> set_vertex(v, P_rot);
                //std::cout<<" to "<<vec2vec(vertex(v))<<std::endl;
            }
        }
    }


    void ExpansionCone::translate_to_move_tip_to_origin(){
        auto tip_pos = vertex(VertexHandle(0));
        for(auto v: vertices()){
            set_vertex(v, vertex(v) - tip_pos);
        }
    }


    int ExpansionCone::total_position_byte_size() const{
        int total_size_B(0);
        for(auto v: vertices()){
            total_size_B += byte_size(vertex(v));
        }
        return total_size_B;
    }


    void ExpansionCone::print_details(int detail_level) const{

        if(detail_level){
            std::cout<<" ================ cone "<<*this<<std::endl;
            std::cout<<" #### VERTICES: "<<std::endl;
            for(auto v: vertices()){
                std::cout<<" - "<<v<<
                           " -> "<<std::setw(5)<<cone_to_mesh_v_handle_prop_[v]<<
                           " at "<<vec2vec(vertex(v))<<
                           (tip_vertices_prop_[v] ? " (tip)" : "")<<
                           (is_boundary(v) ? " (boundary)" : "")<<std::endl;
            }
            std::cout<<" # cone tips: "<<cone_tip_vertices()<<std::endl;

            if(detail_level > 1){

                std::cout<<" #### CELLS: "<<std::endl;
                for(auto c: cells()){
                    std::cout<<" - "<<c<<": "<<get_cell_vertices(c)<<" -> (";
                    for(auto v: get_cell_vertices(c)){
                        std::cout<<cone_to_mesh_handle(v)<<" ";
                    }
                    std::cout<<" )"<<std::endl;
                }

            }

            if(detail_level > 2){
                std::cout<<" #### FACES: "<<std::endl;
                for(auto f: faces()){
                    auto hf = halfface_handle(f,0);
                    std::cout<<" - "<<f<<": "<<get_halfface_vertices(hf)<<" -> (";
                    for(auto v: get_halfface_vertices(hf)){
                        std::cout<<cone_to_mesh_handle(v)<<" ";
                    }
                    std::cout<<" )"<<std::endl;
                }
            }

            if(detail_level > 3){
                std::cout<<" #### EDGES: "<<std::endl;
                for(auto e: edges()){
                    std::cout<<" - "<<e<<": "<<edge(e)<<std::endl;
                }
            }

            /*
            if(detail_level > 3){
                std::cout<<" #### HALFFACES: "<<std::endl;
                for(auto hf: halffaces()){
                    std::cout<<" - "<<hf<<": "<<get_halfface_vertices(hf)<<" -> (";
                    for(auto v: get_halfface_vertices(hf)){
                        std::cout<<cone_to_mesh_handle(v)<<" ";
                    }
                    std::cout<<" )"<<std::endl;
                }
            }*/

            std::cout<<" ====================================== "<<std::endl;
        }

    }

    bool ExpansionCone::base_is_a_single_connected_component(){
        auto visited_prop = request_vertex_property<bool>();

        auto center_vertex = VertexHandle(0);
        visited_prop[center_vertex] = true;

        if(!n_cells() || n_vertices() < 2){
            return true;
        }

        std::vector<VertexHandle> to_visit = {VertexHandle(1)};

        while(!to_visit.empty()){
            auto current_vertex = to_visit.back();
            to_visit.pop_back();

            if(visited_prop[current_vertex]){
                continue;
            }

            visited_prop[current_vertex] = true;

            for(auto out_he_it: outgoing_halfedges(current_vertex)){
                auto neighbor = to_vertex_handle(out_he_it);
                if(!visited_prop[neighbor]){
                    to_visit.push_back(neighbor);
                }
            }
        }

        for(auto v: vertices()){
            if(!visited_prop[v]){
                return false;
            }
        }

        return true;
    }


    int ExpansionCone::tips_euler_characteristic() const{


        TetrahedralMesh cone_tips = static_cast<TetrahedralMesh>(*this);

        for(auto v: vertices()){
            if(!tip_vertices_prop_[v]){
                cone_tips.delete_vertex(v);
             }
        }
        cone_tips.collect_garbage();


        /*if(true){
            std::cout<<" - cone tips mesh: "<<std::endl;
            std::cout<<" VERTICES:"<<std::endl;
            for(auto v: cone_tips.vertices()){
                std::cout<<" -- "<<v<<"->"<<cone_to_mesh_handle(v)<<std::endl;
            }
            std::cout<<" FACES:"<<std::endl;
            for(auto f: cone_tips.faces()){
                std::cout<<" -- "<<f<<": "<<cone_tips.get_halfface_vertices(cone_tips.halfface_handle(f, 0))<<std::endl;
            }
        }*/

        if(cone_tips.n_logical_vertices() != cone_tip_vertices().size()){
            std::cout<<" ERROR - cone tips mesh contains "<<cone_tips.n_logical_vertices()<<
                       " vertices but there are "<<cone_tip_vertices().size()<<" cone tip vertices"<<std::endl;
            std::cout<<" cone tips : "<<cone_tip_vertices()<<std::endl;
            return -1;
        }

        /*std::cout<<" - cone tips verts count = "<<cone_tip_vertices_.size()<<std::endl;
        std::cout<<" - cone tips edges count = "<<cone_tips.n_logical_edges()<<std::endl;
        std::cout<<" - cone tips faces count = "<<cone_tips.n_logical_faces()<<std::endl;
        std::cout<<" - cone tips cells count = "<<cone_tips.n_logical_cells()<<std::endl;
*/

        return (int)cone_tip_vertices_.size() -
                (int)cone_tips.n_logical_edges() +
                (int)cone_tips.n_logical_faces() -
                (int)cone_tips.n_logical_cells();

    }


    EdgeHandle ExpansionCone::opposite_edge_on_face(const VertexHandle& vh,
                                                    const FaceHandle& fh) const{

        //auto hfh = halfface_handle(fh, 0);

        //return next_halfedge_in_halfface(next_halfedge_in_halfface(halfedge() , hfh), hfh);

        for(auto fe_it = fe_iter(fh); fe_it.valid(); fe_it++){
            if(edge(*fe_it).from_vertex() != vh &&
                    edge(*fe_it).to_vertex() != vh){
                return *fe_it;
            }
        }
        return EdgeHandle(-1);
    }



    int ExpansionCone::cell_valence(const VertexHandle& v) const{

        int val(0);
        for(auto vc_it = vc_iter(v); vc_it.valid(); vc_it++){
            val++;
            //std::cout<<" -- cell "<<(*vc_it)<<": "<<get_cell_vertices(*vc_it)<<std::endl;
        }
        return val;
    }


    int ExpansionCone::cell_valence(const EdgeHandle& e) const{

        int val(0);
        for(auto ec_it = ec_iter(e); ec_it.valid(); ec_it++){
            val++;
        }
        return val;
    }



    std::ostream& operator<<(std::ostream& os,
                             const ExpansionCone& cone){
        os<<" EC(V,E,F,C) = ("<<
          cone.n_logical_vertices()<<", "<<
          cone.n_logical_edges()<<", "<<
          cone.n_logical_faces()<<", "<<
          cone.n_logical_cells()<<"), tips (cone->mesh): ";
        for(auto cone_tip_vertex: cone.cone_tip_vertices()){
            os<<cone_tip_vertex<<"->"<<cone.cone_to_mesh_handle(cone_tip_vertex)<<" ";
        }

        return os;

    }


#if ENABLE_EXACT_REPRESENTATION
    VertexPosition vec2vec(const Vec3d& v){
        return {v[0], v[1], v[2]};
    }

    Vec3d vec2vec(const VertexPosition& v){
        return {v[0].to_double(), v[1].to_double(), v[2].to_double()};
    }
#else
    VertexPosition vec2vec(const VertexPosition& v){
    return v;
}
#endif

    ExactType norm(const VertexPosition& pos){

        //std::cout<<" quake norm = "<<CGAL::to_double(quake_norm(pos))<<std::endl;
        //std::cout<<"  linf norm = "<<CGAL::to_double(linf_norm(pos))<<std::endl;

        //NOTE: eigen norm is standard
        return eigen_norm(pos);
        //return halley_norm(pos);
        return l1_norm(pos);
        //return quake_norm(pos);
        return linf_norm(pos);
    }

    ExactType halley_norm(const VertexPosition& pos){

        //std::cout<<" epsilon = "<<eps<<std::endl;

        ExactType a = pos.sqrnorm();
        ExactType x_0 = a;

        //std::cout<<" --> dbl norm = "<<dbl_norm<<std::endl;

        ExactType x_n = x_0;

        //should be enough and if not, it'll take too long
        const int iteration_number(5);

        for(int i(0); i<iteration_number; i++){

            //std::cout<<" ------- iteration "<<i<<std::endl;

            ExactType x_n2 = x_n * x_n;
            ExactType x_n1 = x_n * (3*a + x_n2)/(a + 3*x_n2);

            /*std::cout<<"      x_n = "<<x_n.to_double()<<std::endl;
            std::cout<<"    x_n+1 = "<<x_n1.to_double()<<std::endl;
            std::cout<<" sqr norm = "<<a.to_double()<<std::endl;
            std::cout<<" dbl diff = "<<CGAL::abs(a - x_n1*x_n1).to_double()<<std::endl;*/

            x_n = x_n1;
        }

        return x_n;
    }

    ExactType l1_norm(const VertexPosition& pos){
        return abs(pos[0]) + abs(pos[1]) + abs(pos[2]);
    }

    ExactType linf_norm(const VertexPosition& pos){
        return max(abs(pos.max()), abs(pos.min()));
    }


/* Code taken from wikipedia: https://en.wikipedia.org/wiki/Fast_inverse_square_root
 * Original Quake III source code comments kept for posterity */
    ExactType quake_norm(const VertexPosition& pos){
        ExactType sqr_nrm = pos.sqrnorm();

        //std::cout<<" quake sqrnorm = "<<sqr_nrm<<std::endl;
        float x = CGAL::to_double(sqr_nrm);
        long i;
        float x2, y;
        const float threehalfs = 1.5F;

        x2 = x * 0.5F;
        y  = x;
        i  = * ( long * ) &y;                       // evil floating point bit level hacking
        i  = 0x5f3759df - ( i >> 1 );               // what the fuck?
        y  = * ( float * ) &i;
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 3rd bonus iteration

        //std::cout<<" y = "<<y<<std::endl;

        ExactType test = 1.0/y;

        return test;
    }


    ExactType eigen_norm(const VertexPosition& pos){

        return Eigen::Vector3d(pos[0].to_double(),
                               pos[1].to_double(),
                               pos[2].to_double()).norm();

    }



    VertexPosition normalize(VertexPosition& pos){

        auto norm_ = norm(pos);
        if(norm_ != 0){
            pos /= norm_;
        }else{

            //std::cout<<" --> FOUND ZERO NORM FOR POS "<<pos<<std::endl;
            //exit(EXIT_FAILURE);

            //std::cout<<" using l-inf norm"<<std::endl;

            auto alt_norm = l1_norm(pos);
            if(alt_norm == 0){
            }else{
                pos /= alt_norm;
            }
        }
        return pos;
    }



    std::ostream& operator<<(std::ostream& ss, const Split& split){

        ss<<"("<<split.cone_from_vertex<<"-"<<split.cone_to_vertex<<") -> "<<split.cone_new_vertex<<" at "<<vec2vec(split.new_vertex_position);

        return ss;
    }
}
