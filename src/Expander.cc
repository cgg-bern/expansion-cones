#include "Expander.hh"
#include <random>



namespace OpenVolumeMesh{


#define RAND_SEED 1234
#define DEFAULT_CLUSTER_SIZE_HARD_LIMIT 10
#define DEFAULT_MAX_EXPANSION_VALENCE 20


#define STOP_AT_SIMPLE_EXPANSION 0
#define MANUAL_STOP_ITERATION -1

#define DEFAULT_INTERIOR_SMOOTHING_RING_K 1
#define DEFAULT_INTERIOR_SMOOTHING_MAX_ITERATIONS 100
#define DEFAULT_INTERIOR_SMOOTHING_EPSILON 1
#define DEFAULT_POSITION_BYTE_SIZE_THRESHOLD 400
#define DEFAULT_POSITION_BYTE_SIZE_SMOOTHING_THRESHOLD 100
#define SMOOTHING_ENABLED 1

#define EDGE_COLLAPSE_ENABLED 1
//0: Chebyshev-centroid, 1: some point on the edge
#define EDGE_COLLAPSE_MODE 0
//-1 = as many as there are candidates
#define EDGE_COLLAPSE_MAX_EDGES_TO_COLLAPSE_COUNT 10000
#define EDGE_COLLAPSE_MAX_ITERATION_COUNT 1
#define EDGE_COLLAPSE_ITERATION_TIMEOUT_S (10*60)
#define EDGE_COLLAPSE_MAX_STEP_COUNT 3

#define DEFAULT_MAX_GROWTH_RATIO 50

#define CLUSTER_SS_ENABLED 1

#define ENABLE_ALL_CHECKS 0

#define ENABLE_FINAL_CHECKS 1


#define PRINT_IF_NOT_SILENT(msg) if(!silent_mode_) {std::cout<<msg;}


    Expander::Expander(TetrahedralMesh& mesh,
                       const TetrahedralMesh& input_mesh,
                       const std::string& mesh_name,
                       bool silent_mode,
                       VertexPropertyT<VertexPosition>* exact_vertices_positions,
                       bool expanding_cluster,
                       bool allow_flipped_initial_flipped_tets,
                       const int max_allocated_time_s) :
            initial_n_vertices_(mesh.n_vertices()),
            //initial_non_link_edges_count_(TopoHelper::findNonLinkEdges(mesh).size()),
            global_start_time_(std::chrono::high_resolution_clock::now()),
            max_allocated_time_s_(max_allocated_time_s),
            mesh_(mesh),
            expanded_prop_(mesh.request_vertex_property<bool>("expanded_vertex")),
            vertex_position_prop_(mesh_.request_vertex_property<VertexPosition>("exact_vertex_position")),
            topo_helper_(mesh),
            expanding_cluster_(expanding_cluster),
            domain_vertex_position_prop_(mesh_.request_vertex_property<VertexPosition>("domain_vertex_position")),
            expanse_at_last_iteration_prop_(mesh_.request_vertex_property<VertexExpanse>()),
            unexpanded_neighbors_at_expanse_computation_(mesh_.request_vertex_property<std::vector<VertexHandle>>()),
            moved_during_smoothing_prop_(mesh_.request_vertex_property<bool>()),
            valence_at_expanse_computation_(mesh_.request_vertex_property<std::pair<int,int>>()),
            collapsible_new_halfedge_prop_(mesh_.request_halfedge_property<bool>("collapsible_new_edge")),
            //full_cone_at_last_iteration_prop_(codomain_mesh_.request_vertex_property<std::pair<ExpansionCone, int>>("last_expanse")),
            cluster_index_prop_(mesh_.request_vertex_property<int>()),
            connecting_cluster_indices_prop_(mesh_.request_edge_property<std::pair<int,int>>()),
            data_logger_(mesh_name),
            silent_mode_(silent_mode){

        PRINT_IF_NOT_SILENT(" ===================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" ===================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" INITIALIZING EXPANDER"<<(expanding_cluster_ ? " FOR CLUSTER" : "")<<"..."<<std::endl);


        std::srand(RAND_SEED);

        bool use_exact_pos_from_file(false);

        bool expanded_vertices_already_set(false);
        if(at_least_one_expanded_vertex(mesh_)){

            expanded_vertices_already_set = true;
            PRINT_IF_NOT_SILENT(" ---> expanded prop was already set"<<std::endl);
            PRINT_IF_NOT_SILENT(" ---> USING EXACT POSITIONS FROM FILE"<<std::endl);
            use_exact_pos_from_file = true;
        }

        if(exact_vertices_positions &&
           !use_exact_pos_from_file){
            PRINT_IF_NOT_SILENT(" ---> USING EXACT POSITIONS FROM PROPERTY POINTER "<<exact_vertices_positions<<std::endl);
        }

        auto str_exact_vertices_positions_prop = mesh_.request_vertex_property<std::string>("str_exact_position");

        int max_exact_pos_byte_size(0);

        for(size_t i(0); i<mesh_.n_vertices(); i++){
            const auto v = VertexHandle(i);

            bool is_boundary = mesh_.is_boundary(v);
            //expanded_prop_[v] = is_boundary;
            if(!expanded_vertices_already_set){
                set_expanded_prop(v, is_boundary);
            }

            to_expand_count_ += !expanded_prop_[v];
            cluster_index_prop_[v] = expanded_prop_[v] ? -1 : 0;

            if(use_exact_pos_from_file){
                this->set_vertex(v, PEHelpers::string_to_position(str_exact_vertices_positions_prop[v]));
                max_exact_pos_byte_size = std::max(max_exact_pos_byte_size,
                                                   max_byte_size(this->vertex(v)));

                //std::cout<<" - "<<v<<" is deleted: "<<codomain_mesh_.is_deleted(v)<<", pos size = "<<std::setw(6)<<byte_size(vertex(v))<<" at "<<std::setprecision(30)<<vec2vec(vertex(v))<<std::endl);


            }else{

                if(exact_vertices_positions){
                    this->set_vertex(v, (*exact_vertices_positions)[v]);
                }else{
                    this->set_vertex(v, vec2vec(mesh_.vertex(v)));
                }
            }

        }

        if(use_exact_pos_from_file){
            PRINT_IF_NOT_SILENT(" initial positions max byte size = "<<max_exact_pos_byte_size<<std::endl);
        }

        mesh_.set_shared(expanded_prop_);
        mesh_.set_persistent(expanded_prop_);
        //codomain_mesh_.set_persistent(vertex_position_prop_);
        mesh_.set_shared(collapsible_new_halfedge_prop_);
        mesh_.set_persistent(collapsible_new_halfedge_prop_);
        //codomain_mesh_.set_shared(expanse_at_last_iteration_prop_);
        //codomain_mesh_.set_persistent(expanse_at_last_iteration_prop_);


        for(auto e: mesh_.edges()){
            auto from_v = mesh_.from_vertex_handle(mesh_.halfedge_handle(e, 0));
            auto to_v = mesh_.to_vertex_handle(mesh_.halfedge_handle(e, 0));

            connecting_cluster_indices_prop_[e] = {mesh_.is_boundary(from_v) ? -1 : 0,
                                                   mesh_.is_boundary(to_v)   ? -1 : 0};
        }


        auto bad_tets = ExactBadTetFinder::findBadTets(mesh_,
                                                       vertex_position_prop_);

        PRINT_IF_NOT_SILENT("- bad tets at expander initialization:"<<std::endl);
        PRINT_IF_NOT_SILENT("    "<<bad_tets.first.size()<<" degenerate tets"<<std::endl);
        PRINT_IF_NOT_SILENT("    "<<bad_tets.second.size()<<" flipped tets"<<std::endl);

        if(!bad_tets.second.empty() && !allow_flipped_initial_flipped_tets){
            std::cout<<" ERROR - INITIAL MESH CONTAINS FLIPPED TETS: "<<std::endl;
            for(auto flipped_tet: bad_tets.second){
                PRINT_IF_NOT_SILENT(" - "<<flipped_tet<<": "<<mesh_.get_cell_vertices(flipped_tet)<<std::endl);
            }
            std::cout<<" EXITING NOW"<<std::endl;
            exit(EXIT_FAILURE);
        }

        PRINT_IF_NOT_SILENT(" - vertices to expand at initialization: "<<to_expand_count_<<std::endl);

        if(!to_expand_count_){
            PRINT_IF_NOT_SILENT(" ---> mesh is fully expanded at initialization"<<std::endl);
            fully_expanded_ = true;
        }

        bool cluster_single_cc(false);
        int cluster_euler_charac(0);
        check_main_cluster_topology(cluster_single_cc, cluster_euler_charac);

        PRINT_IF_NOT_SILENT(" - main cluster has a single connected component: "<<cluster_single_cc<<std::endl);
        PRINT_IF_NOT_SILENT(" - main cluster has euler characteristic: "<<cluster_euler_charac<<std::endl);


        if(mesh_.n_vertices() != input_mesh.n_vertices()){
            std::cout<<" ERROR - input meshes don't have the same number of vertices"<<std::endl;
            exit(EXIT_FAILURE);
        }

        for(auto v: mesh_.vertices()){
            domain_vertex_position_prop_[v] = vec2vec(input_mesh.vertex(v));
        }

#if ENABLE_ALL_CHECKS

        if(!expanding_cluster_) {

            auto domain_bad_tets = ExactBadTetFinder::findBadTets(mesh_, domain_vertex_position_prop_);


            if (!domain_bad_tets.first.empty()) {
                std::cout << " ERROR - initial domain mesh contains degenerate tets" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (!domain_bad_tets.second.empty()) {
                std::cout << " ERROR - initial domain mesh contains flipped tets " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
#endif

        PRINT_IF_NOT_SILENT(" ===================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" ===================================="<<std::endl);
    }


    bool Expander::is_fully_expanded() const{
        return fully_expanded_;
    }


    int Expander::to_expand_count() const{
        return to_expand_count_;
    }

    int Expander::left_to_expand_count() const {
        return to_expand_count() - expanded_count();
    }

    int Expander::expanded_count() const{
        return expanded_count_;
    }


    int Expander::star_shapification_count() const{
        return star_shapification_successes_count_;
    }

    int Expander::cluster_expansion_count() const{
        int count(0);
        for(auto c: cluster_size_count_){
            count+=c;
        }
        return count;
    }

    int Expander::cluster_star_shapification_count() const{
        return cluster_star_shapifications_count_;
    }

    const ExpansionDataLogger& Expander::data_logger() const{
        return data_logger_;
    }


    int Expander::fully_expand_mesh(TetrahedralMesh& mesh,
                                    const TetrahedralMesh& input_mesh,
                                    const std::string& mesh_name,
                                    const std::string& output_file_path,
                                    double& expansion_time_s,
                                    bool silent_mode,
                                    EXPANSION_SCHEME expansion_scheme,
                                    VertexPropertyT<VertexPosition>* exact_vertices_positions){

        Expander expander(mesh,
                          input_mesh,
                          mesh_name,
                          silent_mode,
                          exact_vertices_positions);


        return expander.fully_expand_mesh(output_file_path, expansion_scheme, expansion_time_s);
     }




    int Expander::fully_expand_mesh(const std::string& output_file_path,
                                    EXPANSION_SCHEME expansion_scheme,
                                    double& expansion_time_s){



        if(is_fully_expanded()){
            PRINT_IF_NOT_SILENT(" --> mesh is already expanded -> success I guess"<<std::endl);
            return 0;
        }


        if(!expanding_cluster_){
            std::cout<<" ========================================================"<<std::endl;
            std::cout<<" ========================================================"<<std::endl;
            std::cout<<" FULLY EXPANDING MESH WITH SCHEME "<<expansion_scheme<<"..."<<std::endl;
            std::cout<<" max allocated time: "<<max_allocated_time_s_<<" (s) = "<<((double)max_allocated_time_s_/3600.f)<<" (h)"<<std::endl;
        }

        for(auto v: mesh_.vertices()){
            moved_during_smoothing_prop_[v] = true;
        }

        std::vector<double> iteration_times_s;

        const int initial_vertices_count = mesh_.n_logical_vertices();
        const int initial_boundary_vertices_count = n_boundary_vertices();

        iteration_count_ = 0;
        const int max_iteration_count(std::numeric_limits<int>::max());

        int vertices_increase_ratio(1);
        const int max_vertices_increase_ratio(DEFAULT_MAX_GROWTH_RATIO);


        auto start_time = std::chrono::high_resolution_clock::now();

        int expansion_result = 1;
        bool timeout(false);

        do{
            PRINT_IF_NOT_SILENT(" ========================================================"<<std::endl);
            PRINT_IF_NOT_SILENT(" ========================= RUNNING ITERATION N°"<<iteration_count_<<std::endl);
            PRINT_IF_NOT_SILENT(" expanding cluster: "<<expanding_cluster_<<std::endl);

            if(iteration_count_ == MANUAL_STOP_ITERATION){
                PRINT_IF_NOT_SILENT(" WARNING - manual stop at iteration "<<iteration_count_<<std::endl);
                expansion_result = 1;
                break;
            }

            auto iteration_start_time = std::chrono::high_resolution_clock::now();

            try{
            switch(expansion_scheme){
            case FOUR_STAGE_EXPANSION: {
                expansion_result = four_stage_full_expansion(output_file_path);
                break;
            }

            case FIVE_STAGE_EXPANSION: {
                expansion_result = five_stage_full_expansion(output_file_path);
                break;
            }

            case SIX_STAGE_EXPANSION: {
                expansion_result = six_stage_full_expansion(output_file_path);
                break;
            }

            default: {
                PRINT_IF_NOT_SILENT(" ERROR - unhandled expansion scheme "<<expansion_scheme<<std::endl);
                expansion_result = -1;
                break;
            }
            }


            check_for_timeout();

            if(expansion_result){
                PRINT_IF_NOT_SILENT(" --> expansion scheme failed after "<<iteration_count_<<" iterations, stopping"<<std::endl);
                break;
            }

            }catch(TimeOutException e){
                std::cout<<" --> timeout detected"<<std::endl;
                timeout = true;
            }

#if 0 //SMOOTHING_ENABLED
            if(!expanding_cluster_){
                smooth_interior(DEFAULT_INTERIOR_SMOOTHING_MAX_ITERATIONS,
                                1e-1,
                                DEFAULT_POSITION_BYTE_SIZE_SMOOTHING_THRESHOLD);

                if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_)){
                    PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after smoothing"<<std::endl);
                    expansion_result = -1;
                    break;
                }
            }
#endif

            auto iteration_end_time = std::chrono::high_resolution_clock::now();
            double iteration_duration_s = (float)std::chrono::duration_cast<std::chrono::milliseconds>(iteration_end_time - iteration_start_time).count() /1000.0;

            vertices_increase_ratio = mesh_.n_logical_vertices()/initial_n_vertices_;

            iteration_times_s.push_back(iteration_duration_s);
            PRINT_IF_NOT_SILENT(" ========================= DONE WITH ITERATION N°"<<iteration_count_<<" IN "<<iteration_duration_s<<" SECONDS"<<std::endl);
            PRINT_IF_NOT_SILENT(" == expanded vertices: "<<(expanded_count()+1)<<
                       "/"<<to_expand_count()<<
                       " | diff = "<<(to_expand_count() - expanded_count()-1)<<std::endl);
            PRINT_IF_NOT_SILENT(" == current vertices   growth ratio: "<<((float)mesh_.n_logical_vertices()/(float)initial_n_vertices_)<<std::endl);
            PRINT_IF_NOT_SILENT(" == current vertices deletion ratio: "<<((float)(mesh_.n_vertices() - mesh_.n_logical_vertices())/(float)mesh_.n_logical_vertices())<<std::endl);
            PRINT_IF_NOT_SILENT(" == remaining expansion time (s): "<<remaining_seconds_before_timeout()<<std::endl);
            PRINT_IF_NOT_SILENT(" ========================================================"<<std::endl);
            iteration_count_++;
            data_logger_.add_iteration_timing_point(iteration_duration_s * 1e6);

#if ENABLE_ALL_CHECKS

            if(!expanding_cluster_) {
                auto domain_bad_tets = ExactBadTetFinder::findBadTets(mesh_, domain_vertex_position_prop_);

                if (!domain_bad_tets.first.empty()) {
                    std::cout << " ERROR - domain mesh contains degenerate tets after iteration " << iteration_count_
                              << std::endl;
                    return EXPANSION_ERROR;
                }

                if (!domain_bad_tets.second.empty()) {
                    std::cout << " ERROR - domain mesh contains flipped tets after iteration " << iteration_count_
                              << std::endl;
                    return EXPANSION_ERROR;
                }
            }
#endif

            //temp check
            /*if(!expanding_cluster_ && !all_unexpanded_vertices_are_at_the_center(VertexHandle(-1))){
                std::cout<<" ERROR - some unexpanded vertices are not at the origin"<<std::endl);
                return EXPANSION_FAILURE;
            }*/
            //temp check ends here


        }while(iteration_count_ < max_iteration_count &&
               expanded_count() <= to_expand_count() &&
               !expansion_result &&
               !is_fully_expanded() &&
               //(expanding_cluster_ || vertices_increase_ratio < max_vertices_increase_ratio) &&
               !timeout);


        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_s = ((float)std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count())/1000;

        expansion_time_s = duration_s;

#if ENABLE_FINAL_CHECKS
        auto bad_tets = ExactBadTetFinder::findBadTets(mesh_, vertex_position_prop_);
#else
        BadTetFinder::BadTetList bad_tets;
#endif

        int final_boundary_vertices_count = n_boundary_vertices();

        int max_precision, total_precision;
        double average_precision;
        VertexHandle max_precision_vh;
        compute_precision_stats(mesh_,
                                vertex_position_prop_,
                                max_precision,
                                total_precision,
                                average_precision,
                                max_precision_vh);

        if(!expanding_cluster_){

            std::cout<<" DONE! FULLY EXPANDED MESH IN "<<iteration_count_<<" ITERATIONS, IN "<<duration_s<<" SECONDS"<<std::endl;
            std::cout<<"-------------------------------------------"<<std::endl;
            std::cout<<" MESH IS FULLY EXPANDED: "<<is_fully_expanded()<<std::endl;
            std::cout<<" DEGENERATE TETS COUNT = "<<bad_tets.first.size()<<std::endl;
            std::cout<<"    FLIPPED TETS COUNT = "<<bad_tets.second.size()<<std::endl;
            std::cout<<"-------------------------------------------"<<std::endl;
            std::cout<<" INITIAL VERTICES COUNT = "<<initial_vertices_count<<std::endl;
            std::cout<<"   FINAL VERTICES COUNT = "<<mesh_.n_logical_vertices()<<std::endl;
            std::cout<<"           GROWTH RATIO = "<<((double)mesh_.n_logical_vertices()/(double)initial_vertices_count)<<std::endl;
            std::cout<<"-------------------------------------------"<<std::endl;
            std::cout<<" INITIAL BOUNDARY VERTICES COUNT = "<<initial_boundary_vertices_count<<std::endl;
            std::cout<<"   FINAL BOUNDARY VERTICES COUNT = "<<final_boundary_vertices_count<<std::endl;
            std::cout<<"                    GROWTH RATIO = "<<((double)final_boundary_vertices_count/(double)initial_boundary_vertices_count)<<std::endl;
            std::cout<<"-------------------------------------------"<<std::endl;
            VertexPosition empty_pos(M_PI, M_PI, M_PI);
            std::cout<<"   (values for double-equivalent should be "<<(byte_size(empty_pos))<<", "<<(byte_size(empty_pos))<<" and "<<(byte_size(empty_pos)*mesh_.n_logical_vertices())<<")"<<std::endl;
            std::cout<<"       MAX PRECISION (BYTES) = "<<max_precision<<std::endl;
            std::cout<<"   AVERAGE PRECISION (BYTES) = "<<average_precision<<std::endl;
            std::cout<<" TOTAL POSITION SIZE (BYTES) = "<<total_precision<<std::endl;
            std::cout<<"           TOTAL BYTES SAVED = "<<total_saved_position_bytes_<<std::endl;
            std::cout<<"-------------------------------------------------------------------------"<<std::endl;
            //std::cout<<" SINGLE-SV, TOPO-EXPANDABLE EXPANSES WITH NON-EXPANDABLE SV = "<<single_SV_topo_expandable_but_not_expandable_SV_count_<<std::endl;
            //std::cout<<" SINGLE-SV, TOPO-EXPANDABLE EXPANSES WITH     EXPANDABLE SV = "<<single_SV_topo_expandable_and_expandable_SV_count_<<std::endl;
            std::cout<<" OVERVIEW OF CLUSTER SIZES:"<<std::endl;
            for(int i(2); i < (int)cluster_size_count_.size(); i++){
                std::cout<<" - size "<<i<<" : "<<cluster_size_count_[i]<<std::endl;
            }
            std::cout<<"    SIMPLE STAR-SHAPIFICATION SUCCESSES: "<<star_shapification_successes_count_<<std::endl;
            std::cout<<"   CLUSTER STAR-SHAPIFICATION SUCCESSES: "<<cluster_star_shapifications_count_<<std::endl;
            std::cout<<" TOTAL NUMBER OF POST-SS EDGE COLLAPSES: "<<total_collapsed_edges_count_<<"/"<<split_list_.size()<<std::endl;
            std::cout<<"-------------------------------------------------------------------------"<<std::endl;
            std::cout<<" total           simple expansions time (s): "<< total_simple_expansions_time_s_<<std::endl;
            std::cout<<" total          star-shapification time (s): "<< (total_ss_time_s_ - total_edge_collapsing_time_s_)<<std::endl;
            std::cout<<" total          cluster expansions time (s): "<< total_cluster_exp_time_s_<<std::endl;
            std::cout<<" total cluster star-shapifications time (s): "<< total_css_time_s_<<std::endl;
            std::cout<<" total                   smoothing time (s): "<< total_smoothing_time_s_<<std::endl;
            std::cout<<" total             edge collapsing time (s): "<< total_edge_collapsing_time_s_<<std::endl;
            std::cout<<"-------------------------------------------------------------------------"<<std::endl;
        }
        //temporary fix
       /* if(!bad_tets.first.empty()){
            PRINT_IF_NOT_SILENT(" ---> actually remaining tets. Trying to fix those post-op");
            int fixed_cells_count(0);

            if(post_op_fix_degenerate_cells(bad_tets.first)){
                PRINT_IF_NOT_SILENT(" --> couldn't fix those cells post-op");
            }

            bad_tets = ExactBadTetFinder::findBadTets(codomain_mesh_, vertex_position_prop_);
        }*/

        //data_logger().print_out_expansion_histogram();

        auto domain_bad_tets = ExactBadTetFinder::findBadTets(mesh_, domain_vertex_position_prop_);
        PRINT_IF_NOT_SILENT(" DOMAIN DEGENERATE TETS COUNT = "<<domain_bad_tets.first.size()<<std::endl);
        PRINT_IF_NOT_SILENT(" DOMAIN    FLIPPED TETS COUNT = "<<domain_bad_tets.second.size()<<std::endl);

        int final_result = EXPANSION_FAILURE;

        if(/*!expanding_cluster_ && */timeout){
            std::cout<<" --> TIMEOUT! "<<std::endl;
            final_result = 4;

        /*}else if(!expanding_cluster_ && vertices_increase_ratio >= max_vertices_increase_ratio){
            PRINT_IF_NOT_SILENT(" --> REACHED MAX GROWTH RATIO!"<<std::endl);
            final_result = 5;*/
        }else if(expansion_result == -1){
            PRINT_IF_NOT_SILENT(" --> an error occured during the last expansion"<<std::endl);
            final_result = -1;

        }else if(final_result != is_fully_expanded() ||
                 !bad_tets.first.empty() ||
                 !bad_tets.second.empty()){

            PRINT_IF_NOT_SILENT(" - expanding cluster: "<<expanding_cluster_<<std::endl);
            PRINT_IF_NOT_SILENT(" -   degenerate tets: "<<bad_tets.first.size()<<std::endl);
            PRINT_IF_NOT_SILENT(" -      flipped tets: "<<bad_tets.second.size()<<std::endl);

            if(!bad_tets.first.empty() && expanding_cluster_ && is_fully_expanded()){
                PRINT_IF_NOT_SILENT(" --> remaining degenerate tets after expanding cluster, which is normal"<<std::endl);
                final_result = EXPANSION_SUCCESS;
            }else{

                if (iteration_count_ == max_iteration_count) {
                    final_result = 1;
                }  else if(!bad_tets.first.empty()) {

                    PRINT_IF_NOT_SILENT(" remaining degenerate tets after full expansion: ("<<bad_tets.first.size()<<")"<<std::endl);

                    int max_size(20);
                    if((int)bad_tets.first.size() <= max_size){
                        for(auto tet: bad_tets.first){
                            PRINT_IF_NOT_SILENT(" - "<<tet<<" : "<<mesh_.get_cell_vertices(tet)<<std::endl);

                            for(auto v: mesh_.get_cell_vertices(tet)){
                                PRINT_IF_NOT_SILENT("    -- "<<v<<" at "<<std::setprecision(20)<<vec2vec(vertex_position_prop_[v])<<" is boundary: "<<mesh_.is_boundary(v)<<", neighbors:"<<std::endl);
                                for(auto vv_it = mesh_.vv_iter(v); vv_it.valid(); vv_it++){
                                    PRINT_IF_NOT_SILENT("      --- "<<(*vv_it)<<" at "<<vec2vec(vertex_position_prop_[*vv_it])<<" is boundary: "<<mesh_.is_boundary(*vv_it)<<std::endl);
                                }
                            }
                        }
                    }else{
                        PRINT_IF_NOT_SILENT(" more than "<<max_size<<" degenerate tets, not printing the list"<<std::endl);
                    }

                    final_result = 2;
                }else if (expansion_result == 4){
                    final_result = 4;
                }else if (expansion_result) {
                    final_result = 6;
                }else{
                    std::cout << " unknown error !" << std::endl;
                    final_result = -1;
                }
            }

            if(final_result) {
                std::cout << " --> failed to fully expand mesh, result = " << final_result << std::endl;
            }

        }else{
            PRINT_IF_NOT_SILENT(" --> EXPANSION SUCCESFUL!"<<std::endl);
            final_result = EXPANSION_SUCCESS;
        }

        if(!expanding_cluster_ /*&& !output_file_path.empty()*/){
            data_logger().write_data_to_file(output_file_path,
                                             final_result,
                                             duration_s,
                                             max_precision,
                                             total_precision,
                                             average_precision,
                                             total_saved_position_bytes_);


            auto str_exact_position_prop = mesh_.request_vertex_property<std::string>("str_exact_position");
            auto str_domain_exact_position_prop = mesh_.request_vertex_property<std::string>("str_domain_exact_position");

            for(auto v: mesh_.vertices()){
                std::stringstream sstr;
                sstr << vertex(v)[0]<<";"<<vertex(v)[1]<<";"<<vertex(v)[2];
                str_exact_position_prop[v] = sstr.str();

                std::stringstream domain_sstr;
                domain_sstr << domain_vertex_position_prop_[v][0]<<";"<< domain_vertex_position_prop_[v][1]<<";"<< domain_vertex_position_prop_[v][2];
                str_domain_exact_position_prop[v] = domain_sstr.str();
                //std::cout<<" string domain pos for vertex "<<v<<" : "<<str_domain_exact_position_prop[v] <<std::endl);
                //PRINT_IF_NOT_SILENT(" - string pos for vertex "<<v<<": "<<str_exact_position_prop[v]<<std::endl);

            }
            mesh_.set_shared(str_exact_position_prop);
            mesh_.set_persistent(str_exact_position_prop);
            mesh_.set_shared(str_domain_exact_position_prop);
            mesh_.set_persistent(str_domain_exact_position_prop);

            mesh_.collect_garbage();
        }


        return final_result;
    }





    int Expander::four_stage_full_expansion(const std::string& output_file_path){


        PRINT_IF_NOT_SILENT(" ----> RUNNING 4-STAGE SCHEME"<<std::endl);


        //probably deprecated
        VertexExpanse expanse;


        // FIRST STAGE: simple expansions
        //if(expanding_cluster_){
        PRINT_IF_NOT_SILENT(" WARNING - star-shapifying during first stage"<<std::endl);
        int simple_expansion_result = expand_all_expandable_vertices(expanse, true);
        if(simple_expansion_result == -1){
            PRINT_IF_NOT_SILENT(" error while expanding all expandable vertices without star-shapification"<<std::endl);
            return -1;
        }

        PRINT_IF_NOT_SILENT(" - simple expansion result: "<<simple_expansion_result<<std::endl);


        // SECOND STAGE: simple cluster expansions
        if(simple_expansion_result){

            //simple_cluster_expansion_result = 1;

            PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH SIMPLE EXPANSIONS, EXPANDING CLUSTERS"<<std::endl);

#if STOP_AT_SIMPLE_EXPANSION
            PRINT_IF_NOT_SILENT(" WARNING - HARD-CODED BREAK FOR SIMPLE EXPANSION SAVE"<<std::endl);
            simple_expansion_result = 1;
            return 1;
#endif


            int simple_cluster_expansion_result = find_expandable_cluster_and_expand(false);

            PRINT_IF_NOT_SILENT(" - cluster expansion result: "<<simple_cluster_expansion_result<<std::endl);

            if(simple_cluster_expansion_result == -1){
                PRINT_IF_NOT_SILENT(" error while simply expanding cluster"<<std::endl);
                return -1;
            }

            // THIRD STAGE: star-shapification of single vertices
            if(simple_cluster_expansion_result){

                PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH CLUSTER EXPANSION, STAR-SHAPIFYING SINGLE CONES"<<std::endl);

                int ss_expansion_result = expand_best_expandable_vertex_with_memory(true,
                                                                                DEFAULT_MAX_EXPANSION_VALENCE);

                if(ss_expansion_result == -1){
                    PRINT_IF_NOT_SILENT(" error while star-shapifying single cone"<<std::endl);
                    return -1;
                }

                if(!ss_expansion_result){
                    PRINT_IF_NOT_SILENT(" --> done with star-shapification, moving on to next iteration"<<std::endl);
                }

                //already in expand_best_exp...
                /*if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_)) {
                        std::cout << " ERROR - created flipped tets after single EC star-shapification" << std::endl;
                        break;
                    }*/

                //FOURTH STAGE: cluster star-shapification
                if(ss_expansion_result){

                    PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH SINGLE CONE STAR-SHAPIFICATION, STAR-SHAPIFYING CLUSTERS"<<std::endl);

#if !CLUSTER_SS_ENABLED
                    PRINT_IF_NOT_SILENT(" ---> WARNING: REACHED CLUSTER STAR-SHAPIFICATION, WHICH IS NOT IMPLEMENTED YET"<<std::endl);
                    return 1;
#endif

                    int ss_cluster_expansion_result = find_expandable_cluster_and_expand(true);


                    if(ss_cluster_expansion_result == -1){
                        PRINT_IF_NOT_SILENT(" error while star-shapifying expanding cluster"<<std::endl);
                        return -1;
                    }

                    if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                        std::cout << " ERROR - created flipped tets after star-shapifying cluster" << std::endl;
                        return -1;
                    }


                    if(ss_cluster_expansion_result){
                        PRINT_IF_NOT_SILENT(" COULDN'T STAR-SHAPIFY ANY CLUSTER, OUT OF OPTIONS..."<<std::endl);
                        return -1;
                    }
                }
            }else{

#warning TODO: move this in cluster expansion
                if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                    std::cout << " ERROR - created flipped tets after simple cluster expansion" << std::endl;
                    return -1;
                }
            }
        }

        return 0;
    }



    int Expander::five_stage_full_expansion(const std::string& output_file_path){


        PRINT_IF_NOT_SILENT(" ----> RUNNING 5-STAGE SCHEME"<<std::endl);


        //probably deprecated
        VertexExpanse expanse;


        // FIRST STAGE: simple expansions
        //if(expanding_cluster_){
        //PRINT_IF_NOT_SILENT(" WARNING - star-shapifying during first stage"<<std::endl);
        int simple_expansion_result = expand_all_expandable_vertices(expanse, false);
        if(simple_expansion_result == -1){
            PRINT_IF_NOT_SILENT(" error while expanding all expandable vertices without star-shapification"<<std::endl);
            return -1;
        }

        PRINT_IF_NOT_SILENT(" - simple expansion result: "<<simple_expansion_result<<std::endl);

        // SECOND STAGE: simple cluster expansions with k<=3
        if(simple_expansion_result){

            //simple_cluster_expansion_result = 1;

            PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH SIMPLE EXPANSIONS, EXPANDING CLUSTERS WITH k<=3"<<std::endl);

#if STOP_AT_SIMPLE_EXPANSION
            PRINT_IF_NOT_SILENT(" WARNING - HARD-CODED BREAK FOR SIMPLE EXPANSION SAVE"<<std::endl);
            simple_expansion_result = 1;
            return 1;
#endif


            int simple_cluster_expansion_result = find_expandable_cluster_and_expand(false, 3);

            PRINT_IF_NOT_SILENT(" - cluster expansion result: "<<simple_cluster_expansion_result<<std::endl);

            if(simple_cluster_expansion_result == -1){
                PRINT_IF_NOT_SILENT(" error while simply expanding cluster"<<std::endl);
                return -1;
            }

            // THIRD STAGE: star-shapification of single vertices
            if(!simple_cluster_expansion_result){

                if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                    std::cout << " ERROR - created flipped tets after simple cluster expansion" << std::endl;
                    return -1;
                }
            }else{

                PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH CLUSTER EXPANSION, STAR-SHAPIFYING SINGLE CONES"<<std::endl);

                int ss_expansion_result = expand_best_expandable_vertex_with_memory(true,
                                                                                    DEFAULT_MAX_EXPANSION_VALENCE);

                if(ss_expansion_result == -1){
                    PRINT_IF_NOT_SILENT(" error while star-shapifying single cone"<<std::endl);
                    return -1;
                }

                if(!ss_expansion_result){
                    PRINT_IF_NOT_SILENT(" --> done with star-shapification, moving on to next iteration"<<std::endl);
                }


                //FOURTH STAGE: cluster expansion with k<=5
                if(ss_expansion_result){

                    PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH FIRST STAR-SHAPIFICATION EXPANDING CLUSTER WITH k<=5"<<std::endl);

                    int simple_cluster_expansion_result = find_expandable_cluster_and_expand(false, 5);

                    PRINT_IF_NOT_SILENT(" - cluster expansion result: "<<simple_cluster_expansion_result<<std::endl);

                    if(simple_cluster_expansion_result == -1){
                        PRINT_IF_NOT_SILENT(" error while simply expanding cluster"<<std::endl);
                        return -1;
                    }


                    //FIFTH STAGE: cluster star-shapification
                    if(!simple_cluster_expansion_result){
                        if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                            std::cout << " ERROR - created flipped tets after simple cluster expansion" << std::endl;
                            return -1;
                        }
                    }else{
                        PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH SINGLE CONE STAR-SHAPIFICATION, STAR-SHAPIFYING CLUSTERS"<<std::endl);

#if !CLUSTER_SS_ENABLED
                        PRINT_IF_NOT_SILENT(" ---> WARNING: REACHED CLUSTER STAR-SHAPIFICATION, WHICH IS NOT IMPLEMENTED YET"<<std::endl);
                        return 1;
#endif

                        int ss_cluster_expansion_result = find_expandable_cluster_and_expand(true);

                        if(ss_cluster_expansion_result == -1){
                            PRINT_IF_NOT_SILENT(" error while star-shapifying expanding cluster"<<std::endl);
                            return -1;
                        }

                        if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                            std::cout << " ERROR - created flipped tets after star-shapifying cluster" << std::endl;
                            return -1;
                        }

                        if(ss_cluster_expansion_result){
                            PRINT_IF_NOT_SILENT(" COULDN'T STAR-SHAPIFY ANY CLUSTER, OUT OF OPTIONS..."<<std::endl);
                            return -1;
                        }
                    }
                }
            }
        }

        return 0;
    }




    int Expander::six_stage_full_expansion(const std::string& output_file_path){


        PRINT_IF_NOT_SILENT(" ----> RUNNING 6-STAGE SCHEME"<<std::endl);

        LightWeightStopWatch stop_watch;

        //probably deprecated
        VertexExpanse expanse;

        // FIRST STAGE: simple expansions
        //PRINT_IF_NOT_SILENT(" WARNING - star-shapifying during first stage"<<std::endl);
        int simple_expansion_result = expand_all_expandable_vertices(expanse, false);
        if(simple_expansion_result == -1){
            PRINT_IF_NOT_SILENT(" error while expanding all expandable vertices without star-shapification"<<std::endl);
            return -1;
        }
        total_simple_expansions_time_s_ += stop_watch.lap_duration();

        PRINT_IF_NOT_SILENT(" - simple expansion result: "<<simple_expansion_result<<std::endl);

        // SECOND STAGE: simple cluster expansions with k<=3
        if(simple_expansion_result){

            //simple_cluster_expansion_result = 1;

            PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH SIMPLE EXPANSIONS, EXPANDING CLUSTERS WITH k<=3"<<std::endl);

#if STOP_AT_SIMPLE_EXPANSION
            PRINT_IF_NOT_SILENT(" WARNING - HARD-CODED BREAK FOR SIMPLE EXPANSION SAVE"<<std::endl);
            simple_expansion_result = 1;
            return 1;
#endif


            bool DEBUG_direct_cluster_ss_enabled(false);


            int simple_cluster_expansion_result = find_expandable_cluster_and_expand(DEBUG_direct_cluster_ss_enabled, 3);

            PRINT_IF_NOT_SILENT(" - cluster expansion result: "<<simple_cluster_expansion_result<<std::endl);

            if(simple_cluster_expansion_result == -1){
                PRINT_IF_NOT_SILENT(" error while simply expanding cluster"<<std::endl);
                return -1;
            }
            total_cluster_exp_time_s_ += stop_watch.lap_duration();

            // THIRD STAGE: star-shapification of single vertices
            if(!simple_cluster_expansion_result){

#if ENABLE_ALL_CHECKS
                if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                    std::cout << " ERROR - created flipped tets after simple cluster expansion" << std::endl;
                    return -1;
                }
#endif

            }else{

                const int max_valence(30);
                PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH CLUSTER EXPANSION, STAR-SHAPIFYING SINGLE CONES WITH MAX VALENCE "<<max_valence<<std::endl);

                int ss_expansion_result = expand_best_expandable_vertex_with_memory(true,
                                                                                    DEFAULT_MAX_EXPANSION_VALENCE,
                                                                                    max_valence);

                if(ss_expansion_result == -1){
                    PRINT_IF_NOT_SILENT(" error while star-shapifying single cone"<<std::endl);
                    return -1;
                }

                if(!ss_expansion_result){
                    PRINT_IF_NOT_SILENT(" --> done with star-shapification, moving on to next iteration"<<std::endl);
                }

                total_ss_time_s_ += stop_watch.lap_duration();


                //FOURTH STAGE: cluster expansion with k<=5
                if(ss_expansion_result){

                    PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH FIRST STAR-SHAPIFICATION EXPANDING CLUSTER WITH k<=5"<<std::endl);

                    int simple_cluster_expansion_result = find_expandable_cluster_and_expand(false, 5);

                    PRINT_IF_NOT_SILENT(" - cluster expansion result: "<<simple_cluster_expansion_result<<std::endl);

                    if(simple_cluster_expansion_result == -1){
                        PRINT_IF_NOT_SILENT(" error while simply expanding cluster"<<std::endl);
                        return -1;
                    }

                    total_cluster_exp_time_s_ += stop_watch.lap_duration();

                    //FIFTH STAGE: star-shapification with any valence
                    if(!simple_cluster_expansion_result){
#if ENABLE_ALL_CHECKS
                        if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                            std::cout << " ERROR - created flipped tets after simple cluster expansion" << std::endl;
                            return -1;
                        }
#endif
                    }else{

                        //not sure about this
                        /*if(collapse_new_edges_opposite_to_unexpanded_vertices()){
                            PRINT_IF_NOT_SILENT(" error while collapsing new edges"<<std::endl);
                            return EXPANSION_ERROR;
                        }*/

                        PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH CLUSTER EXPANSION, STAR-SHAPIFYING SINGLE CONES WITHOUT MAX VALENCE LIMIT "<<std::endl);

                        int ss_expansion_result = expand_best_expandable_vertex_with_memory(true,
                                                                                            DEFAULT_MAX_EXPANSION_VALENCE);

                        if(ss_expansion_result == -1){
                            PRINT_IF_NOT_SILENT(" error while star-shapifying single cone"<<std::endl);
                            return -1;
                        }

                        if(!ss_expansion_result){
                            PRINT_IF_NOT_SILENT(" --> done with star-shapification, moving on to next iteration"<<std::endl);
                        }

                        total_ss_time_s_ += stop_watch.lap_duration();


                        //SIXTH STAGE: cluster star-shapification
                        if(ss_expansion_result){

                            PRINT_IF_NOT_SILENT(" ============================> COULDN'T FULLY EXPAND MESH WITH SINGLE CONE STAR-SHAPIFICATION, STAR-SHAPIFYING CLUSTERS"<<std::endl);

#if !CLUSTER_SS_ENABLED
                            PRINT_IF_NOT_SILENT(" ---> WARNING: REACHED CLUSTER STAR-SHAPIFICATION, WHICH IS NOT IMPLEMENTED YET"<<std::endl);
                            return 1;
#endif
                            int max_n_vertices(50);
                            int ss_cluster_expansion_result = find_expandable_cluster_and_expand(true, left_to_expand_count()/2, max_n_vertices);


                            if(ss_cluster_expansion_result == -1){
                                PRINT_IF_NOT_SILENT(" error while star-shapifying expanding cluster"<<std::endl);
                                return -1;
                            }

#if ENABLE_ALL_CHECKS
                            if (ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)) {
                                std::cout << " ERROR - created flipped tets after star-shapifying cluster" << std::endl;
                                return -1;
                            }
#endif

                            total_css_time_s_ += stop_watch.lap_duration();

                            if(ss_cluster_expansion_result){
                                PRINT_IF_NOT_SILENT(" COULDN'T STAR-SHAPIFY ANY CLUSTER, OUT OF OPTIONS..."<<std::endl);
                                return -1;
                            }
                        }
                    }
                }
            }
        }

        return 0;
    }




    bool Expander::at_least_one_expanded_vertex(TetrahedralMesh& mesh){

        auto expanded_prop = mesh.request_vertex_property<bool>("expanded_vertex");
        for(auto v: mesh.vertices()){
            if(expanded_prop[v]){
                return true;
            }
        }
        return false;
    }



    int Expander::check_main_cluster_topology(bool& single_connected_component,
                                               int& euler_characteristic){


        PRINT_IF_NOT_SILENT(" ----------------------------"<<std::endl);
        PRINT_IF_NOT_SILENT(" checking main cluster topology..."<<std::endl);
        //first, compute main cluster euler charac.
        int unexpanded_vertices_count(0);
        int unexpanded_edges_count(0);
        int unexpanded_faces_count(0);
        int unexpanded_cells_count(0);

        for(auto v: mesh_.vertices()){
            unexpanded_vertices_count += !expanded_prop_[v];
        }
        PRINT_IF_NOT_SILENT(" - found "<<unexpanded_vertices_count<<" unexpanded vertices"<<std::endl);

        for(auto e: mesh_.edges()){
            auto from_vertex = mesh_.from_vertex_handle(mesh_.halfedge_handle(e,0));
            auto to_vertex   = mesh_.to_vertex_handle(mesh_.halfedge_handle(e,0));

            unexpanded_edges_count += (!expanded_prop_[from_vertex] && !expanded_prop_[to_vertex]);
        }
        PRINT_IF_NOT_SILENT(" - found "<<unexpanded_edges_count<<" unexpanded edges"<<std::endl);


        for(auto f: mesh_.faces()){
            auto f_vertices = mesh_.get_halfface_vertices(mesh_.halfface_handle(f, 0));
            if(f_vertices.size() != 3){
                PRINT_IF_NOT_SILENT("ERROR - face "<<f<<" contains "<<f_vertices.size()<<std::endl);
                return -1;
            }

            unexpanded_faces_count += (!expanded_prop_[f_vertices[0]] &&
                                       !expanded_prop_[f_vertices[1]] &&
                                       !expanded_prop_[f_vertices[2]]);
        }
        PRINT_IF_NOT_SILENT(" - found "<<unexpanded_faces_count<<" unexpanded faces"<<std::endl);


        for(auto c: mesh_.cells()){
            auto c_vertices = mesh_.get_cell_vertices(c);
            if(c_vertices.size() != 4){
                PRINT_IF_NOT_SILENT("ERROR - cell "<<c<<" contains "<<c_vertices.size()<<std::endl);
                return -1;
            }

            unexpanded_cells_count += (!expanded_prop_[c_vertices[0]] &&
                                       !expanded_prop_[c_vertices[1]] &&
                                       !expanded_prop_[c_vertices[2]] &&
                                       !expanded_prop_[c_vertices[3]]);
        }
        PRINT_IF_NOT_SILENT(" - found "<<unexpanded_cells_count<<" unexpanded cells"<<std::endl);

        euler_characteristic = unexpanded_vertices_count -
                               unexpanded_edges_count +
                               unexpanded_faces_count -
                               unexpanded_cells_count;

        PRINT_IF_NOT_SILENT(" --> euler charac. = "<<euler_characteristic<<std::endl);


        auto visited_prop = mesh_.request_vertex_property<bool>();

        std::queue<VertexHandle> to_visit;
        for(auto v: mesh_.vertices()){
            if(!expanded_prop_[v]){
                to_visit.push(v);
                break;
            }
        }

        if(to_visit.empty()){
            PRINT_IF_NOT_SILENT(" ==> main cluster is empty"<<std::endl);
            return -1;
        }

        int visited_count(0);
        while(!to_visit.empty()){
            auto current_vertex = to_visit.front();
            to_visit.pop();
            if(!visited_prop[current_vertex]) {

                //std::cout << " -- current vertex: " << current_vertex << std::endl;

                for (auto out_he: mesh_.outgoing_halfedges(current_vertex)) {
                    auto neighbor = mesh_.to_vertex_handle(out_he);
                    if (!visited_prop[neighbor] && !expanded_prop_[neighbor]) {
                        //std::cout << " -- found unvisited, unexpanded neighbor " << neighbor << std::endl;
                        to_visit.push(neighbor);
                    }
                }
                visited_prop[current_vertex] = true;
                visited_count++;
            }
        }


        if(visited_count != (int)left_to_expand_count()){
            PRINT_IF_NOT_SILENT(" ==> "<<left_to_expand_count()<<" vertices left to expand but only "<<visited_count<<" were visited -> main cluster is not a single-CC"<<std::endl);

            PRINT_IF_NOT_SILENT(" - unexpanded, unvisited vertices: ");
            for(auto v: mesh_.vertices()){
                if(!expanded_prop_[v] && !visited_prop[v]){
                    PRINT_IF_NOT_SILENT(" "<<v);
                }
            }
            PRINT_IF_NOT_SILENT(std::endl);
            single_connected_component = false;
        }else{
            single_connected_component = true;
        }


        PRINT_IF_NOT_SILENT(" - visited "<<visited_count<<" unexpanded vertices"<<std::endl);
        PRINT_IF_NOT_SILENT(" - and there are "<<left_to_expand_count()<<" vertices left to expand"<<std::endl);
        PRINT_IF_NOT_SILENT(" --> single CC: "<<single_connected_component<<std::endl);

        PRINT_IF_NOT_SILENT(" ----------------------------"<<std::endl);
        return euler_characteristic != 1 || !single_connected_component;
    }



    void Expander::reset_EC_memory(){
        for(auto v: mesh_.vertices()){
            expanse_at_last_iteration_prop_[v] = VertexExpanse();
            unexpanded_neighbors_at_expanse_computation_[v].clear();
        }
    }


    EXPANSION_RESULT_STATUS Expander::expand_all_expandable_vertices(VertexExpanse& expanse,
                                                                     bool enable_star_shapification){

        EXPANSION_RESULT_STATUS exp_result(EXPANSION_FAILURE);
        const int mod(20);
        const int max_iterations(10);
        int iteration_count(0);
        bool expanded_something(false);

        do{

            PRINT_IF_NOT_SILENT(" ========================================================================================================================"<<std::endl);
            PRINT_IF_NOT_SILENT(" RUNNING SIMPLE EXPANSIONS LOOP N°"<<(iteration_count+1)<<std::endl);
            int initial_exp_count = expanded_count();
#if SMOOTHING_ENABLED
            if(!expanding_cluster_ &&
                    !is_fully_expanded()){

                smooth_neighborhood_of_unexpanded_vertices(DEFAULT_INTERIOR_SMOOTHING_RING_K,
                                                           DEFAULT_INTERIOR_SMOOTHING_MAX_ITERATIONS,
                                                           DEFAULT_INTERIOR_SMOOTHING_EPSILON,
                                                           DEFAULT_POSITION_BYTE_SIZE_SMOOTHING_THRESHOLD);

                /*smooth_full_interior(DEFAULT_INTERIOR_SMOOTHING_MAX_ITERATIONS,
                                     DEFAULT_INTERIOR_SMOOTHING_EPSILON,
                                     DEFAULT_POSITION_BYTE_SIZE_SMOOTHING_THRESHOLD);*/

#if ENABLE_ALL_CHECKS

                if(ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)){
                    PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after smoothing"<<std::endl);
                    return EXPANSION_ERROR;
                }
#endif
            }
#endif

            LightWeightStopWatch sw;
            expanded_something = false;
            do{
                if(!((expanded_count() + 1) % mod)){
                    PRINT_IF_NOT_SILENT(" ======================================="<<std::endl);
                    PRINT_IF_NOT_SILENT(" ======= EXPANDING VERTEX N°"<<(expanded_count()+1)<<
                               "/"<<to_expand_count()<<
                               " | diff = "<<(to_expand_count() - expanded_count()-1)<<std::endl);
                }
                exp_result = expand_best_expandable_vertex_with_memory(enable_star_shapification,
                                                                       DEFAULT_MAX_EXPANSION_VALENCE);

                if(exp_result == -1){
                    PRINT_IF_NOT_SILENT(" error while expanding best expandable vertex"<<std::endl);
                    return EXPANSION_ERROR;
                }

                expanded_something |= (exp_result == EXPANSION_SUCCESS);

            }while(exp_result == EXPANSION_SUCCESS &&
                   !is_fully_expanded());


            int iteration_exp_count = expanded_count() - initial_exp_count;
            PRINT_IF_NOT_SILENT(" ====> EXPANDED "<<iteration_exp_count<<" VERTICES AT ITERATION "<<(iteration_count + 1)<<std::endl);
            PRINT_IF_NOT_SILENT(" ========================================================================================================================"<<std::endl);

            check_for_timeout();
            //PRINT_IF_NOT_SILENT(" --> expanded something at iteration "<<(iteration_count+1)<<": "<<expanded_something<<std::endl);

            iteration_count++;
        }while(SMOOTHING_ENABLED &&
               expanded_something &&
               iteration_count < max_iterations &&
               !is_fully_expanded());

            PRINT_IF_NOT_SILENT(" DONE TRYING TO EXPAND ALL EXPANDABLE VERTICES."<<std::endl);
        //NOTE: +1 for the -1 (error) result space
        std::vector<int> result_count(EXPANSION_CHECK_RESULT_COUNT + 1);
        for(int i(0); i < (int)mesh_.n_vertices(); i++) {
            auto v = VertexHandle(i);
            if (!mesh_.is_boundary(v) &&
                !expanded_prop_[v]) {

                //auto full_cone_expandability = full_cone_at_last_iteration_prop_[v].second;
                auto full_cone_expandability = expanse_at_last_iteration_prop_[v].exp_result;
                //+1 so the errors go first and everything is shifted
                result_count[full_cone_expandability + 1]++;
            }
        }

        int total(0);
        PRINT_IF_NOT_SILENT(" OVERVIEW OF STATUS OF REMAINING ONES: "<<std::endl);
        PRINT_IF_NOT_SILENT(" ---------------------------------- "<<std::endl);
        PRINT_IF_NOT_SILENT(" ---------------------------------- "<<std::endl);
        for(int i(0); i<(int)result_count.size(); i++) {
            PRINT_IF_NOT_SILENT(" -- "<<(i-1)<<" : "<<result_count[i]<<std::endl);
            total += result_count[i];
        }
        PRINT_IF_NOT_SILENT(" -------------"<<std::endl);
        PRINT_IF_NOT_SILENT(" total: "<<total<<std::endl);
        PRINT_IF_NOT_SILENT(" ---------------------------------- "<<std::endl);
        PRINT_IF_NOT_SILENT(" ---------------------------------- "<<std::endl);

        return exp_result;
    }




    int Expander::merge_expansion_cones(const std::vector<VertexHandle>& unexpanded_vertices,
                                        const std::vector<int>& cluster_indices,
                                        ExpansionCone& cluster_ec){

        /*PRINT_IF_NOT_SILENT(" ------ merging cluster with cones "<<std::endl);
        int i(0);
        for(auto c: cluster_indices){
            PRINT_IF_NOT_SILENT(" - "<<i<<" : "<<c<<" "<<unexpanded_vertices[c]<<std::endl);
            i++;
        }*/

        cluster_ec.clear();

        //first, check that there's a single spanning tree for all tip vertices
        //PRINT_IF_NOT_SILENT(" tip vertices: ";
        auto tip_vertices_prop = mesh_.request_vertex_property<bool>();
        for(auto c: cluster_indices){
            tip_vertices_prop[unexpanded_vertices[c]] = true;
            //PRINT_IF_NOT_SILENT(unexpanded_vertices[c]<<" ";
        }
        //PRINT_IF_NOT_SILENT(std::endl;

#warning TODO: make this optional
        auto visited_tips_prop = mesh_.request_vertex_property<bool>();

        std::queue<VertexHandle> tips_to_visit;
        tips_to_visit.push(unexpanded_vertices[cluster_indices[0]]);

        std::set<VertexHandle> spanned_tips;
        spanned_tips.insert(tips_to_visit.front());

        while(!tips_to_visit.empty() &&
              spanned_tips.size() != cluster_indices.size()){

            //PRINT_IF_NOT_SILENT(" -----"<<std::endl);
            auto current_tip = tips_to_visit.front();
            tips_to_visit.pop();
            //PRINT_IF_NOT_SILENT(" - current tip "<<current_tip<<std::endl);

            for(auto out_he: mesh_.outgoing_halfedges(current_tip)){
                auto neighbor = mesh_.to_vertex_handle(out_he);
                //PRINT_IF_NOT_SILENT(" -- trying neighbor "<<neighbor<<std::endl);
                if(!visited_tips_prop[neighbor] &&
                    tip_vertices_prop[neighbor]){
                    tips_to_visit.push(neighbor);
                    spanned_tips.insert(neighbor);
                    //PRINT_IF_NOT_SILENT(" ---> unvisited tip "<<neighbor<<", updated queue size: "<<tips_to_visit.size()<<std::endl);

                    if(spanned_tips.size() == cluster_indices.size()){
                        //PRINT_IF_NOT_SILENT(" ----> spanned all tip vertices, stopping"<<std::endl);
                        break;
                    }
                }
            }
            visited_tips_prop[current_tip] = true;
        }

        if(spanned_tips.size() != cluster_indices.size()){
            //PRINT_IF_NOT_SILENT(" --> found unvisited tip -> cannot merge"<<std::endl);
            return 1;
        }

        //PRINT_IF_NOT_SILENT(" --> found single spanning tree for all cluster tip vertices"<<std::endl);

        cluster_ec.clear();

        for(auto c: cluster_indices){
            //expanded_prop_[unexpanded_vertices[c]] = true;
            set_expanded_prop(unexpanded_vertices[c], true);
            //PRINT_IF_NOT_SILENT(" - set "<<unexpanded_vertices[c]<<" as artificially expanded"<<std::endl);
        }


        for(auto c: cluster_indices){
            gather_full_expansion_cone(unexpanded_vertices[c], cluster_ec, true);

            /*PRINT_IF_NOT_SILENT(" - cluster after adding cone "<<c
                       <<" with tip vertex "<<unexpanded_vertices[c]
                       <<" : "<<cluster<<std::endl);
            cluster.print_details();*/


            //to avoid duplicate entities
            //expanded_prop_[unexpanded_vertices[c]] = false;
            set_expanded_prop(unexpanded_vertices[c], false);
        }

        for(auto c: cluster_indices){
            //expanded_prop_[unexpanded_vertices[c]] = false;
            set_expanded_prop(unexpanded_vertices[c], false);
        }


        return 0;
    }




    int Expander::cluster_interface_expansion(const ExpansionCone& cluster_ec){


        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" - PERFORMING CLUSTER INTERFACE EXPANSION..."<<std::endl);
        //PRINT_IF_NOT_SILENT(" cluster EC: "; cluster_ec.print_details();
        //std::cout<<" - cluster interface expansion for cluster "; cluster_ec.print_details();

        //gather mesh halfedges to split
        std::vector<HalfEdgeHandle> cluster_connecting_hes;
        for(auto cone_tip: cluster_ec.cone_tip_vertices()){
            auto mesh_tip = cluster_ec.cone_to_mesh_handle(cone_tip);
            for(auto out_he: mesh_.outgoing_halfedges(mesh_tip)){
                auto neighbor(mesh_.to_vertex_handle(out_he));
                auto cone_neighbor = cluster_ec.mesh_to_cone_handle(neighbor);
                if(!expanded_prop_[neighbor] &&
                        !(cone_neighbor.idx() != -1 && cluster_ec.is_cone_tip(cone_neighbor))){
                    PRINT_IF_NOT_SILENT(" - base vertex "<<neighbor<<" is part of the primary cluster, adding edge "<<mesh_.halfedge(out_he)<<" to list"<<std::endl);

                    cluster_connecting_hes.push_back(out_he);
                }
            }
        }

        //PRINT_IF_NOT_SILENT(" temporary failure for checks"<<std::endl);
        //return -1;

        PRINT_IF_NOT_SILENT(" - found "<<cluster_connecting_hes.size()<<" cluster-connecting edges"<<std::endl);
        if(cluster_connecting_hes.empty()){
            return 0;
        }

        std::vector<Split> local_split_list;

        //split them and temporarily move mid-vertices to the midpoint of their edge
        for(auto he: cluster_connecting_hes){
            auto mesh_from_v(mesh_.from_vertex_handle(he));
            auto mesh_to_v(mesh_.to_vertex_handle(he));
            auto mid_vertex = split_edge(he);
            auto mid_pos = (vertex(mesh_from_v) + vertex(mesh_to_v))/2;
            //to update the exact vertex position prop
            this->set_vertex(mid_vertex, mid_pos);

            local_split_list.push_back({mesh_from_v,
                                        mesh_to_v,
                                        mid_vertex,
                                        mid_pos});

            //expanded_prop_[mid_vertex] = false;
            set_expanded_prop(mid_vertex, false);

            //PRINT_IF_NOT_SILENT(" -- added split "<<local_split_list.back()<<std::endl);
        }



        int local_expanded_count(0);
        int i(0);
        while(i < (int)local_split_list.size()){
            //PRINT_IF_NOT_SILENT(" ----------------- expanding vertex n°"<<i<<std::endl);
            for(auto& split: local_split_list){
                auto to_expand = split.cone_new_vertex;
                if(!expanded_prop_[to_expand]){
                    /*if(local_expanded_count >= ((int)local_split_list.size() -1)){
                        PRINT_IF_NOT_SILENT(" -----> only one remaining unexpanded mid-vertex "<<to_expand<<", marking as expanded and adding to split list"<<std::endl);
                        //expanded_prop_[to_expand] = true;
                        set_expanded_prop(to_expand, true);
                        //split_list_.push_back(split);
                        local_expanded_count++;
                        break;
                    }*/
                    //PRINT_IF_NOT_SILENT(" - trying to expand mid-vertex "<<to_expand<<std::endl);

                    ExpansionCone cone;
                    int one_ring_set_res = ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(mesh_,
                                                                                                       vertex_position_prop_,
                                                                                                       to_expand,
                                                                                                       cone);

                    if(one_ring_set_res){
                        PRINT_IF_NOT_SILENT(" error while setting-up 1-ring neighborhood EC for mid-vertex"<<std::endl);
                        return -1;
                    }

                    //PRINT_IF_NOT_SILENT(" - 1-ring neighborhood of "<<to_expand<<": "; cone.print_details();

                    VertexPosition pos;
                    auto exp_result = cone.is_geo_expandable(pos);
                    if(exp_result == -1){
                        PRINT_IF_NOT_SILENT(" --> error while checking expandability of 1-ring neighborhood of mid-vertex "<<to_expand<<std::endl);
                        PRINT_IF_NOT_SILENT(" cone details: "); if(!silent_mode_){cone.print_details();}
                        return -1;
                    }else if(exp_result){
                        //PRINT_IF_NOT_SILENT(" --> couldn't expand 1-ring neighborhood of mid-vertex "<<to_expand<<", result: "<<exp_result<<std::endl);
                    }else{
                        //PRINT_IF_NOT_SILENT(" --> mid-vertex "<<to_expand<<" is expandable!"<<std::endl);

                        set_vertex(to_expand, pos);
                        //expanded_prop_[to_expand] = true;
                        set_expanded_prop(to_expand, true);
                        split.new_vertex_position = pos;
                        local_expanded_count++;
                        //std::cout<<" --> updated split "<<split<<" with expansion position"<<std::endl);
                    }
                }
            }
            i++;
        }

        for(const auto& split: local_split_list){
            split_list_.push_back(split);
            PRINT_IF_NOT_SILENT(" --> added split "<<split_list_.back()<<" to expander split list"<<std::endl);
            auto mesh_from_vertex = split.cone_from_vertex;
            auto mesh_to_vertex   = split.cone_to_vertex;
            auto mesh_new_vertex  = split.cone_new_vertex;

            mark_halfedges_as_new_and_collapsible(split.cone_from_vertex,
                                                  split.cone_new_vertex,
                                                  split.cone_to_vertex);
        }

        if(local_expanded_count != (int)local_split_list.size()){
            PRINT_IF_NOT_SILENT(" ERROR - could only expand "<<local_expanded_count<<" mid-vertices among "<<local_split_list.size()<<std::endl);
            return -1;
        }

        PRINT_IF_NOT_SILENT(" ...done! Expanded "<<local_expanded_count<<"/"<<local_split_list.size()<<" mid-vertices"<<std::endl);
        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);

        return 0;
    }




    int Expander::expand_cluster(const ExpansionCone& cluster,
                                 const VertexPosition& pos){


        PRINT_IF_NOT_SILENT(" ============================================================"<<std::endl);
        PRINT_IF_NOT_SILENT(" =========== EXPANDING CLUSTER "<<cluster<<std::endl);
        PRINT_IF_NOT_SILENT(" ============================================================"<<std::endl);




        //cluster + all the unexpanded neighboring stuff in the primary cluster

        ExpansionCone cluster_submesh;


        //first, set the cluster's tip vertices to the position

        auto mesh_tip_vertices_prop = mesh_.request_vertex_property<bool>();
        //add the tips
        //PRINT_IF_NOT_SILENT(" - added tips: "<<std::endl);
        for(auto cone_tip_vertex: cluster.cone_tip_vertices()) {
            auto mesh_tip_vertex = cluster.cone_to_mesh_handle(cone_tip_vertex);
            auto extended_tip_vertex = cluster_submesh.add_vertex(mesh_tip_vertex,
                                                                   vertex(mesh_tip_vertex),
                                                                   true);

            this->set_vertex(mesh_tip_vertex, pos);
            cluster_submesh.set_vertex(extended_tip_vertex, pos);

            mesh_tip_vertices_prop[mesh_tip_vertex] = true;
            //PRINT_IF_NOT_SILENT("    - "<<cone_tip_vertex<<"/"<<mesh_tip_vertex<<"/"<<extended_tip_vertex<<std::endl);
        }


        //do the cluster interface expansion
        int interface_exp_result = cluster_interface_expansion(cluster);
        if(interface_exp_result){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't expand cluster interface"<<std::endl);
            return -1;
        }


#if ENABLE_ALL_CHECKS
        if(ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)){
            PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after performing cluster interface expansion. Flipped tets:"<<std::endl);
            auto bad_tets = ExactBadTetFinder::findBadTets(mesh_, vertex_position_prop_);
            for(auto flipped_tet: bad_tets.second){
                auto volume = OVMtetToCGALtet(mesh_, vertex_position_prop_, flipped_tet).volume();
                PRINT_IF_NOT_SILENT(" - "<<flipped_tet<<" : "<<mesh_.get_cell_vertices(flipped_tet)<<", volume = "<<volume<<std::endl);
            }
            return -1;
        }
#endif


        //and then all cells to the cluster submesh
        //PRINT_IF_NOT_SILENT(" - adding cells... "<<std::endl);
        auto added_mesh_cell_prop = mesh_.request_cell_property<bool>();

        for(auto cone_tip_vertex: cluster.cone_tip_vertices()) {
            auto mesh_tip_vertex = cluster.cone_to_mesh_handle(cone_tip_vertex);

            //PRINT_IF_NOT_SILENT(" ----- adding cells around tip vertex "<<mesh_tip_vertex<<"/"<<cone_tip_vertex<<": "<<std::endl);

            for(auto vc_it = mesh_.vc_iter(mesh_tip_vertex); vc_it.valid(); vc_it++){

                if(!added_mesh_cell_prop[*vc_it]) {
                    auto extended_cell = add_cell_to_cone(*vc_it, cluster_submesh);
                    if(!extended_cell.is_valid()){
                        PRINT_IF_NOT_SILENT(" ERROR - couldn't add cell "<<*vc_it<<" : "<<mesh_.get_cell_vertices(*vc_it)<<std::endl);
                        return -1;
                    }
                    added_mesh_cell_prop[*vc_it] = true;
                    //PRINT_IF_NOT_SILENT("   - "<<*vc_it<<"/"<<extended_cell<<" : "<<extended_cluster.get_cell_vertices(extended_cell)<<std::endl);
                }
            }
        }

        cluster_submesh.collect_garbage();

        //then set-up the secondary Expander


        /*PRINT_IF_NOT_SILENT(" - cluster submesh: "<<std::endl);
        cluster_submesh.print_details();*/

        /*PRINT_IF_NOT_SILENT(" cluster interior vertices: "<<std::endl);
        for(auto v: cluster_submesh.vertices()) {
            if (!cluster_submesh.is_boundary(v)) {
                std::cout << "  - " << v <<
                             " at " << vec2vec(cluster_submesh.vertex(v)) <<
                             " is boundary: " << cluster_submesh.is_boundary(v) << std::endl;
            }
        }*/

        /*PRINT_IF_NOT_SILENT(" cluster boundary vertices: "<<std::endl);
        for(auto v: cluster_submesh.vertices()) {
            if (cluster_submesh.is_boundary(v)) {
                std::cout << "  - " << v <<
                             " at " << vec2vec(cluster_submesh.vertex(v)) <<
                             " is boundary: " << cluster_submesh.is_boundary(v) << std::endl;
            }
        }*/

        /*PRINT_IF_NOT_SILENT(" cluster boundary halffaces: "<<std::endl);
        for(auto hf: extended_cluster.halffaces()){
            if(extended_cluster.is_boundary(hf)){
                PRINT_IF_NOT_SILENT(" - "<<hf<<" : "<<extended_cluster.get_halfface_vertices(hf)<<std::endl);
            }
        }*/




        TetrahedralMesh cluster_mesh = cluster_submesh;
        auto vertex_position_prop = cluster_submesh.vertex_position_prop();

        PRINT_IF_NOT_SILENT(" remaining seconds before timeout: "<<remaining_seconds_before_timeout()<<std::endl);
        Expander cluster_expander(cluster_mesh,
                                  cluster_mesh,
                                  "cluster_mesh",
                                  silent_mode_,
                                  &vertex_position_prop,
                                  true,
                                  false,
                                  remaining_seconds_before_timeout() + 1);

        /*if(cluster_submesh.n_vertices() < 50){
            PRINT_IF_NOT_SILENT(" -- cluster submesh: "; cluster_submesh.print_details();
        }else{
            PRINT_IF_NOT_SILENT(" -- cluster submesh contains more than 50 vertices, change this line if you really want to print it out"<<std::endl);
        }*/


        //expand the cluster (recursion magic happening here)
        double cluster_expansion_time_s(0);
        auto cluster_expansion_result = cluster_expander.fully_expand_mesh("temp", DEFAULT_EXPANSION_SCHEME, cluster_expansion_time_s);

        if(cluster_expansion_result){
            PRINT_IF_NOT_SILENT(" ERROR - coudln't expand cluster"<<std::endl);

            check_for_timeout();

            return -1;
        }

        PRINT_IF_NOT_SILENT(" --> done expanding cluster, applying expansion..."<<std::endl);

        //update the cluster index
        latest_cluster_index_++;

        //then update the positions and status
        for(auto extended_cone_tip_vertex: cluster_submesh.cone_tip_vertices()){
            auto mesh_tip_vertex = cluster_submesh.cone_to_mesh_handle(extended_cone_tip_vertex);
            move_vertex_and_mark_as_expanded(mesh_tip_vertex, cluster_expander.vertex(extended_cone_tip_vertex));
            PRINT_IF_NOT_SILENT("  --> moved cluster tip vertex "<<extended_cone_tip_vertex<<
                       " to "<<vec2vec(vertex(mesh_tip_vertex))<<std::endl);


            if(cluster_index_prop_[mesh_tip_vertex] > 0){
                PRINT_IF_NOT_SILENT(" ERROR - cluster tip vertex "<<mesh_tip_vertex<<" is already part of cluster "<<cluster_index_prop_[mesh_tip_vertex]<<std::endl);
                return -1;
            }

            cluster_index_prop_[mesh_tip_vertex] = latest_cluster_index_;
        }


        if(!cluster_expander.get_split_list().empty()){
            PRINT_IF_NOT_SILENT(" - CLUSTER REQUIRED SPLITS, REPRODUCING THEM IN BASE MESH"<<std::endl);
            auto submesh_to_mesh_prop = cluster_mesh.request_vertex_property<VertexHandle>();
            for(auto v: cluster_submesh.vertices()){
                submesh_to_mesh_prop[v] = cluster_submesh.cone_to_mesh_handle(v);
            }

            /*PRINT_IF_NOT_SILENT(" initial new submesh->mesh prop for cluster mesh : "<<std::endl);
            for(auto v: cluster_mesh.vertices()){
                //submesh_to_mesh_prop[v] = extended_cluster.cone_to_mesh_handle(v);
                PRINT_IF_NOT_SILENT(" - "<<v<<" -> "<<submesh_to_mesh_prop[v]<<std::endl);
            }*/

            int split_list_start_index = split_list_.size();

            auto cluster_application_result = apply_split_list(cluster_expander.get_split_list(),
                                                               submesh_to_mesh_prop);

            if(cluster_application_result){
                PRINT_IF_NOT_SILENT(" ERROR - couldn't apply cluster split list to base mesh"<<std::endl);
                return -1;
            }


            reduce_mid_vertices_precision(split_list_start_index, DEFAULT_POSITION_BYTE_SIZE_THRESHOLD);

        }


#warning eventually update data logger stuff
        data_logger_.add_evolution_point({(to_expand_count() - expanded_count()),
                                          expanded_count(),
                                          (int)mesh_.n_logical_vertices(),
                                          (int)cluster.cone_tip_vertices().size(),
                                          (int)cluster.cone_tip_vertices().size(),
                                          last_single_sv_ec_count_,
                                          0,
                                          last_unexp_valence_,
                                          PEHelpers::DBCI_vertices_count(mesh_),
                                          0});

        last_single_sv_ec_count_ = 0;
        last_unexp_valence_ = 0;

#warning until here


#if 0
        auto updated_edge_prop = codomain_mesh_.request_edge_property<bool>();

        for(auto extended_cone_tip_vertex: cluster_submesh.cone_tip_vertices()) {
            auto mesh_tip_vertex = cluster_submesh.cone_to_mesh_handle(extended_cone_tip_vertex);

            for(auto out_he: codomain_mesh_.outgoing_halfedges(mesh_tip_vertex)) {
                auto eh = codomain_mesh_.edge_handle(out_he);

                if(!updated_edge_prop[eh]){

                    auto e = codomain_mesh_.edge(eh);

                    auto from_vertex = e.from_vertex();
                    auto to_vertex   = e.to_vertex();

                    auto from_vertex_cluster_index = cluster_index_prop_[from_vertex];
                    auto to_vertex_cluster_index   = cluster_index_prop_[to_vertex];

                    if ((connecting_cluster_indices_prop_[eh].first > 0 && from_vertex_cluster_index != connecting_cluster_indices_prop_[eh].first) ||
                            (connecting_cluster_indices_prop_[eh].second > 0 && to_vertex_cluster_index != connecting_cluster_indices_prop_[eh].second)) {
                        std::cout << " ERROR - edge " << codomain_mesh_.edge(eh) <<
                                  " is already connecting clusters " << connecting_cluster_indices_prop_[eh].first <<
                                  " and " << connecting_cluster_indices_prop_[eh].second << std::endl;
                        return -1;
                    }

                    connecting_cluster_indices_prop_[codomain_mesh_.edge_handle(out_he)] = {from_vertex_cluster_index,
                                                                                   to_vertex_cluster_index};

                    //PRINT_IF_NOT_SILENT(" - set edge "<<codomain_mesh_.edge(eh)<<" as connecting clusters "<<from_vertex_cluster_index<<" and "<<to_vertex_cluster_index<<std::endl);

                    updated_edge_prop[eh] = true;
                }
            }
        }
#endif

        PRINT_IF_NOT_SILENT(" ======== CLUSTER EXPANSION RESULT: "<<cluster_expansion_result<<std::endl);

#if ENABLE_ALL_CHECKS
        auto degenerate_tets_count = ExactBadTetFinder::findBadTets(mesh_, vertex_position_prop_).first.size();
        PRINT_IF_NOT_SILENT(" ======== updated degenerates count: "<<degenerate_tets_count<<std::endl);
        if(ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)){
            PRINT_IF_NOT_SILENT(" ERROR: Found flipped tets after expanding cluster"<<std::endl);
            return -1;
        }
#endif
        PRINT_IF_NOT_SILENT(" ============================================================"<<std::endl);
        PRINT_IF_NOT_SILENT(" ============================================================"<<std::endl);



#warning high-jacking EC-break data logging field
        last_single_sv_ec_count_ = cluster.cone_tip_vertices().size();

        return cluster_expansion_result;

    }





    int Expander::find_expandable_cluster_and_expand(bool enable_star_shapification,
                                                     int max_k,
                                                     int max_n_vertices){
        PRINT_IF_NOT_SILENT(" ---------"<<std::endl);
        PRINT_IF_NOT_SILENT(" trying to merge unexpandable cones with max k = "<<max_k<<" and max #vertices = "<<max_n_vertices<<std::endl);

        LightWeightStopWatch stop_watch;

        bool expanded_at_least_one_cluster(false);

        auto unexpanded_vertices_index_prop = mesh_.request_vertex_property<int>();
        auto candidate_vertices_prop = mesh_.request_vertex_property<bool>();

        std::vector<VertexHandle> unexpanded_vertices;
        //gather unexpanded vertices
        int positive_exp_valence_unexp_vertices_count(0);
        int unexp_vertices_count(0);
        int computed_cones_count(0);
        for(auto v: mesh_.vertices()){
            if(!expanded_prop_[v]){
                unexp_vertices_count++;

                if(expansion_valence(v)) {
                    positive_exp_valence_unexp_vertices_count++;
                    unexpanded_vertices.push_back(v);



                    candidate_vertices_prop[v] = true;

                    //if(full_cone_at_last_iteration_prop_[v].first.cone_to_mesh_handle(VertexHandle(0)) != v) {
                    if(expanse_at_last_iteration_prop_[v].expansion_cone.cone_to_mesh_handle(VertexHandle(0)) != v) {
                        computed_cones_count++;
                        ExpansionCone cone;
                        auto result = gather_full_expansion_cone(v, cone);
                        if (result) {
                            std::cout << " ERROR - could not gather cone for vertex " << v << std::endl;
                            return -1;
                        }
                    }

                    unexpanded_vertices_index_prop[v] = unexpanded_vertices.size()-1;
                }
            }
        }
        PRINT_IF_NOT_SILENT(" -      found " << unexp_vertices_count);
        PRINT_IF_NOT_SILENT(" unexpanded vertices, among which "<<positive_exp_valence_unexp_vertices_count);
        PRINT_IF_NOT_SILENT(" have positive exp-valence"<< std::endl);
        PRINT_IF_NOT_SILENT(" -   computed "<<computed_cones_count<<" cones"<<std::endl);

        const int cluster_size_hard_limit(DEFAULT_CLUSTER_SIZE_HARD_LIMIT);
        PRINT_IF_NOT_SILENT(" WARNING - hard-coded max cluster size = "<<cluster_size_hard_limit<<std::endl);
        const int max_cluster_size(std::min(max_k, std::min((int)unexpanded_vertices.size()/2, cluster_size_hard_limit)));
        const int unexp_v_count(unexpanded_vertices.size());

        PRINT_IF_NOT_SILENT(" - max cluster size = "<<max_cluster_size<<std::endl);
        PRINT_IF_NOT_SILENT(" - unexp_v_count = "<<unexp_v_count<<std::endl);

        bool found_expandable_cluster(false);

        //initial, 1-clusters
        auto* p_previous_iterator = new ConnectedVertexSubsetRecursiveIterator(1, mesh_, candidate_vertices_prop, nullptr, false);
        while(p_previous_iterator->is_valid()){
            p_previous_iterator->next();
        }

        if(!p_previous_iterator){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't allocate 1-cluster iterator"<<std::endl);
            return -1;
        }

        //outer loop for the number of cones to cluster
        for(int cluster_size(2); cluster_size <= max_cluster_size; cluster_size++){

            if(found_expandable_cluster){
                break;
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            std::vector<int> expansion_result_count(EXPANSION_CHECK_RESULT_COUNT, 0);

            PRINT_IF_NOT_SILENT(" ------------- size "<<cluster_size<<" clusters: "<<std::endl);

            int expandable_cluster_count(0);
            int topo_expandable_cluster_count(0);
            int connected_cluster_count(0);

            long int iteration_count(0);

            //TODO: set a limit based on the cluster and unexpanded set sizes
            const long int safety_guard = std::numeric_limits<long int>::max();


            auto p_cluster_iterator = new ConnectedVertexSubsetRecursiveIterator(cluster_size,
                                                                                 mesh_,
                                                                                 candidate_vertices_prop,
                                                                                 p_previous_iterator,
                                                                                 false);

            if(!p_cluster_iterator){
                PRINT_IF_NOT_SILENT(" ERROR - couldn't allocate "<<cluster_size<<"-cluster iterator"<<std::endl);
                return -1;
            }

            //PRINT_IF_NOT_SILENT(" -- max for first index : "<<(unexp_v_count - max_cluster_size)<<std::endl);
            while(iteration_count < safety_guard &&
                  !found_expandable_cluster &&
                  p_cluster_iterator->is_valid()){

                iteration_count++;

                check_for_timeout();

                auto cluster_vertices = p_cluster_iterator->next();

                std::vector<int> cluster_indices;
                for(auto v: cluster_vertices){
                    cluster_indices.push_back(unexpanded_vertices_index_prop[v]);
                }


                //PRINT_IF_NOT_SILENT(" --- "<<cluster_vertices<<std::endl);

                //check cluster
                ExpansionCone cluster;
                auto merge_start = std::chrono::high_resolution_clock::now();
                int cluster_result = merge_expansion_cones(unexpanded_vertices,
                                                           cluster_indices,
                                                           cluster);
                auto merge_end = std::chrono::high_resolution_clock::now();
                float merge_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(merge_end - merge_start).count() / 1000000;


                if(cluster_result){
                    PRINT_IF_NOT_SILENT(" ERROR - couldn't merge cluster ");
                    for(auto c: cluster_indices){
                        PRINT_IF_NOT_SILENT(unexpanded_vertices[c]<<" ");
                    }
                    PRINT_IF_NOT_SILENT(std::endl);
                    return -1;
                }


                connected_cluster_count++;



                //VertexHandle tip_vertex(0);
                VertexPosition new_position;
                auto exp_check_start = std::chrono::high_resolution_clock::now();
                auto cluster_expandability = cluster.is_expandable(new_position);
                auto exp_check_end = std::chrono::high_resolution_clock::now();
                float exp_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(exp_check_end - exp_check_start).count() / 1000000;

                check_for_timeout();


                if(merge_duration_s > 1 || exp_check_duration_s > 1){
                    PRINT_IF_NOT_SILENT(" ------------------------------ checked cluster "<<cluster_vertices<<std::endl);
                    PRINT_IF_NOT_SILENT(" - cluster merge duration: "<<merge_duration_s<<std::endl);
                    PRINT_IF_NOT_SILENT(" - expandability check duration: "<<exp_check_duration_s<<std::endl);
                    int max, total;
                    double avg;
                    VertexHandle max_precision_vh;
                    compute_precision_stats(cluster, cluster.vertex_position_prop(), max, total, avg, max_precision_vh);
                    PRINT_IF_NOT_SILENT(" - cluster size = "<<cluster.n_vertices()<<
                                        ", max precision = "<<max<<
                                        ", total precision = "<<total<<
                                        ", average precision = "<<avg<<std::endl);

                    data_logger_.add_cluster_computation_point((merge_duration_s + exp_check_duration_s)*1e6,
                                                               max,
                                                               total,
                                                               cluster.n_vertices());
                }

                if(cluster_expandability == IS_NOT_GEO_EXPANDABLE && enable_star_shapification){
                    for(auto tip: cluster.cone_tip_vertices()){
                        auto manifold = TopoHelper::manifoldVertex(cluster, tip);
                        if(!manifold){
                            cluster_expandability = CLUSTER_TIP_IS_NOT_BOUNDARY_2_MANIFOLD;
                            PRINT_IF_NOT_SILENT(" --> one of the tips is actually not boundary 2-manifold"<<std::endl);
                            PRINT_IF_NOT_SILENT(" cluster: "<<cluster<<std::endl);
                            break;
                        }
                    }


                    if((int)cluster.n_vertices() > max_n_vertices){
                        PRINT_IF_NOT_SILENT(" --> cluster has "<<cluster.n_vertices()<<" vertices, which is more than "<<max_n_vertices<<", skipping"<<std::endl);
                        //kind of a hack there
                        cluster_expandability = IS_NOT_TOPO_EXPANDABLE;

                    }
                }


                if(cluster_expandability == IS_NOT_GEO_EXPANDABLE && enable_star_shapification){

                    if(split_edges_between_cone_vertices_but_not_part_of_the_cone(cluster)){
                        std::cout<<" error while splitting 'extra-cone' edges"<<std::endl;
                        return -1;
                    }

                    if(split_edges_of_tets_with_three_faces_on_cone_base(cluster)){
                        std::cout<<" error while splitting 'tricky tets' edges"<<std::endl;
                        return -1;
                    }

                    check_for_timeout();

                    auto cluster_ss_result = star_shapify_cluster(cluster);

                    if(cluster_ss_result){
                        PRINT_IF_NOT_SILENT(" ERROR - couldn't star-shapify cluster "<<cluster<<std::endl);
                        //cluster.print_details(3);

                        return -1;
                    }
                    found_expandable_cluster = true;
                    cluster_star_shapifications_count_++;


                    /*ExpansionCone cluster;
                    cluster_result = merge_expansion_cones(unexpanded_vertices,
                                                           cluster_indices,
                                                           cluster);

                    if(cluster_result){
                        PRINT_IF_NOT_SILENT(" ERROR - couldn't merge cluster ");
                        for(auto c: cluster_indices){
                            PRINT_IF_NOT_SILENT(unexpanded_vertices[c]<<" ");
                        }
                        PRINT_IF_NOT_SILENT(std::endl);
                        return -1;
                    }

                    cluster_expandability = cluster.is_expandable(new_position);

                    if(cluster_expandability){
                        PRINT_IF_NOT_SILENT(" ERROR - cluster is still not expandable after star-shapification, code: "<<cluster_expandability<<std::endl);
                        cluster.print_details();
                        return -1;
                    }*/
                }


                //PRINT_IF_NOT_SILENT(" - cluster "<<cluster<<" expandability: "<<cluster_expandability<<std::endl);
                if (!cluster_expandability) {

                    PRINT_IF_NOT_SILENT(" ==> cluster " << iteration_count << " : ");
                    for (auto c: cluster_indices) {
                        PRINT_IF_NOT_SILENT(unexpanded_vertices[c] << " ");
                    }
                    PRINT_IF_NOT_SILENT(" IS EXPANDABLE. Position: " << vec2vec(new_position) << std::endl);
                    expandable_cluster_count++;

                    if(split_edges_between_cone_vertices_but_not_part_of_the_cone(cluster)){
                        std::cout<<" error while splitting 'extra-cone' edges"<<std::endl;
                        return -1;
                    }

                    if(split_edges_of_tets_with_three_faces_on_cone_base(cluster)){
                        std::cout<<" error while splitting 'tricky tets' edges"<<std::endl;
                        return -1;
                    }

                    auto cluster_expansion_result = expand_cluster(cluster, new_position);
                    if(!cluster_expansion_result){
                        found_expandable_cluster = true;
                    }else{
                        PRINT_IF_NOT_SILENT(" ERROR - couldn't expand cluster... "<<std::endl);
                        return -1;
                    }

                    //temp check
                    /*bool single_cc;
                    int euler_charac;
                    if(check_main_cluster_topology(single_cc, euler_charac)){
                        PRINT_IF_NOT_SILENT(" =====> main cluster is no longer ball topology after expanding cluster "<<cluster<<std::endl);
                        cluster.print_details();
                        cluster.is_topo_expandable(true);
                        return -1;
                    }*/
                    //temp check ends here

                    while(cluster_size_count_.size() <= cluster.cone_tip_vertices().size()){
                        cluster_size_count_.push_back(0);
                    }
                    cluster_size_count_[cluster.cone_tip_vertices().size()]++;
                    break;

                } else if (cluster_expandability == IS_NOT_GEO_EXPANDABLE) {
                    if(enable_star_shapification){
                        PRINT_IF_NOT_SILENT(" --> cluster already expanded during cluster star-shapification"<<std::endl);
                    }else{
                        topo_expandable_cluster_count++;
                    }

                } else if(cluster_expandability == EXPANDABILITY_CHECK_ERROR) {

                    PRINT_IF_NOT_SILENT(" an error occurred while looking for expandable cluster"<<std::endl);
                    found_expandable_cluster = true;
                    return -1;
                }else{

                    expansion_result_count[cluster_expandability]++;

                    //This shouldn't happen because the ConnectedVertexSubsetRecursiveIterator should only
                    //produce cluster that are a single connected component
                    if (cluster_expandability == CLUSTER_TIPS_ARE_NOT_A_SINGLE_CONNECTED_COMPONENT) {

                        std::cout << " ERROR - cluster tips are not a single connected component" << std::endl;
                        std::cout << " cluster: " << cluster << std::endl;
                        std::cout << " tip vertices: " << cluster.cone_tip_vertices() << std::endl;

                        auto result = cluster.is_topo_expandable(true);

                        PRINT_IF_NOT_SILENT(" ---> result = "<<result<<std::endl);

                        return -1;
                    }
                }
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            float duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000000;

            delete p_previous_iterator;
            p_previous_iterator = p_cluster_iterator;

            PRINT_IF_NOT_SILENT(" - done with size-"<<cluster_size<<" clusters. Tried "<<iteration_count<<" of them in "<<duration_s<<" seconds"<<std::endl);
            PRINT_IF_NOT_SILENT("     fully expandable clusters: "<<expandable_cluster_count<<std::endl);
            PRINT_IF_NOT_SILENT("      topo-expandable clusters: "<<topo_expandable_cluster_count<<std::endl);
            PRINT_IF_NOT_SILENT("            connected clusters: "<<connected_cluster_count<<std::endl);
            /*PRINT_IF_NOT_SILENT(" scanned clusters overview: "<<std::endl);
            for(size_t i(0); i<expansion_result_count.size(); i++){
                PRINT_IF_NOT_SILENT(" - "<<i<<" -> "<<expansion_result_count[i]<<std::endl);
            }*/
            PRINT_IF_NOT_SILENT(" -------------------------------"<<std::endl);
            PRINT_IF_NOT_SILENT(" -------------------------------"<<std::endl);

        }

        PRINT_IF_NOT_SILENT(" ... done! Found expandable cluster: "<<found_expandable_cluster<<std::endl);
        PRINT_IF_NOT_SILENT(" ---------"<<std::endl);

        if(found_expandable_cluster) {
            PRINT_IF_NOT_SILENT(" EXPANDED CLUSTER"<<std::endl);
        }
        return found_expandable_cluster ? 0 : 1;
    }



    int Expander::brute_force_find_expandable_cluster_and_expand(){


        PRINT_IF_NOT_SILENT(" BRUTE FORCE CLUSTER EXPANSION IS DEPRECATED "<<std::endl);
        return -1;

        PRINT_IF_NOT_SILENT(" ---------"<<std::endl);
        PRINT_IF_NOT_SILENT(" trying to merge unexpandable cones..."<<std::endl);

        bool expanded_at_least_one_cluster(false);

        std::vector<VertexHandle> unexpanded_vertices;
        //gather unexpanded vertices
        int positive_exp_valence_unexp_vertices_count(0);
        int unexp_vertices_count(0);
        int computed_cones_count(0);
        for(auto v: mesh_.vertices()){
            if(!expanded_prop_[v]){
                unexp_vertices_count++;

                if(expansion_valence(v)) {
                    positive_exp_valence_unexp_vertices_count++;
                    unexpanded_vertices.push_back(v);

                    if(expanse_at_last_iteration_prop_[v].expansion_cone.cone_to_mesh_handle(VertexHandle(0)) != v) {

                        computed_cones_count++;
                        ExpansionCone cone;
                        auto result = gather_full_expansion_cone(v, cone);
                        if (result) {
                            std::cout << " ERROR - could not gather cone for vertex " << v << std::endl;
                            return -1;
                        }
                    }
                }
            }
        }
        PRINT_IF_NOT_SILENT(" -      found " << unexp_vertices_count);
        PRINT_IF_NOT_SILENT(" unexpanded vertices, among which "<<positive_exp_valence_unexp_vertices_count);
        PRINT_IF_NOT_SILENT(" have positive exp-valence"<< std::endl);
        PRINT_IF_NOT_SILENT(" -   computed "<<computed_cones_count<<" cones"<<std::endl);



        const int cluster_size_hard_limit(DEFAULT_CLUSTER_SIZE_HARD_LIMIT);
        PRINT_IF_NOT_SILENT(" WARNING - hard-coded max cluster size = "<<cluster_size_hard_limit<<std::endl);
        const int max_cluster_size(std::min((int)unexpanded_vertices.size()/2, cluster_size_hard_limit));
        const int unexp_v_count(unexpanded_vertices.size());

        PRINT_IF_NOT_SILENT(" - max cluster size = "<<max_cluster_size<<std::endl);
        PRINT_IF_NOT_SILENT(" - unexp_v_count = "<<unexp_v_count<<std::endl);

        bool found_expandable_cluster(false);

        //outer loop for the number of cones to cluster
        for(int cluster_size(4); cluster_size <= max_cluster_size; cluster_size++){

            if(found_expandable_cluster){
                break;
            }

            std::vector<int> expansion_result_count(EXPANSION_CHECK_RESULT_COUNT, 0);

            PRINT_IF_NOT_SILENT(" ------------- size "<<cluster_size<<" clusters: "<<std::endl);

            int expandable_cluster_count(0);
            int topo_expandable_cluster_count(0);
            int connected_cluster_count(0);

            //cluster initialization (indices)
            std::vector<int> cluster_indices(cluster_size);
            for(int j(0); j<(int)cluster_indices.size(); j++){
                cluster_indices[j] = j;
            }

            long int iteration_count(0);

            //TODO: set a limit based on the cluster and unexpanded set sizes
            const long int safety_guard = std::numeric_limits<long int>::max();

            //PRINT_IF_NOT_SILENT(" -- max for first index : "<<(unexp_v_count - max_cluster_size)<<std::endl);
            while(iteration_count < safety_guard && !found_expandable_cluster){
                iteration_count++;

                //check cluster
                ExpansionCone cluster = ExpansionCone();

                int cluster_result = merge_expansion_cones(unexpanded_vertices,
                                                           cluster_indices,
                                                           cluster);

                if(cluster_result == -1){
                    PRINT_IF_NOT_SILENT(" ERROR - couldn't merge cluster ");
                    for(auto c: cluster_indices){
                        PRINT_IF_NOT_SILENT(c<<" ");
                    }
                    PRINT_IF_NOT_SILENT(std::endl);
                }else if(!cluster_result) {

                    connected_cluster_count++;
                    /*PRINT_IF_NOT_SILENT(" --- ";
                    for(auto idx: cluster_indices){
                        PRINT_IF_NOT_SILENT(" "<<unexpanded_vertices[idx];
                    }
                    PRINT_IF_NOT_SILENT(std::endl;*/

                    VertexPosition new_position;
                    auto cluster_expandability = cluster.is_expandable(new_position);

                    if (!cluster_expandability) {
                        PRINT_IF_NOT_SILENT(" ==> cluster " << iteration_count << " : ");
                        for (auto c: cluster_indices) {
                            PRINT_IF_NOT_SILENT(c << " ");
                        }
                        PRINT_IF_NOT_SILENT(" IS EXPANDABLE. Position: " << vec2vec(new_position) << std::endl);
                        expandable_cluster_count++;

                        auto cluster_expansion_result = expand_cluster(cluster, new_position);
                        if(!cluster_expansion_result){
                            found_expandable_cluster = true;
                        }else{
                            PRINT_IF_NOT_SILENT(" ERROR - couldn't expand cluster... "<<std::endl);
                            return -1;
                        }

                        //temp check
                        /*bool single_cc;
                        int euler_charac;
                        if(check_main_cluster_topology(single_cc, euler_charac)){
                            PRINT_IF_NOT_SILENT(" =====> main cluster is no longer ball topology after expanding cluster "<<cluster<<std::endl);
                            cluster.print_details();
                            cluster.is_topo_expandable(true);
                            return -1;
                        }*/

                        while(cluster_size_count_.size() <= cluster.cone_tip_vertices().size()){
                            cluster_size_count_.push_back(0);
                        }
                        cluster_size_count_[cluster.cone_tip_vertices().size()]++;

                        break;

                    } else if (cluster_expandability == IS_NOT_GEO_EXPANDABLE) {
                        PRINT_IF_NOT_SILENT(" ==> cluster " << iteration_count << " : ");
                        for (auto c: cluster_indices) {
                            PRINT_IF_NOT_SILENT(c << " ");
                        }
                        PRINT_IF_NOT_SILENT(" IS TOPO-EXPANDABLE." << std::endl);
                        topo_expandable_cluster_count++;
                    } else if(cluster_expandability == EXPANDABILITY_CHECK_ERROR) {

                        PRINT_IF_NOT_SILENT(" an error occurred while looking for expandable cluster"<<std::endl);
                        return -1;
                    }else{

                        expansion_result_count[cluster_expandability]++;
                    }
                }



                //update cluster indices
                if(cluster_indices[0] == (unexp_v_count - cluster_size)){
                    break;
                }

                for(int j(0); j<cluster_size; j++){
                    //PRINT_IF_NOT_SILENT(" -- j = "<<j<<std::endl);
                    //PRINT_IF_NOT_SILENT(" -- cluster["<<cluster_size - 1 - j<<"] = "<<cluster_indices[cluster_size - 1 - j]<<std::endl);

                    if(cluster_indices[cluster_size - 1 - j] < unexp_v_count - 1 - j){
                        cluster_indices[cluster_size - 1 - j]++;
                        for(int k(cluster_size - 1 - j + 1); k<cluster_size; k++){
                            cluster_indices[k] = cluster_indices[k - 1] + 1;
                        }
                        break;
                    }else{

                    }
                }

            }

            PRINT_IF_NOT_SILENT(" - done with size-"<<cluster_size<<" clusters. Tried "<<iteration_count<<" of them"<<std::endl);
            PRINT_IF_NOT_SILENT("     fully expandable clusters: "<<expandable_cluster_count<<std::endl);
            PRINT_IF_NOT_SILENT("      topo-expandable clusters: "<<topo_expandable_cluster_count<<std::endl);
            PRINT_IF_NOT_SILENT("            connected clusters: "<<connected_cluster_count<<std::endl);
            PRINT_IF_NOT_SILENT(" scanned clusters overview: "<<std::endl);
            for(size_t i(0); i<expansion_result_count.size(); i++){
                PRINT_IF_NOT_SILENT(" - "<<i<<" -> "<<expansion_result_count[i]<<std::endl);
            }
            PRINT_IF_NOT_SILENT(" -------------------------------"<<std::endl);
            PRINT_IF_NOT_SILENT(" -------------------------------"<<std::endl);

        }

        PRINT_IF_NOT_SILENT(" ... done!"<<std::endl);
        PRINT_IF_NOT_SILENT(" ---------"<<std::endl);

        if(!found_expandable_cluster) {

            PRINT_IF_NOT_SILENT("--------- unexpanded vertices stuff: " << std::endl);
            for (auto v: mesh_.vertices()) {
                if (!expanded_prop_[v]) {
                    PRINT_IF_NOT_SILENT(" - unexpanded vertex " << v << ": " << std::endl);

                }
            }

            PRINT_IF_NOT_SILENT(" --------- vertices out of the primary cluster: "<<std::endl);
            for (auto v: mesh_.vertices()) {
                if(cluster_index_prop_[v] > 0) {
                    PRINT_IF_NOT_SILENT(" - vertex " << v << ", expanded: " << expanded_prop_[v] << ", in cluster: "<< cluster_index_prop_[v] << std::endl);
                }
            }

        }
        return found_expandable_cluster ? 0 : 1;

    }








    int Expander::break_best_unexpandable_cone_with_memory(VertexExpanse& best_expanse,
                                                           bool stop_at_first_feasible_max_cone){

        PRINT_IF_NOT_SILENT(" ---------"<<std::endl);
        PRINT_IF_NOT_SILENT(" breaking best expandable cone with memory..."<<std::endl);

        //PRINT_IF_NOT_SILENT(" - stopping at first feasible max cone: "<<stop_at_first_feasible_max_cone<<std::endl);

        int previously_breakable_and_now_expanded_count(0);
        std::vector<VertexHandle> unexpanded_vertices;
        //gather unexpanded vertices
        for(auto v: mesh_.vertices()){
            if(!expanded_prop_[v]){
                unexpanded_vertices.push_back(v);

                //if the vertex is now expanded but there was a single-steiner vertex
                //EC around it previously then we remove this expanse
            }else if(expanse_at_last_iteration_prop_[v].center_vertex == v){
                expanse_at_last_iteration_prop_[v] = VertexExpanse();
                unexpanded_neighbors_at_expanse_computation_[v].clear();
                previously_breakable_and_now_expanded_count++;
            }
        }
        //PRINT_IF_NOT_SILENT(" cones to break: "<<unexpanded_vertices.size()<<std::endl);
        //PRINT_IF_NOT_SILENT(" previously breakable and now expanded count: "<<previously_breakable_and_now_expanded_count<<std::endl);



        int computed_cones_count(0);

        //gather all unmodified expanses from last turn and
        //re-compute the modified ones
        for(auto v: unexpanded_vertices){

            //first, check if this vertex already has a
            //potential expanse around it from last iteration
            if(expanse_at_last_iteration_prop_[v].center_vertex == v){

                //PRINT_IF_NOT_SILENT(" - vertex "<<v<<" already has EC "<<expanse_at_last_iteration_prop_[v].expansion_cone<<std::endl);

                //if yes, we check that the neighborhood of the vertex hasn't changed
                bool found_one_neighbor_expanded_this_turn(false);
                for(auto neighbor: unexpanded_neighbors_at_expanse_computation_[v]){
                    if(expanded_prop_[neighbor]){
                        found_one_neighbor_expanded_this_turn = true;
                        //PRINT_IF_NOT_SILENT(" --> neighbor "<<neighbor<<" was unexpanded before and is now expanded"<<std::endl);
                        break;
                    }
                }

                //and if yes or if we were stuck before, then we have to re-build its expansion cone
                if(found_one_neighbor_expanded_this_turn ||
                   !stop_at_first_feasible_max_cone){
                    expanse_at_last_iteration_prop_[v] = VertexExpanse();
                    unexpanded_neighbors_at_expanse_computation_[v].clear();
                    //PRINT_IF_NOT_SILENT(" --> but neighborhood has changed, recomputing..."<<std::endl);
                    //but if not, then we can skip it since it's still the same
                }else{
                    //PRINT_IF_NOT_SILENT(" --> and neighborhood hasn't changed, skipping"<<std::endl);
                    continue;
                }
            }

            //if this vertex doesn't have an expanse from last turn
            //or if its neighborhood has changed, then we compute its new expanse

            ExpansionCone cone;
            auto result = find_maximal_expansion_cone(v,
                                                      cone,
                                                      stop_at_first_feasible_max_cone);

            computed_cones_count++;

            //set the expanse at last iteration
            expanse_at_last_iteration_prop_[v] = {v, cone, VertexPosition(0), EXPANSION_RESULT_STATUS(result)};
            for(auto out_he: mesh_.outgoing_halfedges(v)){
                auto neighbor = mesh_.to_vertex_handle(out_he);
                if(!expanded_prop_[neighbor]){
                    unexpanded_neighbors_at_expanse_computation_[v].push_back(neighbor);
                }
            }
            //PRINT_IF_NOT_SILENT(" --> computed EC "<<cone<<" for vertex "<<v<<std::endl);
            //PRINT_IF_NOT_SILENT(" --> computed count = "<<computed_cones_count<<std::endl);

            //PRINT_IF_NOT_SILENT(" ---> "<<unexpanded_neighbors_at_expanse_computation_[v].size()<<" unexpanded neighbors"<<std::endl);


            if(result == -1){
                PRINT_IF_NOT_SILENT(" error while trying to find maximal expansion cone for vertex "<<v<<std::endl);
                return -1;
            }
        }

        //PRINT_IF_NOT_SILENT(" computed "<<computed_cones_count<<" cones."<<std::endl);
        //PRINT_IF_NOT_SILENT(" looking for best expanse"<<std::endl);

        std::vector<int> expansion_result_count(FAILURE_CODES_COUNT);
        int empty_cones_count(0);

        //then we go through all existing expanses to find the best one
        int min_steiner_vertices_count(std::numeric_limits<int>::max());
        std::vector<VertexExpanse> single_steiner_vertex_expanses;

        std::vector<VertexExpanse> single_steiner_vertex_and_topo_expanses;
        bool at_least_one_success(false);

        for(auto v: unexpanded_vertices){

            int init_expansion_valence = expansion_valence(v);

            auto expanse = expanse_at_last_iteration_prop_[v];
            expansion_result_count[expanse.result_status]++;

            if(expanse.result_status == EXPANSION_ERROR){
                PRINT_IF_NOT_SILENT( " ERROR - expanse at previous iteration with bad status"<<std::endl);
                return -1;
            }else if(expanse.result_status != EXPANSION_SUCCESS){

                if(!expanse.expansion_cone.n_vertices()){
                    empty_cones_count++;
                }

                continue;
            }

            //-1 to avoid counting the center vertex
            int max_cone_expansion_valence = expanse.expansion_cone.n_vertices() ?
                                             (int)expanse.expansion_cone.n_vertices() -1 : 0;


            auto steiner_vertices_created_count = init_expansion_valence - max_cone_expansion_valence;



            if(steiner_vertices_created_count <= min_steiner_vertices_count &&
               expanse.expansion_cone.n_cells() > 1){
                //PRINT_IF_NOT_SILENT(" -- new expansion cone will create "<<steiner_vertices_created_count<<
                //           " Steiner vertices, updating best expanse as "<<expanse.expansion_cone<<
                //           " around center vertex "<<v<<std::endl);

                best_expanse = expanse;
                min_steiner_vertices_count = steiner_vertices_created_count;
                at_least_one_success = true;

                if(steiner_vertices_created_count == 1){
                    single_steiner_vertex_expanses.push_back(best_expanse);

                    ExpansionCone full_cone;
                    gather_full_expansion_cone(best_expanse.center_vertex, full_cone);
                    VertexPosition pos;
                    if(!full_cone.is_topo_expandable()) {
                        single_steiner_vertex_and_topo_expanses.push_back(best_expanse);
                        PRINT_IF_NOT_SILENT(" --> found single-SV, topo-expandable expanse"<<std::endl);
                    }
                }
            }
        }

        PRINT_IF_NOT_SILENT(" -------------------------------"<<std::endl);
        PRINT_IF_NOT_SILENT(" scanned all ECs, overview: "<<std::endl);
        for(size_t i(0); i<expansion_result_count.size(); i++){
            PRINT_IF_NOT_SILENT(" - "<<i<<" -> "<<expansion_result_count[i]<<std::endl);
        }
        PRINT_IF_NOT_SILENT(" - empty cones: "<<empty_cones_count<<std::endl);
        PRINT_IF_NOT_SILENT(" -------------------------------"<<std::endl);


        const bool use_min_unexp_valence(true);

        //by now, the best expanse should be either the expanse with the lowest number of Steiner vertices
        //or the last one with one Steiner vertex.
        //if it's the latter case, then we take a random one in the list
        if(!single_steiner_vertex_and_topo_expanses.empty()){

            if(use_min_unexp_valence) {
                std::cout << " ---> found " << single_steiner_vertex_and_topo_expanses.size() <<
                          " single-SV, topo-expandable expanses, taking the one with less many unexpanded vertices"
                          << std::endl;

                int min_unexp_valence(std::numeric_limits<int>::max());
                for (const auto &single_sv_exp: single_steiner_vertex_expanses) {
                    auto tip_vertex = single_sv_exp.center_vertex;
                    int unexp_valence = (int) mesh_.valence(tip_vertex) - expansion_valence(tip_vertex);

                    if (unexp_valence < min_unexp_valence) {
                        min_unexp_valence = unexp_valence;
                        best_expanse = single_sv_exp;
                    }
                }
                std::cout << " --- best expanse has " << min_unexp_valence << " unexpanded neghbors" << std::endl;

                last_unexp_valence_ = min_unexp_valence;

            }else {

                std::cout << " ---> found " << single_steiner_vertex_and_topo_expanses.size()
                          << " single-SV, topo-expandable expanses, taking a random one" << std::endl;

                int random_index = std::rand() % single_steiner_vertex_and_topo_expanses.size();

                int unexp_valence =
                        (int) mesh_.valence(best_expanse.center_vertex) - expansion_valence(best_expanse.center_vertex);

                //last_unexp_valence_ = random_index;
                last_unexp_valence_ = unexp_valence;

                best_expanse = single_steiner_vertex_and_topo_expanses[random_index];
                std::cout << " --- random index = " << random_index <<
                          " -> selected cone: " << best_expanse.expansion_cone << std::endl;
            }


        }else if(!single_steiner_vertex_expanses.empty()){

            if(use_min_unexp_valence) {
                std::cout << " ---> found " << single_steiner_vertex_expanses.size() <<
                          " single-SV expanses, taking the one with less many unexpanded vertices"
                          << std::endl;

                int min_unexp_valence(std::numeric_limits<int>::max());
                for (const auto &single_sv_exp: single_steiner_vertex_expanses) {
                    auto tip_vertex = single_sv_exp.center_vertex;
                    int unexp_valence = (int) mesh_.valence(tip_vertex) - expansion_valence(tip_vertex);

                    if (unexp_valence < min_unexp_valence) {
                        min_unexp_valence = unexp_valence;
                        best_expanse = single_sv_exp;
                    }
                }
                std::cout << " --- best expanse has " << min_unexp_valence << " unexpanded neghbors" << std::endl;

                last_unexp_valence_ = min_unexp_valence;

            }else {

                std::cout << " ---> found " << single_steiner_vertex_expanses.size()
                          << " single Steiner vertex expanses, taking a random one" << std::endl;

                //std::srand(time(NULL));
                int random_index = std::rand() % single_steiner_vertex_expanses.size();

                best_expanse = single_steiner_vertex_expanses[random_index];

                int unexp_valence =
                        (int) mesh_.valence(best_expanse.center_vertex) - expansion_valence(best_expanse.center_vertex);

                //last_unexp_valence_ = random_index;
                last_unexp_valence_ = unexp_valence;


                std::cout << " --- random index = " << random_index <<
                          " -> selected cone: " << best_expanse.expansion_cone << std::endl;
            }


        }else if(!at_least_one_success){
            PRINT_IF_NOT_SILENT(" - couldn't find a single expansion cone to break"<<std::endl);
            if(!got_stuck_){
                PRINT_IF_NOT_SILENT(" ---> marking Expander as stuck"<<std::endl);

                //PRINT_IF_NOT_SILENT(" WARNING - MANUAL STOP"<<std::endl);
                //return 1;

                got_stuck_ = true;
                //returning as success to run for another iteration but now with the "got stuck" status
                return 0;
            }else{
                PRINT_IF_NOT_SILENT(" ---> already stuck, stopping"<<std::endl);
                best_expanse.result_status = EXPANSION_FAILURE;
                return 1;
            }
        }

        //if we succesfully found something to break, then we're no longer stuck
        got_stuck_ = false;


        last_single_sv_ec_count_ = single_steiner_vertex_expanses.size();
        auto break_result = perform_expansion_with_max_cone(best_expanse.center_vertex,
                                                            best_expanse.expansion_cone);

        PRINT_IF_NOT_SILENT(" --> broke and expanded cone "<<best_expanse.expansion_cone<<std::endl);

        if(break_result){
            PRINT_IF_NOT_SILENT(" ERROR - cannot redo maximal break, result = "<<break_result<<std::endl);
            return -1;
        }



        //PRINT_IF_NOT_SILENT(" ... done, expanded count: "<<expanded_count_<<"/"<<to_expand_count_<<std::endl);
        //PRINT_IF_NOT_SILENT(" ---------"<<std::endl);

        return 0;

    }



    int Expander::n_boundary_vertices() const{
        int count(0);
        for(auto v: mesh_.vertices()){
            count += mesh_.is_boundary(v);
        }
        return count;
    }



    EXPANSION_RESULT_STATUS Expander::expand_best_expandable_vertex_with_memory(bool enable_star_shapification,
                                                                                int target_expansion_valence,
                                                                                int max_expansion_valence){
        //PRINT_IF_NOT_SILENT(" ---------"<<std::endl);
        //PRINT_IF_NOT_SILENT(" expanding best expandable vertex with memory..."<<std::endl);
        //PRINT_IF_NOT_SILENT(" - max exp valence = "<<max_expansion_valence<<std::endl);

        //PRINT_IF_NOT_SILENT(" - star-shapification enabled: "<<enable_star_shapification<<std::endl);

        if(is_fully_expanded()){
            PRINT_IF_NOT_SILENT(" --> already fully expanded"<<std::endl);
            return EXPANSION_SUCCESS;
        }

        LightWeightStopWatch stop_watch;

        std::vector<VertexHandle> unexpanded_vertices;
        //gather unexpanded vertices
        for(auto v: mesh_.vertices()){
            if(!expanded_prop_[v]){
                unexpanded_vertices.push_back(v);

                //if the vertex is now expanded then we remove this expanse
            }else if(expanse_at_last_iteration_prop_[v].center_vertex == v){
                expanse_at_last_iteration_prop_[v] = VertexExpanse();
                unexpanded_neighbors_at_expanse_computation_[v].clear();
                valence_at_expanse_computation_[v] = {0,0};
            }
        }

        if(expanded_count_ == (to_expand_count_-1)){
            PRINT_IF_NOT_SILENT(" - found last unexpanded vertex, which is already expanded"<<std::endl);
            //expanded_prop_[unexpanded_vertices.front()] = true;
            set_expanded_prop(unexpanded_vertices.front(), true);
            fully_expanded_ = true;
            return EXPANSION_SUCCESS;
        }


        int computed_cones_count(0);
        int star_shapifyable_count(0);

        auto set_up_duration_s = stop_watch.lap_duration();

        //std::cout<<" [t] initial set-up duration: "<<set_up_duration_s<<std::endl);

        float average_iteration_duration_s(0);
        int lap_count(0);

        float flow_handling_duration_s(0);
        float one_ring_setup_duration_s(0);
        float expandability_check_duration_s(0);

        //gather all unmodified expanses from last turn and
        //re-compute the modified ones
        for(auto v: unexpanded_vertices){

            LightWeightStopWatch iteration_sw;

            //first, check if this vertex already has a
            //potential expanse around it from last iteration
            if(expanse_at_last_iteration_prop_[v].center_vertex == v){

                //PRINT_IF_NOT_SILENT(" - vertex "<<v<<" already has EC "<<expanse_at_last_iteration_prop_[v].expansion_cone<<
                //           " with expandability status :"<<expanse_at_last_iteration_prop_[v].exp_result<<std::endl);

                //if yes, we check that the neighborhood of the vertex hasn't changed
                bool found_one_modified_neighbor(false);
                bool valence_changed_this_turn(false);


                for(auto neighbor: unexpanded_neighbors_at_expanse_computation_[v]){
                    if(expanded_prop_[neighbor]){
                        found_one_modified_neighbor = true;
                        //PRINT_IF_NOT_SILENT(" --> neighbor "<<neighbor<<" was unexpanded before and is now expanded"<<std::endl);
                        break;
                    }
                }

                int exp_val = expansion_valence(v);
                if(exp_val != valence_at_expanse_computation_[v].first){
                    valence_changed_this_turn = true;
                    valence_at_expanse_computation_[v].first = exp_val;
                }

                int unexp_val = unexpansion_valence(v);
                if(unexp_val != valence_at_expanse_computation_[v].second){
                    valence_changed_this_turn = true;
                    valence_at_expanse_computation_[v].second = unexp_val;
                }

                //PRINT_IF_NOT_SILENT(" warning - manual cone recomputation"<<std::endl);
                //found_one_modified_neighbor = true;

                //and if yes or if we were stuck before, then we have to re-build its expansion cone
                if(found_one_modified_neighbor ||
                        valence_changed_this_turn){

                    expanse_at_last_iteration_prop_[v] = VertexExpanse();
                    unexpanded_neighbors_at_expanse_computation_[v].clear();

                    //PRINT_IF_NOT_SILENT(" --> but neighborhood has changed, recomputing..."<<std::endl);
                    //but if not, then we can skip it since it's still the same
                }else{
                    //PRINT_IF_NOT_SILENT(" --> and neighborhood hasn't changed, skipping"<<std::endl);
                    flow_handling_duration_s += iteration_sw.lap_duration();
                    continue;
                }
            }
            flow_handling_duration_s += iteration_sw.lap_duration();

            //if this vertex doesn't have an expanse from last turn
            //or if its neighborhood has changed, then we compute its new expanse

            //1. gather full cone
            ExpansionCone cone;
            auto full_cone_gathering_result = gather_full_expansion_cone(v, cone);
            if(full_cone_gathering_result){
                PRINT_IF_NOT_SILENT(" ERROR - couldn't gather full cone for vertex "<<v<<std::endl);
                return EXPANSION_ERROR;
            }
            one_ring_setup_duration_s += iteration_sw.lap_duration();

            computed_cones_count++;

            //2. check expandability
            VertexPosition pos;
            auto exp_result = cone.is_expandable(pos);

            auto exp_check_duration = iteration_sw.lap_duration();
            //std::cout<<" -- check took "<<exp_check_duration<<" seconds"<<std::endl);
            //std::cout<<" --------------------"<<std::endl);
            expandability_check_duration_s += exp_check_duration;

            VertexExpanse new_expanse = {v,
                                         cone,
                                         pos,
                                         EXPANSION_SUCCESS,
                                         exp_result};

            //3. if geo-expandable only, then check for star-shapifyability
            if(exp_result){
                //PRINT_IF_NOT_SILENT(" - EC of vertex "<<v<<" is not expandable, code "<<exp_result<<std::endl);
                if(exp_result == IS_NOT_GEO_EXPANDABLE){
                    //PRINT_IF_NOT_SILENT(" --> not geo-expandable, checking star-shapifyability"<<std::endl);
#warning TODO: remove this
                    StarShapifyableExpansionCone ss_cone(cone, 0);
                    auto ss_result = ss_cone.is_star_shapifyable();
                    new_expanse.ss_result = ss_result;

                    if(ss_result){
                        if(ss_result == SS_MIN_VALENCE_GREATER_THAN_TWO){
                            //PRINT_IF_NOT_SILENT(" ---> not star-shapifyable"<<std::endl);
                        }else{
                            //PRINT_IF_NOT_SILENT(" ERROR - this shouldn't happen"<<std::endl);
                            return EXPANSION_ERROR;
                        }
                    }else{
                        //PRINT_IF_NOT_SILENT(" ---> star-shapifyable!"<<std::endl);
                        star_shapifyable_count++;
                    }
                }
            }



            //set the expanse at last iteration
            expanse_at_last_iteration_prop_[v] = new_expanse;
            for(auto out_he: mesh_.outgoing_halfedges(v)){
                auto neighbor = mesh_.to_vertex_handle(out_he);
                if(!expanded_prop_[neighbor]){
                    unexpanded_neighbors_at_expanse_computation_[v].push_back(neighbor);
                }
            }

            average_iteration_duration_s += stop_watch.lap_duration();
            lap_count++;
            //PRINT_IF_NOT_SILENT(" --> computed EC "<<cone<<" for vertex "<<v<<std::endl);
            //PRINT_IF_NOT_SILENT(" --> computed count = "<<computed_cones_count<<std::endl);

        }

        /*std::cout<<" - computed cones: "<<computed_cones_count<<std::endl);

        std::cout<<" [t] total memory set-up duration = "<<average_iteration_duration_s<<" s"<<std::endl);
        std::cout<<" [t]       flow-handling duration = "<<flow_handling_duration_s<<" s"<<std::endl);
        std::cout<<" [t]            1-set-up duration = "<<one_ring_setup_duration_s<<" s"<<std::endl);
        std::cout<<" [t] expandability check duration = "<<expandability_check_duration_s<<" s"<<std::endl);

        average_iteration_duration_s /= lap_count;
        std::cout<<" [t] avg memory set-up iteration duration = "<<average_iteration_duration_s<<" s"<<std::endl);*/


        /*PRINT_IF_NOT_SILENT(" ---------------------"<<std::endl);
        PRINT_IF_NOT_SILENT(" computed "<<computed_cones_count<<" cones, among which "<<star_shapifyable_count<<" were star-shapifyable"<<std::endl);
        PRINT_IF_NOT_SILENT(" unexpanded vertices count = "<<unexpanded_vertices.size()<<std::endl);*/
        //PRINT_IF_NOT_SILENT(" looking for best expanse..."<<std::endl);

        //NOTE: if star-shapification is enabled, it's normal that no new cones would be computed
        //since nothing has changed since.
        if(!computed_cones_count && !enable_star_shapification){
            PRINT_IF_NOT_SILENT(" --> no new cones computed"<<std::endl);
            return EXPANSION_FAILURE;
            //PRINT_IF_NOT_SILENT(" ERROR - No new cones computed"<<std::endl);
            //return EXPANSION_ERROR;
        }


        int min_expansion_valence(std::numeric_limits<int>::max());
        int min_new_pos_size(std::numeric_limits<int>::max());
        VertexExpanse best_expanse;

        int expandable_vertices_count(0);
        int tried_vertices_count(0);

        stop_watch.lap_duration();
        average_iteration_duration_s = 0;
        //then we go through all existing expanses to find the best one
        for(auto v: unexpanded_vertices){

            tried_vertices_count++;

            auto expanse = expanse_at_last_iteration_prop_[v];

            if(expanse.result_status == EXPANSION_ERROR){
                PRINT_IF_NOT_SILENT( " ERROR - expanse at previous iteration with bad status"<<std::endl);
                return EXPANSION_ERROR;

                //can't remember why this is here.
            }else if(expanse.result_status != EXPANSION_SUCCESS){

                continue;
            }

            if(expanse.exp_result == EXPANDABILITY_CHECK_ERROR){
                PRINT_IF_NOT_SILENT( " ERROR - expanse at previous iteration with bad status"<<std::endl);
                return EXPANSION_ERROR;

                //can't remember why this is here either.
            }else if(expanse.exp_result != IS_EXPANDABLE &&
                     (expanse.exp_result != IS_NOT_GEO_EXPANDABLE ||
                      expanse.ss_result != SS_SUCCESS)){
                continue;
            }

            //if the cone is not geo-expandable and we're not star-shapifying cones we can skip this one
            if(expanse.exp_result == IS_NOT_GEO_EXPANDABLE &&
                    !enable_star_shapification){
                continue;
            }


            expandable_vertices_count++;


            //if(enable_star_shapification){
                int exp_valence = expansion_valence(expanse.center_vertex);

                if(exp_valence < max_expansion_valence){
                    if(exp_valence < min_expansion_valence){
                        //PRINT_IF_NOT_SILENT(" --> updated best expandable vertex to "<<v<<" with exp valence = "<<exp_valence<<std::endl);

                        min_expansion_valence = exp_valence;
                        best_expanse = expanse;

                        if(enable_star_shapification &&
                                exp_valence < target_expansion_valence){
                            //PRINT_IF_NOT_SILENT(" ---> which is lower than the maximum, using this one"<<std::endl);
                            break;
                        }
                    }
                }else{
                    //PRINT_IF_NOT_SILENT(" --> exp valence greater than max, skipping"<<std::endl);
                }

            average_iteration_duration_s += stop_watch.lap_duration();

        }

        //std::cout<<" [t] total cone selection duration = "<<average_iteration_duration_s<<" s"<<std::endl);
        //average_iteration_duration_s /= lap_count;
        //std::cout<<" [t] average cone selection iteration duration = "<<average_iteration_duration_s<<" s"<<std::endl);


        //std::cout<<" - tried "<<tried_vertices_count<<" unexpanded vertices"<<std::endl);


        if(best_expanse.center_vertex.idx() == -1){
            PRINT_IF_NOT_SILENT(" - found no vertex to expand"<<std::endl);
            return EXPANSION_FAILURE;
        }


        //then perform the expanse
        /*PRINT_IF_NOT_SILENT(" --> found best cone for vertex "<<best_expanse.center_vertex<<" with expansion valence = "<<expansion_valence(best_expanse.center_vertex)<<std::endl);
        PRINT_IF_NOT_SILENT(" tried "<<tried_vertices_count<<" vertices, among which "<<expandable_vertices_count<<" were expandable"<<std::endl);
        PRINT_IF_NOT_SILENT(" (remaining unexpanded vertices: ("<<unexpanded_vertices.size()<<")"<<std::endl);*/


        /*PRINT_IF_NOT_SILENT(" ---> best vertex is "<<best_expanse.center_vertex<<" with valences (val, unexp, exp) = "<<codomain_mesh_.valence(best_expanse.center_vertex)<<", "<<unexpansion_valence(best_expanse.center_vertex)<<", "<<expansion_valence(best_expanse.center_vertex)<<std::endl);
            PRINT_IF_NOT_SILENT(" ---> best cone: "<<best_expanse.expansion_cone<<std::endl);
            PRINT_IF_NOT_SILENT(" ---> best out of "<<expandable_vertices_count<<" expandable vertices found and "<<tried_vertices_count<<" vertices tried"<<std::endl);
            PRINT_IF_NOT_SILENT(" ---> result status: "<<best_expanse.result_status<<std::endl);
            PRINT_IF_NOT_SILENT(" ---> expansability: "<<best_expanse.exp_result<<std::endl);
            PRINT_IF_NOT_SILENT(" ---> star-shapifyability: "<<best_expanse.ss_result<<std::endl);
            PRINT_IF_NOT_SILENT(" ---> position = "<<best_expanse.new_tip_vertex_pos<<std::endl);*/


        stop_watch.lap_duration();

        int split_list_start_index = split_list_.size();
#warning make this a function
        if(best_expanse.exp_result == IS_NOT_GEO_EXPANDABLE){
            //PRINT_IF_NOT_SILENT(" - best cone is not geo-expandable"<<std::endl);

            if(best_expanse.ss_result != SS_SUCCESS){
                PRINT_IF_NOT_SILENT(" ERROR -  BEST EC is not geo-expandable AND not star-shapifyable"<<std::endl);
                return EXPANSION_ERROR;
            }

            PRINT_IF_NOT_SILENT(" - cone is not geo-expandable, star-shapifying it..."<<std::endl);

            //recompute expansion cone because the positions might have changed due to the smoothing
            /*auto ec_gather_result = gather_full_expansion_cone(best_expanse.center_vertex,
                                                               best_expanse.expansion_cone);*/


            if(split_edges_between_cone_vertices_but_not_part_of_the_cone(best_expanse.expansion_cone)){
                std::cout<<" error while splitting 'extra-cone' edges"<<std::endl;
                return EXPANSION_ERROR;
            }

            if(split_edges_of_tets_with_three_faces_on_cone_base(best_expanse.expansion_cone)){
                std::cout<<" error while splitting 'tricky tets' edges"<<std::endl;
                return EXPANSION_ERROR;
            }


            StarShapifyableExpansionCone ss_cone(best_expanse.expansion_cone,
                                                 remaining_seconds_before_timeout());

            //auto cone_copy = ss_cone;


            PRINT_IF_NOT_SILENT(" WARNING - STAR-SHAPIFIYING WITH SAFE METHOD ONLY"<<std::endl);
            int ss_result = 1;
            //int ss_result = ss_cone.star_shapify(true);
            /*if(ss_result == -1){
                PRINT_IF_NOT_SILENT(" error while running star-shapification"<<std::endl);
                return EXPANSION_ERROR;
            }else */if(ss_result){
                //PRINT_IF_NOT_SILENT(" --> cone couldn't be star-shapified with fast method, using safe one"<<std::endl);

                //this should reset the mesh
                //std::cout<<" TODO: try to fix that mesh reset thing"<<std::endl);

                //ss_cone = StarShapifyableExpansionCone(best_expanse.expansion_cone,
                //                                       remaining_seconds_before_timeout());

                //and then we can star-shapify from scratch again
                ss_result = ss_cone.star_shapify(false, true);

                if(ss_result){

                    PRINT_IF_NOT_SILENT(" -- SS with precision reduction failed, trying without"<<std::endl);
                    ss_cone = StarShapifyableExpansionCone(best_expanse.expansion_cone,
                                                           remaining_seconds_before_timeout());
                    ss_result = ss_cone.star_shapify(false, false);
                }


            }else{
                PRINT_IF_NOT_SILENT(" --> cone could be star-shapified with the fast method"<<std::endl);
            }


            //export_cone(cone_copy, ss_cone);


            if(ss_result){
                if(ss_result == SS_NO_VERTICES_TO_COLLAPSE){
                    PRINT_IF_NOT_SILENT("  WARNING - found cone where all vertices are neighbors to a single witness but is not geo-expandable"<<std::endl);
                    return EXPANSION_ERROR;
                }else if(ss_result == SS_MIN_VALENCE_GREATER_THAN_TWO){

                }else{
                    PRINT_IF_NOT_SILENT(" - ERROR: star-shapification failed"<<std::endl);

                    /*IO::FileManager file_manager;
                    ss_cone.collect_garbage();
                    TetrahedralMesh cone_copy = ss_cone;
                    file_manager.writeFile("bad_EC_ss.ovm", cone_copy);*/
                    return EXPANSION_ERROR;
                }
            }else{


                auto submesh_to_mesh_prop = ss_cone.cone_to_mesh_v_handle_prop();
                int ss_app_result = apply_split_list(ss_cone.get_split_list(),
                                                     submesh_to_mesh_prop);

                if(ss_app_result){
                    PRINT_IF_NOT_SILENT(" - ERROR: couldn't apply star-shapification to mesh"<<std::endl);
                    if(ss_cone.n_cells() > 50){
                        PRINT_IF_NOT_SILENT(" --> cone contains "<<ss_cone.n_cells()<<" cells. Change this line if you really want to print it"<<std::endl);
                    }else{
                        ss_cone.print_details();
                    }

                    /*IO::FileManager file_manager;
                    ss_cone.collect_garbage();
                    TetrahedralMesh cone_copy = best_expanse.expansion_cone;
                    file_manager.writeFile("bad_EC.ovm", cone_copy);

                    cone_copy = ss_cone;
                    file_manager.writeFile("bad_EC_ss.ovm", cone_copy);*/
                    return EXPANSION_ERROR;
                }

#if ENABLE_ALL_CHECKS

                if(ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(mesh_, vertex_position_prop_, best_expanse.center_vertex)){
                    PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after applying star-shapification"<<std::endl);
                    return EXPANSION_ERROR;
                }
#endif


                PRINT_IF_NOT_SILENT(" -- applied star-shapification, re-computing cone for vertex "<<best_expanse.center_vertex<<std::endl);

                ExpansionCone cone;
                auto result = gather_full_expansion_cone(best_expanse.center_vertex, cone);
                PRINT_IF_NOT_SILENT(" -- new cone: "<<cone<<std::endl); //cone.print_details();
                best_expanse.expansion_cone = cone;

                if(result == -1){
                    PRINT_IF_NOT_SILENT(" error while gathering expansion cone"<<std::endl);
                    return EXPANSION_ERROR;
                }

                //cone.print_details();

                VertexPosition new_pos;
                auto expandability_result = cone.is_expandable(new_pos);
                best_expanse.new_tip_vertex_pos = new_pos;
                PRINT_IF_NOT_SILENT(" -- expandability result: "<<expandability_result<<" with pos = "<<vec2vec(new_pos)<<std::endl);

                if(expandability_result){
                    PRINT_IF_NOT_SILENT(" ERROR - cone is still not expandable (code "<<expandability_result<<") after star-shapification"<<std::endl);

                    PRINT_IF_NOT_SILENT(" new cone details: "<<std::endl);
                    cone.print_details(2);

                    VertexPosition new_pos;
                    cone.is_geo_expandable(new_pos, true);

                    return EXPANSION_ERROR;
                }


                total_saved_position_bytes_ += ss_cone.get_saved_position_bytes();
                star_shapification_successes_count_++;
                PRINT_IF_NOT_SILENT(" =========================================================================================================="<<std::endl);
                PRINT_IF_NOT_SILENT(" ================================================= SUCCESFULLY STAR-SHAPIFIED EXPANSION CONE FOR VERTEX "<<best_expanse.center_vertex<<std::endl);
                PRINT_IF_NOT_SILENT("  - saved position bytes = "<<ss_cone.get_saved_position_bytes()<<std::endl);
                PRINT_IF_NOT_SILENT("  - updated total saved position bytes = "<<total_saved_position_bytes_<<std::endl);
                PRINT_IF_NOT_SILENT(" star-shapification successes: "<<star_shapification_successes_count_<<std::endl);
                PRINT_IF_NOT_SILENT(" expansion status: "<<(expanded_count()+1)<<
                           "/"<<to_expand_count()<<
                           " | diff = "<<(to_expand_count() - expanded_count()-1)<<
                           ", current vertices increase ratio: "<<((float)mesh_.n_logical_vertices()/(float)initial_n_vertices_)<<std::endl);

                PRINT_IF_NOT_SILENT(" =========================================================================================================="<<std::endl);
            }
        }

        auto expanse_handling_duration = stop_watch.lap_duration();
        //std::cout<<" [t] expanse handling duration: "<<expanse_handling_duration<<std::endl);


        /*if((to_expand_count() - expanded_count()-1) < 100){
            PRINT_IF_NOT_SILENT(" - finalizing expanse... max precision: "<<max_vertex_precision_B_<<std::endl);
        }*/
        auto final_result = finalize_expanse(best_expanse.center_vertex,
                                             best_expanse.expansion_cone,
                                             best_expanse.new_tip_vertex_pos);


        auto expanse_finalization_duration = stop_watch.lap_duration();
        //std::cout<<" [t] expanse finalization duration: "<<expanse_finalization_duration<<std::endl);


        /*if((to_expand_count() - expanded_count()-1) < 100){
            PRINT_IF_NOT_SILENT(" - checking for flipped tets after expansion... max precision: "<<max_vertex_precision_B_<<std::endl);
        }*/

#if ENABLE_ALL_CHECKS

        if(ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(mesh_, vertex_position_prop_, best_expanse.center_vertex)){
            PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after expanding vertex "<<best_expanse.center_vertex<<std::endl);
            auto flipped_tets = ExactBadTetFinder::findBadTets(mesh_, vertex_position_prop_).second;



            for(auto tet: flipped_tets){
                PRINT_IF_NOT_SILENT(" tet "<<tet<<" : "<<std::endl);
                for(auto v: mesh_.get_cell_vertices(tet)){
                    PRINT_IF_NOT_SILENT(" - "<<v<<" at "<<vec2vec(this->vertex(v))<<", expanded: "<<expanded_prop_[v]<<std::endl);
                }
            }

            return EXPANSION_ERROR;
        }
#endif

        auto flipped_tets_check_duration = stop_watch.lap_duration();
        //std::cout<<" [t] flipped tets check duration: "<<flipped_tets_check_duration<<std::endl);

        /*PRINT_IF_NOT_SILENT(" --> no flipped tets after finalizing expanse for vertex "<<best_expanse.center_vertex<<std::endl);
        if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_)){
            PRINT_IF_NOT_SILENT(" ERROR - actually yes: "<<std::endl);
            auto bad_tets = ExactBadTetFinder::findBadTets(codomain_mesh_, vertex_position_prop_);
            for(auto flipped: bad_tets.second){
                PRINT_IF_NOT_SILENT(" - "<<codomain_mesh_.get_cell_vertices(flipped)<<std::endl);
            }
            return EXPANSION_ERROR;
        }*/

        reduce_mid_vertices_precision(split_list_start_index,
                                      DEFAULT_POSITION_BYTE_SIZE_THRESHOLD);


        auto mid_vertices_precision_reduction_duration = stop_watch.lap_duration();
        //std::cout<<" [t] mid-vertices precision reduction duration: "<<mid_vertices_precision_reduction_duration<<std::endl);

        /*if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_)){
            PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after reducing mid-vertices positions"<<std::endl);
            return EXPANSION_ERROR;
        }*/


#if EDGE_COLLAPSE_ENABLED
        if(!expanding_cluster_ && enable_star_shapification){

            /*smooth_mid_vertices_1_ring_neighborhood(split_list_start_index,
                                                    DEFAULT_INTERIOR_SMOOTHING_MAX_ITERATIONS,
                                                    DEFAULT_INTERIOR_SMOOTHING_EPSILON,
                                                    DEFAULT_POSITION_BYTE_SIZE_SMOOTHING_THRESHOLD);*/

            if(collapse_new_edges_opposite_to_unexpanded_vertices()){
                PRINT_IF_NOT_SILENT(" error while collapsing new edges"<<std::endl);
                return EXPANSION_ERROR;
            }
        }
        /*if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_)){
            PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after collapsing new edges"<<std::endl);
            return EXPANSION_ERROR;
        }*/
#endif

        //std::cout<<" [t] total iteration duration = "<<stop_watch.total_duration()<<" s"<<std::endl);


        return EXPANSION_SUCCESS;
    }




    int Expander::gather_full_expansion_cone(const VertexHandle& center_vertex,
                                             ExpansionCone& cone,
                                             bool accumulating_cluster) const {

        //PRINT_IF_NOT_SILENT(" -------"<<std::endl);
        //PRINT_IF_NOT_SILENT(" gathering full expansion cone around vertex "<<center_vertex<<std::endl);
        //erase the current cone
        if(!accumulating_cluster) {
            cone.clear();
        }


        //to force the tip to be the first vertex (convention)
        if(!accumulating_cluster || !cone.mesh_to_cone_handle(center_vertex).is_valid()) {
            cone.add_vertex(center_vertex,
                            this->vertex(center_vertex),
                            true);
        }

        if(accumulating_cluster){
            cone.set_as_tip(center_vertex);
        }

        //PRINT_IF_NOT_SILENT(" ============================= ADDING CELLS "<<std::endl);

        //then add all the neighboring tets consisting of three
        //expanded vertices and the center vertex
        for(auto vc_it = mesh_.vc_iter(center_vertex); vc_it.valid(); vc_it++){
            //auto cell_vertices = codomain_mesh_.get_cell_vertices(*vc_it);
            //PRINT_IF_NOT_SILENT(" - checking cell "<<*vc_it<<": "<<codomain_mesh_.get_cell_vertices(*vc_it)<<std::endl);

            //when accumulating stuff on a cluster, we temporarily set the tip vertices as
            //expanded, as a trick to work with this function.
            //therefore, we don't need 3 expanded vertices + the tip, but
            //4 expanded vertices, since the tips are temporarily considered as expanded
            if(expanded_vertices_count(*vc_it) == (3 + accumulating_cluster)){
                auto cone_cell = add_cell_to_cone(*vc_it, cone);
                if(cone_cell.idx() == -1){
                    PRINT_IF_NOT_SILENT(" error whille adding cell"<<std::endl);
                    return -1;
                }
                //PRINT_IF_NOT_SILENT(" --> added cell "<<*vc_it<<": "<<codomain_mesh_.get_cell_vertices(*vc_it)<<" to the expansion cone as cell "<<cone_cell<<", positions: "<<std::endl);
            }else{
                //PRINT_IF_NOT_SILENT(" -- expanded count = "<<expanded_vertices_count(*vc_it)<<std::endl);
            }
        }

        //PRINT_IF_NOT_SILENT(" current cone: "<<cone<<std::endl);
        //cone.print_details();
        //PRINT_IF_NOT_SILENT(" ============================= ADDING FACES "<<std::endl);

        //then add all neighboring faces consisting of two
        //expanded vertices (that are not already part of the cone)
        //and the center vertex
        for(auto vf_it = mesh_.vf_iter(center_vertex); vf_it.valid(); vf_it++){
            //auto face_vertices = codomain_mesh_.get_halfface_vertices(codomain_mesh_.halfface_handle(*vf_it, 0));
            //PRINT_IF_NOT_SILENT(" - checking face "<<*vf_it<<": "<<codomain_mesh_.get_halfface_vertices(codomain_mesh_.halfface_handle(*vf_it, 0))<<std::endl);

            //same comment as above regardin clusters
            if(expanded_vertices_count(*vf_it) == (2 + accumulating_cluster)){

                //PRINT_IF_NOT_SILENT(" --> expanded face, adding it to the expansion cone"<<std::endl);
                auto cone_face = add_face_to_cone(*vf_it, cone);

                if(cone_face.idx() == -1){
                    PRINT_IF_NOT_SILENT(" error whille adding face"<<std::endl);
                    return -1;
                }
                /*PRINT_IF_NOT_SILENT(" --> added face "<<*vf_it<<
                           ": "<<codomain_mesh_.get_halfface_vertices(codomain_mesh_.halfface_handle(*vf_it, 0))<<
                           " to the expansion cone as face "<<cone_face<<std::endl);*/

            }else{
                //PRINT_IF_NOT_SILENT(" -- expanded count = "<<expanded_vertices_count(*vf_it)<<std::endl);
            }
        }

        //PRINT_IF_NOT_SILENT(" current cone: "<<cone<<std::endl);
        //cone.print_details();
        //PRINT_IF_NOT_SILENT(" ============================= ADDING EDGES "<<std::endl);

        //and finally add all neighboring edges connecting expanded
        //vertices (that are not already part of the cone) and the center vertex
        for(auto out_he_it: mesh_.outgoing_halfedges(center_vertex)){

            auto to_vertex = mesh_.to_vertex_handle(out_he_it);
            //PRINT_IF_NOT_SILENT(" - checking edge "<<codomain_mesh_.find_halfedge(out_he_it)<<std::endl);
            //PRINT_IF_NOT_SILENT(" - to vertex : "<<to_vertex<<std::endl);
           // PRINT_IF_NOT_SILENT(" - expanded : "<<expanded_prop_[to_vertex]<<std::endl);
            if(expanded_prop_[to_vertex]){
                //PRINT_IF_NOT_SILENT(" --> expanded edge, adding it to the expansion cone"<<std::endl);

                auto cone_to_vertex = cone.mesh_to_cone_handle(to_vertex);

                if(cone_to_vertex.idx() == -1){
                    cone_to_vertex = cone.add_vertex(to_vertex, this->vertex(to_vertex));

                    if(cone_to_vertex.idx() == -1){
                        PRINT_IF_NOT_SILENT(" ERROR - couldn't add vertex "<<to_vertex<<" to expansion cone"<<std::endl);
                        return -1;
                    }
                    //PRINT_IF_NOT_SILENT(" added vertex "<<to_vertex<<" to the expansion cone"<<std::endl);
                }

                auto cone_edge = cone.add_edge(center_vertex,
                                               to_vertex);

                if(cone_edge.idx() == -1){
                    //PRINT_IF_NOT_SILENT(" ERROR - couldn't add edge "<<out_he_it<<": "<<codomain_mesh_.find_halfedge(out_he_it)<<" to the expansion cone"<<std::endl);
                    return -1;
                }
                //PRINT_IF_NOT_SILENT(" --> added edge "<<out_he_it<<" to the expansion cone as edge "<<cone_edge<<std::endl);
            }
        }

        cone.collect_garbage();

        /*PRINT_IF_NOT_SILENT(" ...done. Cone contains "<<cone.n_logical_cells()<<" cells, "<<std::endl);
        PRINT_IF_NOT_SILENT("                        "<<cone.n_logical_faces()<<" faces, "<<std::endl);
        PRINT_IF_NOT_SILENT("                        "<<cone.n_logical_edges()<<" edges"<<std::endl);
        PRINT_IF_NOT_SILENT("                    and "<<cone.n_logical_vertices()<<" vertices"<<std::endl);
        PRINT_IF_NOT_SILENT(" -------"<<std::endl);
        */

        return 0;
    }



    void Expander::export_cone(ExpansionCone& cone_before_SS,
                               const ExpansionCone& cone_after_SS,
                               int sv_threshold,
                               int sv_ratio_threshold) const{


        int steiner_vertices_count = cone_after_SS.n_vertices() - cone_before_SS.n_vertices();
        int steiner_vertices_ratio = cone_after_SS.n_vertices() / cone_before_SS.n_vertices();

        if(steiner_vertices_count < sv_threshold && steiner_vertices_ratio < sv_ratio_threshold){
            PRINT_IF_NOT_SILENT(" -> SVs below threshold, skipping export"<<std::endl);
            return;
        }

        cone_before_SS.project_base_stereographically(cone_before_SS.north_pole());

        VertexPosition bbox_low(cone_before_SS.vertex(VertexHandle(1)));
        VertexPosition bbox_high(bbox_low);

        for(auto v: cone_before_SS.vertices()){
            if(!v.idx()){
                continue;
            }
            const auto& pos = cone_before_SS.vertex(v);
            bbox_low[0] = CGAL::min(bbox_low[0], pos[0]);
            bbox_low[1] = CGAL::min(bbox_low[1], pos[1]);

            bbox_high[0] = CGAL::max(bbox_high[0], pos[0]);
            bbox_high[1] = CGAL::max(bbox_high[1], pos[1]);

            /*PRINT_IF_NOT_SILENT(" - vertex "<<v<<" pos = "<<vec2vec(pos)<<std::endl);
            PRINT_IF_NOT_SILENT("             bbox low = "<<vec2vec(bbox_low)<<std::endl);
            PRINT_IF_NOT_SILENT("            bbox high = "<<vec2vec(bbox_high)<<std::endl<<std::endl);*/
        }
        //PRINT_IF_NOT_SILENT(" - bbox  low = "<<vec2vec(bbox_low)<<std::endl);
        //PRINT_IF_NOT_SILENT(" - bbox high = "<<vec2vec(bbox_high)<<std::endl);

        VertexPosition new_tip_pos = {(bbox_high[0] + bbox_low[0]) / 2,
                                      (bbox_high[1] + bbox_low[1]) / 2,
                                      1-CGAL::min(bbox_high[0] - bbox_low[0], bbox_high[1] - bbox_low[1]) / 2};

        cone_before_SS.set_vertex(VertexHandle(0), new_tip_pos);
        //PRINT_IF_NOT_SILENT(" new tip pos = "<<vec2vec(new_tip_pos)<<std::endl);


        cone_before_SS.collect_garbage();
        TetrahedralMesh mesh_copy = cone_before_SS;
        std::string mesh_name = data_logger_.mesh_name() +
                "_idx" + std::to_string(star_shapification_successes_count_) +
                "_nverts" + std::to_string(cone_before_SS.n_vertices()) +
                "_SV" + std::to_string(steiner_vertices_count) +
                "_SVratio" + std::to_string(steiner_vertices_ratio) +
                ".ovm";

        IO::FileManager file_manager;
        file_manager.writeFile(mesh_name, mesh_copy);

        PRINT_IF_NOT_SILENT(" --> exported cone to file "<<mesh_name<<std::endl);


    }

    int Expander::break_expansion_cone_and_expand(const VertexHandle& center_vertex,
                                                  ExpansionCone& maximal_cone,
                                                  bool stop_at_first_feasible_max_cone){

        PRINT_IF_NOT_SILENT(" ---------"<<std::endl);
        PRINT_IF_NOT_SILENT(" - breaking expansion cone "<<maximal_cone<<
                 " of vertex "<<center_vertex<<std::endl);
        //maximal_cone.print_details();



        auto max_expansion_result = find_maximal_expansion_cone(center_vertex,
                                                                maximal_cone,
                                                                stop_at_first_feasible_max_cone);
        if(max_expansion_result){
            return max_expansion_result;
        }



        max_expansion_result = perform_expansion_with_max_cone(center_vertex,
                                                               maximal_cone);
        if(max_expansion_result){
            return max_expansion_result;
        }


        //PRINT_IF_NOT_SILENT(" ....done. It is expanded"<<std::endl);
        //PRINT_IF_NOT_SILENT(" ---------"<<std::endl);
        return 0;
    }




    int Expander::find_maximal_expansion_cone(const VertexHandle& center_vertex,
                                              ExpansionCone& best_maximal_cone,
                                              bool stop_at_first_feasible){


        //PRINT_IF_NOT_SILENT(" -----"<<std::endl);
        //PRINT_IF_NOT_SILENT(" looking for maximal expansion cone around vertex "<<center_vertex<<std::endl);

        //erase the best max cone
        best_maximal_cone = ExpansionCone();

        int best_result(NO_MIN_EXPANSION);

        int exp_valence = expansion_valence(center_vertex);

        int min_steiner_vertices_count = exp_valence;

        int tried_spokes_count(0);

        //PRINT_IF_NOT_SILENT(" expansion valence = "<<exp_valence<<std::endl);

        for(auto spoke_edge: mesh_.outgoing_halfedges(center_vertex)){
            auto neighbor = mesh_.to_vertex_handle(spoke_edge);

            if(expanded_prop_[neighbor]){

                tried_spokes_count++;

                //PRINT_IF_NOT_SILENT(" - testing spoke "<<codomain_mesh_.find_halfedge(spoke_edge)<<std::endl);

                //then take the neighbor edge (connecting an expanded vertex) with the most
                //expanded tets around
                ExpansionCone minimal_cone;
                //auto result = find_best_minimal_expansion_cone(center_vertex,
                //                                               minimal_cone);

                auto result = extract_expansion_cone_around_spoke_edge(spoke_edge,
                                                                       minimal_cone);

                //PRINT_IF_NOT_SILENT(" - minimal cone around spoke "<<codomain_mesh_.find_halfedge(spoke_edge)<<
                //           " : "<<minimal_cone<<std::endl);

                //if there's no minimal cone around this spoke, we skip it
                if(result == -1){
                    PRINT_IF_NOT_SILENT(" error while extracting expansion cone around spoke edge"<<std::endl);
                    return -1;
                }else if(result){
                    PRINT_IF_NOT_SILENT(" couldn't find minimal expansion cone around spoke edge "<<mesh_.halfedge(spoke_edge)<<std::endl);
                    continue;
                }

                //updating the best result if it's the new best result
                if(best_result && best_result < UNEXP_MIN_CONE) {
                    best_result = UNEXP_MIN_CONE;
                }


                //PRINT_IF_NOT_SILENT(" - checking minimal expansion cone expandability..."<<std::endl);
                //check that this minimal expansion cone is indeed expandable
                VertexPosition new_pos;
                auto expansion_result = minimal_cone.is_expandable(new_pos);

                if(expansion_result){
                    //PRINT_IF_NOT_SILENT(" minimal expansion cone is not expandable. error code = "<<expansion_result<<std::endl);
                    if(expansion_result == IS_NOT_GEO_EXPANDABLE){
                        PRINT_IF_NOT_SILENT(" ERROR - minimal expansion cone is not expandable for geometrical reasons"<<std::endl);
                        return -1;
                    }else if(expansion_result == EXPANDABILITY_CHECK_ERROR){
                        PRINT_IF_NOT_SILENT(" an error occurred while looking for maximal expansion cone"<<std::endl);
                        return -1;
                    }
                    continue;
                }

                if(best_result && best_result < MAX_CONE_EQUALS_MIN) {
                    best_result = MAX_CONE_EQUALS_MIN;
                }


                if(ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)){
                    PRINT_IF_NOT_SILENT(" ==> found flipped tets after expanding the minimal cone"<<std::endl);
                    return -1;
                }

                //PRINT_IF_NOT_SILENT(" epxansion valence before growing = "<<expansion_valence(center_vertex)<<std::endl);

                //then grow the expansion cone to find the biggest expandable possible one
                ExpansionCone maximal_cone;
                auto growth_result = grow_maximal_expansion_cone(center_vertex,
                                                                 minimal_cone,
                                                                 maximal_cone);

                //PRINT_IF_NOT_SILENT(" - max cone: "<<maximal_cone<<std::endl);

                if(growth_result == -1){
                    PRINT_IF_NOT_SILENT(" error while growing expansion cone"<<std::endl);
                    PRINT_IF_NOT_SILENT(" center vertex at "<<this->vertex(center_vertex)<<std::endl);
                    return -1;
                }

                if(ExactBadTetFinder::meshContainsFlippedTets(mesh_,
                                                              vertex_position_prop_)){
                    PRINT_IF_NOT_SILENT(" ==> found flipped tets after growing the maximal cone"<<std::endl);
                    return -1;
                }

                //auto steiner_vertices_created_count = exp_valence - (maximal_cone.n_vertices() - 1);


                if(maximal_cone.same_entities_count(minimal_cone)/* &&
                            steiner_vertices_created_count == 1*/){
                    //PRINT_IF_NOT_SILENT(" maximal cone is equal to the minimal one"<<std::endl);
                    //PRINT_IF_NOT_SILENT(" exp valence = "<<exp_valence<<
                    //           ", maximal cone exp valence = "<<(maximal_cone.n_vertices() - 1)<<std::endl);
                    continue;
                }


                if(best_result && best_result < MAX_CONE_SMALLER_THAN_4_CELLS) {
                    best_result = MAX_CONE_SMALLER_THAN_4_CELLS;
                }

                //NOTE: this is necessary to catch cases where e.g. 1 start tet + 2 tets
                //would create the same cone as the min cone around another spoke
                //since this is impossible if we have more than 3 tets, we use this as a necessary condition
                if(maximal_cone.n_cells() < 4){
                    //PRINT_IF_NOT_SILENT(" max cone has less than 4 cells, skipping"<<std::endl);
                    continue;
                }


                int steiner_vertices_count = exp_valence - (maximal_cone.n_vertices() - 1);

                //PRINT_IF_NOT_SILENT(" - found valid max cone creating "<<steiner_vertices_count<<" SVs"<<std::endl);

                if(steiner_vertices_count <= min_steiner_vertices_count &&
                   maximal_cone.n_cells() > best_maximal_cone.n_cells()){

                    //PRINT_IF_NOT_SILENT(" ----> less than min, updating best max cone"<<std::endl);
                    min_steiner_vertices_count = steiner_vertices_count;

                    best_maximal_cone = maximal_cone;

                    best_result = EXPANSION_SUCCESS;

                    if(min_steiner_vertices_count == 1
                        /*stop_at_first_feasible*/){

                        break;
                        //ExpansionCone full_cone = full_cone_at_last_iteration_prop_[center_vertex].first;
                    }
                }

            }
        }

        /*PRINT_IF_NOT_SILENT(" done! Best result: "<<best_result<<
                   ", min SVs: "<<min_steiner_vertices_count<<
                   " after trying "<<tried_spokes_count<<
                   " spokes (exp valence = "<<exp_valence<<")"<<std::endl);
        PRINT_IF_NOT_SILENT(" -----"<<std::endl);
        */
        return best_result;
    }




    int Expander::find_best_minimal_expansion_cone(const VertexHandle& center_vertex,
                                                   ExpansionCone& minimal_cone){


        //PRINT_IF_NOT_SILENT(" looking for best minimal cone around spoke edge"<<std::endl);

        ExpansionCone best_minimal_cone;
        VertexPosition dummy_pos;

        for(auto out_he_it: mesh_.outgoing_halfedges(center_vertex)){
            auto neighbor = mesh_.to_vertex_handle(out_he_it);

            if(expanded_prop_[neighbor]){
                //PRINT_IF_NOT_SILENT(" -----"<<std::endl);
                //PRINT_IF_NOT_SILENT(" - checking spoke edge "<<codomain_mesh_.find_halfedge(out_he_it)<<std::endl);

                /*auto spoke_edge = topo_helper_.halfedge_exists(center_vertex,
                                                               neighbor);

    #warning probably useless check. remove eventually
                if(spoke_edge.idx() == -1){
                    PRINT_IF_NOT_SILENT(" ERROR - edge from center vertex "<<center_vertex<<
                               " to expanded neighbor "<<neighbor<<
                               " doesn't exist"<<std::endl);
                    return -1;
                }*/

                //PRINT_IF_NOT_SILENT(" ----------"<<std::endl);
                //PRINT_IF_NOT_SILENT(" - checking spoke edge "<<codomain_mesh_.find_halfedge(spoke_edge)<<" with valence "<<codomain_mesh_.valence(codomain_mesh_.edge_handle(spoke_edge))<<std::endl);

                auto result = extract_expansion_cone_around_spoke_edge(out_he_it,
                                                                       minimal_cone);

                //PRINT_IF_NOT_SILENT(" - minimal cone around spoke "<<codomain_mesh_.find_halfedge(spoke_edge)<<
                //           " : "<<minimal_cone<<std::endl);
                if(result == -1){
                    PRINT_IF_NOT_SILENT(" error while extracting expansion cone around spoke edge"<<std::endl);
                    return -1;
                }

                //PRINT_IF_NOT_SILENT(" - cone around spoke edge "<<codomain_mesh_.find_halfedge(out_he_it)<<": "<<minimal_cone<<std::endl);

                /*auto min_expandability = minimal_cone.is_expandable(center_vertex, dummy_pos);
                if(min_expandability == -1){
                    PRINT_IF_NOT_SILENT(" error while checking expandability of minimal cone"<<std::endl);
                    return -1;
                }*/

                if(minimal_cone.n_cells() > best_minimal_cone.n_cells()){

                    best_minimal_cone = minimal_cone;
                    //PRINT_IF_NOT_SILENT(" --> updated best minimal cone as "<<best_minimal_cone<<std::endl);
                }
            }
        }


        if(!minimal_cone.n_cells()){
            //PRINT_IF_NOT_SILENT(" failure - Found no spoke with at least one expanded tet around it"<<std::endl);
            return 1;
        }

        minimal_cone = best_minimal_cone;

        return 0;
    }




    int Expander::extract_expansion_cone_around_spoke_edge(const HalfEdgeHandle& spoke_edge,
                                                           ExpansionCone& cone){

        //PRINT_IF_NOT_SILENT(" ------"<<std::endl);
        //PRINT_IF_NOT_SILENT(" extracting expansion cone around spoke edge "<<codomain_mesh_.find_halfedge(spoke_edge)<<std::endl);
        cone.clear();

        //to make the center vertex the first of the cone
        auto center_vertex = mesh_.from_vertex_handle(spoke_edge);
        cone.add_vertex(center_vertex, this->vertex(center_vertex));

        for(auto hec_it = mesh_.hec_iter(spoke_edge); hec_it.valid(); hec_it++){
            //PRINT_IF_NOT_SILENT(" - checking cell "<<*hec_it<<": "<<codomain_mesh_.get_cell_vertices(*hec_it)<<std::endl);
            if(expanded_vertices_count(*hec_it) == 3){
                //PRINT_IF_NOT_SILENT(" --> which contains 3 expanded vertices, adding it to the cone"<<std::endl);
                add_cell_to_cone(*hec_it, cone);
            }
        }
        return 0;
    }





    int Expander::grow_maximal_expansion_cone(const VertexHandle& center_vertex,
                                              const ExpansionCone& minimal_cone,
                                              ExpansionCone& cone){



        /*PRINT_IF_NOT_SILENT(" -------- "<<std::endl);
        PRINT_IF_NOT_SILENT(" growing expansion cone of "<<minimal_cone.n_cells()<<
                   " cells around vertex "<<center_vertex<<std::endl);
                   */

        cone = minimal_cone;

        //PRINT_IF_NOT_SILENT(" - adding minimal cone cells to visitation list"<<std::endl);
        //mark the neighboring tet already in the cone as such
        std::queue<CellHandle> to_visit;
        auto visited_prop = mesh_.request_cell_property<bool>();
        for(auto c: minimal_cone.cells()){
            auto mesh_cell = minimal_cone.cone_to_mesh_handle(c);
            if(c.idx() == -1){
                PRINT_IF_NOT_SILENT(" ERROR - cone cell "<<c<<" doesn't match any mesh cell"<<std::endl);
                return -1;
            }
            visited_prop[mesh_cell] = true;
            to_visit.push(mesh_cell);
            //PRINT_IF_NOT_SILENT(" -- added cell "<<mesh_cell<<std::endl);
        }

#warning temp check
        if(to_visit.size() != minimal_cone.n_cells()){
            PRINT_IF_NOT_SILENT(" ERROR - not as many cells in the cone as in the visitation list"<<std::endl);
            return -1;
        }

        VertexPosition new_pos;
        //add the neighboring tets not in the cone yet, growing from the minimal cone
        while(!to_visit.empty()){
            auto current_cell = to_visit.front();
            to_visit.pop();
            //PRINT_IF_NOT_SILENT(" -- current cell: "<<current_cell<<": "<<codomain_mesh_.get_cell_vertices(current_cell)<<std::endl);
            visited_prop[current_cell] = true;

            auto current_cone = cone;
            //PRINT_IF_NOT_SILENT(" -- current cone: "<<current_cone<<std::endl);

            for(auto chf_it = mesh_.chf_iter(current_cell); chf_it.valid(); chf_it++){
                //PRINT_IF_NOT_SILENT(" --- checking hf "<<*chf_it<<": "<<codomain_mesh_.get_halfface_vertices(*chf_it)<<std::endl);
                auto opposite_hf = mesh_.opposite_halfface_handle(*chf_it);
                if(opposite_hf.idx() != -1){
                    auto neighbor_cell = mesh_.incident_cell(opposite_hf);

                    if(neighbor_cell.idx() != -1 &&
                       !visited_prop[neighbor_cell] &&
                       expanded_vertices_count(neighbor_cell) == 3){
                        //PRINT_IF_NOT_SILENT(" ---- adding unvisited neighbor cell "<<neighbor_cell<<": "<<codomain_mesh_.get_cell_vertices(neighbor_cell)<<std::endl);
                        auto cone_cell = add_cell_to_cone(neighbor_cell, current_cone);
                        if(cone_cell.idx() == -1){
                            PRINT_IF_NOT_SILENT(" ERROR - couldn't add cell "<<neighbor_cell<<" to cone"<<std::endl);
                            return -1;
                        }



                        add_side_elements_to_expansion_cone(center_vertex,
                                                            current_cone);

                        //PRINT_IF_NOT_SILENT(" added side-elements, checking expandability..."<<std::endl);
                        auto result = current_cone.is_expandable(new_pos);

                        if(!result){
                            /*PRINT_IF_NOT_SILENT(" -----> cone is still expandable,"<<
                                       " updating maximal cone and adding cell "<<neighbor_cell<<
                                       " to visitation list"<<std::endl);
                                       */
                            cone = current_cone;
                            to_visit.push(neighbor_cell);
                        }else if(result == EXPANDABILITY_CHECK_ERROR){
                            /*PRINT_IF_NOT_SILENT(" error while adding cell "<<neighbor_cell<<
                                       ": "<<codomain_mesh_.get_cell_vertices(neighbor_cell)<<
                                       " to cone"<<std::endl);*/
                            return -1;
                        }else{
                            //PRINT_IF_NOT_SILENT(" -----> cone is no longer expandable"<<std::endl);
                        }
                    }
                }
            }
        }


        //PRINT_IF_NOT_SILENT(" ...done, maximal cone contains "<<cone.n_cells()<<" cells"<<std::endl);
        //PRINT_IF_NOT_SILENT(" -------- "<<std::endl);


        return 0;
    }


    int Expander::add_side_elements_to_expansion_cone(const VertexHandle& center_vertex,
                                                      ExpansionCone& cone){

        //PRINT_IF_NOT_SILENT(" ------ "<<std::endl);
        //PRINT_IF_NOT_SILENT(" adding side-elements to expansion cone..."<<std::endl);

        //check cells
        for(auto vc_it = mesh_.vc_iter(center_vertex); vc_it.valid(); vc_it++){
            auto mesh_vertices = mesh_.get_cell_vertices(*vc_it);
            //PRINT_IF_NOT_SILENT(" - checking cell "<<*vc_it<<": "<<mesh_vertices<<std::endl);
            std::vector<VertexHandle> cone_vertices;
            for(auto v: mesh_vertices){
                auto cone_vertex = cone.mesh_to_cone_handle(v);
                //PRINT_IF_NOT_SILENT(" -- "<<v<<" -> "<<cone_vertex<<std::endl);
                if(cone_vertex.idx() != -1){
                    cone_vertices.push_back(cone_vertex);
                }
            }

            if(cone_vertices.size() == 4 &&
               TopoHelper::cell_exists(cone, cone_vertices).idx() == -1
                /*cone.halfface(cone_vertices).idx() == -1*/){
                //PRINT_IF_NOT_SILENT(" --> all four vertices are in the cone but not the tet, adding it"<<std::endl);
                cone.add_cell(*vc_it, mesh_vertices);
            }
        }


        //and then check faces
        for(auto vf_it = mesh_.vf_iter(center_vertex); vf_it.valid(); vf_it++){
            auto mesh_vertices = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it,0));
            //PRINT_IF_NOT_SILENT(" - checking face "<<*vf_it<<": "<<mesh_vertices<<std::endl);
            std::vector<VertexHandle> cone_vertices;
            for(auto v: mesh_vertices){
                auto cone_vertex = cone.mesh_to_cone_handle(v);
                //PRINT_IF_NOT_SILENT(" -- "<<v<<" -> "<<cone_vertex<<std::endl);
                if(cone_vertex.idx() != -1){
                    cone_vertices.push_back(cone_vertex);
                }
            }

            if(cone_vertices.size() == 3 &&
               cone.find_halfface(cone_vertices).idx() == -1){
                //PRINT_IF_NOT_SILENT(" --> all three vertices are in the cone but not the face, adding it"<<std::endl);
                cone.add_face(mesh_vertices);
            }
        }

        //PRINT_IF_NOT_SILENT(" ... done"<<std::endl);
        //PRINT_IF_NOT_SILENT(" ------ "<<std::endl);

        return 0;
    }


    int Expander::perform_expansion_with_max_cone(const VertexHandle& center_vertex,
                                                  const ExpansionCone& maximal_cone){

        //PRINT_IF_NOT_SILENT(" --> re-performing max expansion, initial to-expand count = "<<to_expand_count_<<std::endl);
        size_t to_expand_count_before_spoke_splitting = to_expand_count_;

        //split-and-shrink the spokes not included in the expansion cone
        std::vector<VertexHandle> steiner_vertices;
        auto spoke_splitting_result = split_expanded_spokes_out_of_expansion_cone(center_vertex,
                                                                                  maximal_cone,
                                                                                  steiner_vertices);
        size_t to_expand_count_diff = to_expand_count_ - to_expand_count_before_spoke_splitting;

        //PRINT_IF_NOT_SILENT(" --> to-expand count after splitting spokes = "<<to_expand_count_<<std::endl);
        if(spoke_splitting_result == -1){
            PRINT_IF_NOT_SILENT(" error while splitting spokes out of expansion cone"<<std::endl);
            return -1;
        }

        if(!all_unexpanded_vertices_are_at_the_center(center_vertex)){
            PRINT_IF_NOT_SILENT(" ==> found bad vertices after splitting spokes"<<std::endl);
            return -1;
        }

        if(ExactBadTetFinder::meshContainsFlippedTets(mesh_,
                                                      vertex_position_prop_)){
            PRINT_IF_NOT_SILENT(" ==> found flipped tets after splitting spokes"<<std::endl);
            return -1;
        }

        auto max_cone_copy = maximal_cone;

        //and finally re-do the expanse with the maximal cone
        VertexPosition new_pos;
        auto expansion_result = max_cone_copy.is_expandable(new_pos);
        if(expansion_result){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't redo the maximal expansion. error code = "<<expansion_result<<std::endl);
            return -1;
        }


        finalize_expanse(center_vertex,
                         maximal_cone,
                         new_pos,
                         to_expand_count_diff);

        if(!all_unexpanded_vertices_are_at_the_center(center_vertex)){
            PRINT_IF_NOT_SILENT(" ==> found bad vertices after expanding vertex using maximal cone"<<std::endl);
            return -1;
        }

        if(ExactBadTetFinder::meshContainsFlippedTets(mesh_,
                                                      vertex_position_prop_)){
            PRINT_IF_NOT_SILENT(" ==> found flipped tets after expanding vertex using maximal cone"<<std::endl);
            return -1;
        }


        if(steiner_vertices.size() == 1){
            PRINT_IF_NOT_SILENT(" --> single SV expansion, trying to expand it right-away"<<std::endl);

            auto sv = steiner_vertices[0];

            ExpansionCone cone;
            auto result = gather_full_expansion_cone(sv, cone);

            if(result == -1){
                PRINT_IF_NOT_SILENT(" error while gathering expansion cone for single SV"<<std::endl);
                return -1;
            }

            VertexPosition new_pos;
            result = cone.is_expandable(new_pos);


            if(!result) {
                finalize_expanse(sv,
                                 cone,
                                 new_pos);
                PRINT_IF_NOT_SILENT(" ===> SV could succesfully be expanded!"<<std::endl);
                PRINT_IF_NOT_SILENT(" ===> expanded secondary cone "<<cone<<std::endl);
                single_SV_topo_expandable_and_expandable_SV_count_++;
            }else{
                PRINT_IF_NOT_SILENT(" ===> failed to expand SV..."<<std::endl);
                PRINT_IF_NOT_SILENT(" original expanse : "<<maximal_cone<<std::endl);
                PRINT_IF_NOT_SILENT(" SV cone : "<<cone<<std::endl);
                single_SV_topo_expandable_but_not_expandable_SV_count_++;
            }
        }


        return 0;
    }





    int Expander::apply_split_list(const SplitList& split_list,
                                   VertexPropertyT<VertexHandle>& submesh_to_mesh_prop,
                                   bool star_shapifying_cluster){

        PRINT_IF_NOT_SILENT(" =============================================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" =============== APPLYING SPLIT LIST TO MESH..."<<std::endl);

        PRINT_IF_NOT_SILENT(" submesh -> mesh prop size = "<<submesh_to_mesh_prop.size()<<std::endl);
        PRINT_IF_NOT_SILENT(" split list size = "<<split_list.size()<<std::endl);

        int initial_split_list_size(split_list_.size());

        int split_count(0);
        //PRINT_IF_NOT_SILENT(" applying split list to original mesh"<<std::endl);
        for(auto& split: split_list){
            check_for_timeout();
            //PRINT_IF_NOT_SILENT(" ---------- applying split "<<split<<std::endl);
            auto mesh_from_vertex = submesh_to_mesh_prop[split.cone_from_vertex];
            if(mesh_from_vertex.idx() == -1){
                PRINT_IF_NOT_SILENT(" -- ERROR - couldn't find from vertex in submesh->mesh map"<<std::endl);
                return -1;
            }

            auto mesh_to_vertex = submesh_to_mesh_prop[split.cone_to_vertex];
            if(mesh_to_vertex.idx() == -1){
                PRINT_IF_NOT_SILENT(" -- ERROR - couldn't find to vertex in submesh->mesh map"<<std::endl);
                return -1;
            }

            auto mesh_heh = mesh_.find_halfedge(mesh_from_vertex,
                                           mesh_to_vertex);

            if(mesh_heh.idx() == -1){
                PRINT_IF_NOT_SILENT(" ERROR - couldn't find halfedge ("<<mesh_from_vertex<<", "<<mesh_to_vertex<<") in mesh"<<std::endl);
                PRINT_IF_NOT_SILENT(" ---------------------- neighbors of "<<mesh_from_vertex<<": "<<std::endl);
                for(auto out_he: mesh_.outgoing_halfedges(mesh_from_vertex)){
                    PRINT_IF_NOT_SILENT(" -- "<<mesh_.to_vertex_handle(out_he)<<std::endl);
                }
                PRINT_IF_NOT_SILENT(" ---------------------- neighbors of "<<mesh_to_vertex<<": "<<std::endl);
                for(auto out_he: mesh_.outgoing_halfedges(mesh_to_vertex)){
                    PRINT_IF_NOT_SILENT(" -- "<<mesh_.to_vertex_handle(out_he)<<std::endl);
                }
                return -1;
            }
            //PRINT_IF_NOT_SILENT(" - halfedge: "<<mesh_heh<<codomain_mesh_.find_halfedge(mesh_heh)<<std::endl);

            auto mesh_new_vertex = split_edge(mesh_heh);
            //PRINT_IF_NOT_SILENT(" - new mesh vertex: "<<mesh_new_vertex<<" moved to "<<vec2vec(split.new_vertex_position)<<std::endl);
            if(split.cone_new_vertex.idx() >= (int)submesh_to_mesh_prop.size()){
                PRINT_IF_NOT_SILENT(" ERROR - submesh to mesh prop has size "<<submesh_to_mesh_prop.size()<<" but trying to set vertex "<<split.cone_new_vertex<<std::endl);
                return -1;
            }
            submesh_to_mesh_prop[split.cone_new_vertex] = mesh_new_vertex;



            //PRINT_IF_NOT_SILENT(" - inserted ("<<split.cone_new_vertex<<"->"<<mesh_new_vertex<<") in map"<<std::endl);

            //move_vertex_and_mark_as_expanded(mesh_new_vertex, split.new_vertex_position);
            //expanded_prop_[mesh_new_vertex] = true;
            set_expanded_prop(mesh_new_vertex, true);
            this->set_vertex(mesh_new_vertex, split.new_vertex_position);

            //PRINT_IF_NOT_SILENT(" - marked vertex "<<mesh_new_vertex<<" as expanded "<<std::endl);

            //appending the split to the list for nested expansion
            split_list_.push_back({mesh_from_vertex,
                                   mesh_to_vertex,
                                   mesh_new_vertex,
                                   split.new_vertex_position});

            /*PRINT_IF_NOT_SILENT(" - forwarded submesh split n°"<<split_count<<
                                " "<<split<<
                                " into split "<<split_list_.back()<<
                                ", position size = "<<byte_size(split_list_.back().new_vertex_position)<<std::endl);*/


            if(split.cone_collapsed_to_vertex.is_valid()){
                //PRINT_IF_NOT_SILENT(" -> split "<<split<<" was followed by a collapse, forwarding it to mesh"<<std::endl);
                auto mesh_collapsed_to_vertex = submesh_to_mesh_prop[split.cone_collapsed_to_vertex];
                auto mesh_to_collapse_he = mesh_.find_halfedge(mesh_new_vertex,
                                                          mesh_collapsed_to_vertex);
                if(!mesh_to_collapse_he.is_valid()){
                    PRINT_IF_NOT_SILENT(" ERROR - couldn't recover mesh halfedge ("<<mesh_collapsed_to_vertex<<", "<<mesh_new_vertex<<")"<<std::endl);
                    return -1;
                }

                //std::cout<<" --> link condition: "<<link_condition(codomain_mesh_, mesh_to_collapse_he)<<std::endl);


                //PRINT_IF_NOT_SILENT(" - collapsing edge "<<codomain_mesh_.find_halfedge(mesh_to_collapse_he)<<" of valence "<<codomain_mesh_.valence(codomain_mesh_.edge_handle(mesh_to_collapse_he))<<std::endl);

                mesh_.collapse_edge(mesh_to_collapse_he);
                split_list_.back().cone_collapsed_to_vertex = mesh_collapsed_to_vertex;
            }else{

                mark_halfedges_as_new_and_collapsible(mesh_from_vertex,
                                                      mesh_new_vertex,
                                                      mesh_to_vertex);

            }

            //temp check
            /*if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_)){
                PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after applying split"<<std::endl);
                auto bad_tets = ExactBadTetFinder::findBadTets(codomain_mesh_, vertex_position_prop_);
                for(auto flipped_tet: bad_tets.second){
                    PRINT_IF_NOT_SILENT(" - "<<flipped_tet<<": "<<codomain_mesh_.get_cell_vertices(flipped_tet)<<std::endl);
                }
                return -1;
            }*/
            //temp check ends here

            split_count++;
        }


#if ENABLE_ALL_CHECKS
        if(!star_shapifying_cluster){
            std::vector<CellHandle> bad_cells;
            for(int i(initial_split_list_size); i<(int)split_list_.size(); i++){
                auto mid_v = split_list_[i].cone_new_vertex;
                if(!mesh_.is_deleted(mid_v)){
                    for(auto vc_it = mesh_.vc_iter(mid_v); vc_it.valid(); vc_it++){
                        auto c_vertices = mesh_.get_cell_vertices(*vc_it);
                        auto unexp_count(0);
                        for(auto cv: c_vertices){
                            unexp_count += !expanded_prop_[cv];
                        }
                        if(unexp_count <= 1){
                            if(OVMtetToCGALtet(mesh_, vertex_position_prop_, *vc_it).is_degenerate()){
                                bad_cells.push_back(*vc_it);
                            }
                        }
                    }
                }
            }


            if(!bad_cells.empty()){

                PRINT_IF_NOT_SILENT(" ERROR/WARNING: 0-volume cells with only one (or less) unexpanded vertex: "<<std::endl);

                if(bad_cells.size() < 50){
                    for(auto c: bad_cells){
                        auto c_vertices = mesh_.get_cell_vertices(c);
                        auto volume = OVMtetToCGALtet(mesh_, vertex_position_prop_, c).volume();
                        PRINT_IF_NOT_SILENT(" -- "<<c<<": "<<c_vertices<<", volume = "<<volume<<std::endl);
                    }
                }else{
                    PRINT_IF_NOT_SILENT(" --> too many bad cells, skipping output"<<std::endl);
                }

                //PRINT_IF_NOT_SILENT(" --> remaining unexpanded vertices count = "<<left_to_expand_count()<<std::endl);

                return -1;
            }
        }
#endif

        //temp check
        /*if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_)){
            std::cout<<" ERROR - mesh contains flipped tets right after split list application"<<std::endl);
            return -1;
        }*/

        /*PRINT_IF_NOT_SILENT(" -------- cone->mesh prop after applying split list: "<<std::endl);
        for(int i(0); i<(int)submesh_to_mesh_prop.size(); i++){
            VertexHandle v(i);
            PRINT_IF_NOT_SILENT(" - "<<v<<" -> "<<submesh_to_mesh_prop[v]<<std::endl);
        }*/
        PRINT_IF_NOT_SILENT(" ...done with split list application"<<std::endl);
        PRINT_IF_NOT_SILENT(" =============================================================="<<std::endl);

        return 0;
    }






    void Expander::mark_halfedges_as_new_and_collapsible(const VertexHandle& from_vertex,
                                                         const VertexHandle& mid_vertex,
                                                         const VertexHandle& to_vertex){
        auto first_new_he = mesh_.find_halfedge(mid_vertex, from_vertex);
        if(!first_new_he.is_valid()){
            std::cout<<" ERROR - couldn't recover new halfedge "<<mid_vertex<<" - "<<from_vertex<<std::endl;
            return;
        }
        collapsible_new_halfedge_prop_[first_new_he] = true;

        auto second_new_he = mesh_.find_halfedge(mid_vertex, to_vertex);
        if(!first_new_he.is_valid()){
            std::cout<<" ERROR - couldn't recover new halfedge "<<mid_vertex<<" - "<<to_vertex<<std::endl;
            return;
        }
        collapsible_new_halfedge_prop_[second_new_he] = true;

        //PRINT_IF_NOT_SILENT(" - set edges "<<first_new_he<<": "<<codomain_mesh_.find_halfedge(first_new_he)<<" and "<<second_new_he<<": "<<codomain_mesh_.find_halfedge(second_new_he)<<" as collapsible new edges"<<std::endl);
    }



    void Expander::reduce_mid_vertices_precision(int split_list_start_index,
                                                 int position_byte_size_threshold){

        if(split_list_start_index >= (int)split_list_.size()){
            return;
        }

        //0 -> binary search
        //1 -> chebyshev centroid
        const int mode(0);

        PRINT_IF_NOT_SILENT(" =========================================================================================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" =============== MOVING SPLIT LIST MID-VERTICES WITH MODE: "<<(mode ? " BINARY-SEARCH": " 1-RING NEIGHBORHOOD'S CHEBYSHEV CENTER")<<"..."<<std::endl);

        int max, total;
        double avg;
        VertexHandle max_precision_vh;
        compute_precision_stats(mesh_, vertex_position_prop_, max, total, avg, max_precision_vh);
        PRINT_IF_NOT_SILENT(" before: max = "<<max<<", avg = "<<avg<<", total = "<<total<<" (max precision vertex: "<<max_precision_vh<<")"<<std::endl);
        int initial_total = total;

        auto start_time = std::chrono::high_resolution_clock::now();
        int reduction_count(0);
        for(int i(split_list_start_index); i<(int)split_list_.size(); i++){

            check_for_timeout();

            auto mid_v = split_list_[i].cone_new_vertex;

            if(mesh_.is_boundary(mid_v) || mesh_.is_deleted(mid_v)){
                //PRINT_IF_NOT_SILENT(" -> mid-vertex "<<mid_v<<" is boundary"<<std::endl);
                continue;
            }

            if(byte_size(vertex(mid_v)) < position_byte_size_threshold){
                continue;
            }

            //PRINT_IF_NOT_SILENT(" -> checking mid-vertex "<<mid_v<<std::endl);

            ExpansionCone cone;
            ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(mesh_,
                                                                        vertex_position_prop_,
                                                                        mid_v,
                                                                        cone);

            VertexPosition new_pos;
            if(mode){
                cone.is_geo_expandable(new_pos);
            }else{

                new_pos = binary_search_minimum_precision_position(mid_v,
                                                                   cone,
                                                                   this->vertex(mid_v),
                                                                   position_byte_size_threshold);
            }

            if(byte_size(new_pos) < byte_size(split_list_[i].new_vertex_position)){
                //PRINT_IF_NOT_SILENT(" --> reduced position size from "<<byte_size(split_list_[i].new_vertex_position)<<" to "<<byte_size(new_pos)<<" for vertex "<<mid_v<<std::endl);
                split_list_[i].new_vertex_position = new_pos;
                this->set_vertex(mid_v, new_pos);
                reduction_count++;
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        float duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000000;
        compute_precision_stats(mesh_, vertex_position_prop_, max, total, avg, max_precision_vh);
        PRINT_IF_NOT_SILENT(" after "<<reduction_count<<" reductions in "<<duration_s<<"s "
                   ": max = "<<max<<
                   ", avg = "<<avg<<
                   ", total = "<<total<<
                   " ==> total saved = "<<(initial_total - total)<<
                   " (max precision vertex: "<<max_precision_vh<<")"<<std::endl);

        total_saved_position_bytes_ += (initial_total - total);
        PRINT_IF_NOT_SILENT(" =========================================================================================================="<<std::endl);
    }



    void Expander::smooth_neighborhood_of_unexpanded_vertices(int ring_k,
                                                              int max_iterations,
                                                              double delta_epsilon,
                                                              int position_byte_size_threshold){

        PRINT_IF_NOT_SILENT(" - smoothing "<<ring_k<<"-ring expanded neighborhood of unexpanded vertices"<<std::endl);
        auto added_prop = mesh_.request_vertex_property<bool>();
        std::vector<VertexHandle> current_ring;
        for(auto v: mesh_.vertices()){
            if(!expanded_prop_[v]){
                current_ring.push_back(v);
                added_prop[v] = true;
            }
        }

        std::vector<VertexHandle> k_ring_neighborhood;
        std::vector<VertexHandle> next_ring;
        for(int i(0); i<ring_k; i++){
            next_ring.clear();

            for(auto v: current_ring){
                for(auto vv_it = mesh_.vv_iter(v); vv_it.valid(); vv_it++){
                    if(!added_prop[*vv_it] &&
                            expanded_prop_[*vv_it] &&
                            !mesh_.is_boundary(*vv_it)){

                        added_prop[*vv_it] = true;
                        next_ring.push_back(*vv_it);
                        k_ring_neighborhood.push_back(*vv_it);
                    }
                }
            }

            current_ring = next_ring;
        }

        //PRINT_IF_NOT_SILENT(" vertices to smooth: "<<k_ring_neighborhood<<std::endl);

        smooth_vertices(k_ring_neighborhood,
                        max_iterations,
                        delta_epsilon,
                        position_byte_size_threshold);
    }



    void Expander::smooth_mid_vertices_1_ring_neighborhood(int split_list_start_index,
                                                           int max_iterations,
                                                           double delta_epsilon,
                                                           int position_byte_size_threshold){

        PRINT_IF_NOT_SILENT(" - smoothing 1-ring expanded neighborhood of mid-vertices"<<std::endl);
        auto added_prop = mesh_.request_vertex_property<bool>();
        std::vector<VertexHandle> to_smooth;
        for(int i(split_list_start_index); i<(int)split_list_.size(); i++){
            auto mid_v = split_list_[i].cone_new_vertex;

            if(!mesh_.is_deleted(mid_v) &&
                    expanded_prop_[mid_v] &&
                    !mesh_.is_boundary(mid_v)){

                //checking here so we still ciruculate through neighbors, even if it was already added
                if(!added_prop[mid_v]){
                    to_smooth.push_back(mid_v);
                    added_prop[mid_v] = true;
                }

                for(auto vv_it = mesh_.vv_iter(mid_v); vv_it.valid(); vv_it++){
                    if(expanded_prop_[*vv_it] && !mesh_.is_boundary(*vv_it) && !added_prop[*vv_it]){
                        to_smooth.push_back(*vv_it);
                        added_prop[*vv_it] = true;
                    }
                }
            }
        }

        //PRINT_IF_NOT_SILENT(" vertices to smooth: "<<k_ring_neighborhood<<std::endl);

        smooth_vertices(to_smooth,
                        max_iterations,
                        delta_epsilon,
                        position_byte_size_threshold);
    }

    void Expander::smooth_full_interior(int max_iterations,
                                        double delta_epsilon,
                                        int position_byte_size_threshold){

        PRINT_IF_NOT_SILENT(" smoothing full expanded interior"<<std::endl);
        std::vector<VertexHandle> interior_expanded_vertices;

        for(auto v: mesh_.vertices()){
            if(!mesh_.is_boundary(v) && expanded_prop_[v]){
                interior_expanded_vertices.push_back(v);
            }
        }

        smooth_vertices(interior_expanded_vertices,
                        max_iterations,
                        delta_epsilon,
                        position_byte_size_threshold);
    }


    void Expander::smooth_vertices(const std::vector<VertexHandle>& vertices_to_smooth,
                                   int max_iterations,
                                   double delta_epsilon,
                                   int position_byte_size_threshold){


        bool reduce_precision_after_each_smoothing_iterations(true);


        PRINT_IF_NOT_SILENT(" =========================================================================================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" =============== UNIFORMLY SMOOTHING EXPANDED INTERIOR VERTICES..."<<std::endl);

        int max, total;
        double avg;
        VertexHandle max_precision_vh;
        compute_precision_stats(mesh_, vertex_position_prop_, max, total, avg, max_precision_vh);
        PRINT_IF_NOT_SILENT(" before: max = "<<max<<", avg = "<<avg<<", total = "<<total<<" (max precision vertex: "<<max_precision_vh<<")"<<std::endl);
        int initial_total = total;

        //first pass for smoothing
        //int reduction_count(0);
        double max_delta(delta_epsilon + 1);
        auto next_pos_prop = mesh_.request_vertex_property<VertexPosition>();
        auto smoothing_start_time = std::chrono::high_resolution_clock::now();
        //first, check which tets are already degenerates
        auto deg_tet_prop = mesh_.request_cell_property<bool>();
        for(auto c: mesh_.cells()){
            deg_tet_prop[c] = OVMtetToCGALtet(mesh_, vertex_position_prop_, c).is_degenerate();
        }
        auto shortest_length_prop = mesh_.request_vertex_property<double>("");

        int iteration_count(0);
        while(iteration_count < max_iterations && max_delta > delta_epsilon){
            int smoothed_vertices_count(0);
            int interior_vertices_count(0);
            max_delta = 0;
            PRINT_IF_NOT_SILENT(" --------- running smoothing iteration "<<(iteration_count+1)<<"..."<<std::endl);
            auto iteration_start_time = std::chrono::high_resolution_clock::now();
            int max_val(0);
            double avg_val(0);


            for(auto v: vertices_to_smooth){

                check_for_timeout();

                if(!iteration_count){
                    moved_during_smoothing_prop_[v] = false;
                }

                if(mesh_.is_boundary(v) || !expanded_prop_[v]){
                    continue;
                }

                next_pos_prop[v] = {0,0,0};

                //PRINT_IF_NOT_SILENT(" -> checking mid-vertex "<<mid_v<<std::endl);

                float n(0);
                bool found_origin(false);
                shortest_length_prop[v] = std::numeric_limits<double>::max();

                for(auto vv_it = mesh_.vv_iter(v); vv_it.valid(); vv_it++){
                    next_pos_prop[v] += vertex(*vv_it);
                    if(vertex(*vv_it) == VertexPosition(0,0,0)){
                        found_origin = true;
                    }else{
                        n++;
                    }
                    shortest_length_prop[v] = std::min(shortest_length_prop[v],
                                                       CGAL::to_double(linf_norm(vertex(*vv_it) - vertex(v))));
                }
                next_pos_prop[v] /= (n + found_origin);
            }


            for(auto v: vertices_to_smooth){

                check_for_timeout();

                if(mesh_.is_boundary(v) || !expanded_prop_[v]){
                    //PRINT_IF_NOT_SILENT(" -> mid-vertex "<<mid_v<<" is boundary"<<std::endl);
                    continue;
                }

                interior_vertices_count++;

                auto current_pos = vertex(v);
                set_vertex(v, next_pos_prop[v]);
                if(ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(mesh_, vertex_position_prop_, v)){
                    set_vertex(v, current_pos);
                    continue;
                }
                for(auto vc_it = mesh_.vc_iter(v); vc_it.valid(); vc_it++){
                    bool is_deg = OVMtetToCGALtet(mesh_, vertex_position_prop_, *vc_it).is_degenerate();
                    if(deg_tet_prop[*vc_it] != is_deg){
                        set_vertex(v, current_pos);
                        continue;
                    }
                }

                if(reduce_precision_after_each_smoothing_iterations){
                    const auto initial_byte_size = byte_size(vertex(v));
                    //PRINT_IF_NOT_SILENT(" - vertex "<<v<<" position byte size = "<<initial_byte_size<<std::endl);
                    if(initial_byte_size < position_byte_size_threshold){
                        continue;
                    }

                    ExpansionCone cone;
                    ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(mesh_,
                                                                                vertex_position_prop_,
                                                                                v,
                                                                                cone);

                    VertexPosition new_pos = binary_search_minimum_precision_position(v,
                                                                                      cone,
                                                                                      this->vertex(v),
                                                                                      position_byte_size_threshold);


                    if(byte_size(new_pos) < initial_byte_size){
                        //PRINT_IF_NOT_SILENT(" --> reduced position size from "<<initial_byte_size<<" to "<<byte_size(new_pos)<<" for vertex "<<v<<std::endl);
                        this->set_vertex(v, new_pos);
                        moved_during_smoothing_prop_[v] = true;
                    }
                }

                if(shortest_length_prop[v] > DBL_EPSILON){
                    max_delta = std::max(max_delta, CGAL::to_double(linf_norm(current_pos - vertex(v))) / shortest_length_prop[v]);
                    //PRINT_IF_NOT_SILENT(" shortest length = "<<shortest_length_prop[v]<<std::endl);
                }
                //std::cout<<" -- smoothed vertex "<<v<<std::endl);
                smoothed_vertices_count++;
                moved_during_smoothing_prop_[v] = true;
                max_val = std::max(max_val, (int)mesh_.valence(v));
                avg_val += mesh_.valence(v);

            }
            avg_val /= smoothed_vertices_count;
            auto iteration_end_time = std::chrono::high_resolution_clock::now();
            float iteration_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(iteration_end_time - iteration_start_time).count() / 1000000;
            PRINT_IF_NOT_SILENT(" --> smoothed "<<smoothed_vertices_count<<"/"<<interior_vertices_count<<" expanded interior vertices in "<<iteration_duration_s<<"s"<<
                       ", max delta = "<<max_delta<<std::endl);
            PRINT_IF_NOT_SILENT(" --> max valence = "<<max_val<<", average valence = "<<avg_val<<std::endl);

            if(smoothed_vertices_count){
                data_logger_.add_smoothing_computation_point(iteration_duration_s * 1e6, max_val, avg_val * smoothed_vertices_count, smoothed_vertices_count);
            }

            iteration_count++;
        }
        auto smoothing_end_time = std::chrono::high_resolution_clock::now();
        float smoothing_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(smoothing_end_time - smoothing_start_time).count() / 1000000;
        compute_precision_stats(mesh_, vertex_position_prop_, max, total, avg, max_precision_vh);
        int smoothing_total = total;
        PRINT_IF_NOT_SILENT(" after "<<iteration_count<<" smoothing iterations in "<<smoothing_duration_s<<"s "
                   ": max = "<<max<<
                   ", avg = "<<avg<<
                   ", total = "<<total<<
                   " ==> total increase = "<<(smoothing_total - initial_total)<<
                   " (max precision vertex: "<<max_precision_vh<<")"<<std::endl);

        //second pass for precision reduction
        auto start_time = std::chrono::high_resolution_clock::now();
        int reduction_count(0);
        for(auto v: vertices_to_smooth){

            check_for_timeout();
            if(mesh_.is_boundary(v) || !expanded_prop_[v]){
                //PRINT_IF_NOT_SILENT(" -> vertex "<<v<<" is boundary or expanded"<<std::endl);
                continue;
            }

            const auto initial_byte_size = byte_size(vertex(v));
            //PRINT_IF_NOT_SILENT(" - vertex "<<v<<" position byte size = "<<initial_byte_size<<std::endl);
            if(initial_byte_size < position_byte_size_threshold){
                continue;
            }

            ExpansionCone cone;
            ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(mesh_,
                                                                        vertex_position_prop_,
                                                                        v,
                                                                        cone);

            VertexPosition new_pos = binary_search_minimum_precision_position(v,
                                                                              cone,
                                                                              this->vertex(v),
                                                                              position_byte_size_threshold);


            if(byte_size(new_pos) < initial_byte_size){
                //PRINT_IF_NOT_SILENT(" --> reduced position size from "<<initial_byte_size<<" to "<<byte_size(new_pos)<<" for vertex "<<v<<std::endl);
                this->set_vertex(v, new_pos);
                reduction_count++;
                moved_during_smoothing_prop_[v] = true;
            }else{
                //PRINT_IF_NOT_SILENT(" --> position size went from "<<initial_byte_size<<" to "<<byte_size(new_pos)<<", skipping"<<std::endl);
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        float duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000000;


        compute_precision_stats(mesh_, vertex_position_prop_, max, total, avg, max_precision_vh);
        PRINT_IF_NOT_SILENT(" after "<<reduction_count<<" reductions in "<<duration_s<<"s "
                   ": max = "<<max<<
                   ", avg = "<<avg<<
                   ", total = "<<total<<
                   " ==> total variation = "<<(smoothing_total - initial_total)<<" "<<(total - smoothing_total)<<" = "<<(total - initial_total)<<
                   " (max precision vertex: "<<max_precision_vh<<")"<<std::endl);

        total_saved_position_bytes_ += (initial_total - total);
        PRINT_IF_NOT_SILENT(" =========================================================================================================="<<std::endl);

        reset_EC_memory();

        total_smoothing_time_s_ += duration_s;
    }


#if 0
    int Expander::collapse_new_edges(const int split_list_start_index){

        if(split_list_start_index >= (int)split_list_.size()){
            return 0;
        }

        PRINT_IF_NOT_SILENT(" - collapsing all new edges starting from index "<<split_list_start_index<<std::endl);
        std::vector<std::pair<VertexHandle, VertexHandle>> to_collapse;

        for(int i(split_list_start_index); i<(int)split_list_.size(); i++){

            auto const mid_v  = split_list_[i].cone_new_vertex;
            auto const from_v = split_list_[i].cone_from_vertex;
            auto const to_v   = split_list_[i].cone_to_vertex;

            if(codomain_mesh_.is_boundary(mid_v) ||
                    codomain_mesh_.is_deleted(mid_v) ||
                    codomain_mesh_.is_deleted(from_v) ||
                    codomain_mesh_.is_deleted(to_v)){
                continue;
            }
            to_collapse.push_back({mid_v, from_v});
            to_collapse.push_back({mid_v, to_v});
        }

        return collapse_as_many_edges_as_possible(to_collapse, 100);
    }
#endif


    int Expander::collapse_new_edges_opposite_to_unexpanded_vertices(){

        PRINT_IF_NOT_SILENT(" - collapsing all edges opposite to unexpanded vertices"<<std::endl);
        std::vector<std::pair<VertexHandle, VertexHandle>> to_collapse;


        for(auto v: mesh_.vertices()){
            if(expanded_prop_[v]){
                continue;
            }

            for(auto vf_it = mesh_.vf_iter(v); vf_it.valid(); vf_it++){
                auto hf_vertices = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it,0));
                HalfEdgeHandle op_he1(-1);
                if(hf_vertices[0] == v){
                    op_he1 = mesh_.find_halfedge(hf_vertices[1], hf_vertices[2]);
                }else if(hf_vertices[1] == v){
                    op_he1 = mesh_.find_halfedge(hf_vertices[0], hf_vertices[2]);
                }else{
                    op_he1 = mesh_.find_halfedge(hf_vertices[0], hf_vertices[1]);
                }

                if(mesh_.is_deleted(mesh_.halfedge(op_he1).from_vertex()) ||
                        mesh_.is_deleted(mesh_.halfedge(op_he1).to_vertex())){
                    continue;
                }

                if(collapsible_new_halfedge_prop_[op_he1]){
                    to_collapse.push_back({mesh_.from_vertex_handle(op_he1),
                                           mesh_.to_vertex_handle(op_he1)});
                    //std::cout<<" - added edge "<<op_he1<<": "<<codomain_mesh_.find_halfedge(op_he1)<<" for collapse"<<std::endl);
                    /*if(!is_candidate_prop[op_he1]){
                        std::cout<<"  --> he1 not an old candidate!"<<std::endl);
                        return -1;
                    }*/

                }else{
                    auto op_he2 = mesh_.opposite_halfedge_handle(op_he1);
                    if(collapsible_new_halfedge_prop_[op_he2]){
                        to_collapse.push_back({mesh_.from_vertex_handle(op_he2),
                                               mesh_.to_vertex_handle(op_he2)});
                        //std::cout<<" - added edge "<<op_he2<<": "<<codomain_mesh_.find_halfedge(op_he2)<<" for collapse"<<std::endl);
                        /*if(!is_candidate_prop[op_he2]){
                            std::cout<<"  --> he2 not an old candidate!"<<std::endl);
                            return -1;
                        }*/
                    }
                }
            }
        }



        return collapse_as_many_edges_as_possible(to_collapse);

    }




    int Expander::collapse_as_many_edges_as_possible(const std::vector<std::pair<VertexHandle, VertexHandle>>& to_collapse){

        if(to_collapse.empty()){
            return 0;
        }

        //0 -> Chebyshev centroid
        //1 -> mid-point
        int mode(EDGE_COLLAPSE_MODE);

        int max_collapsed_edges_count = EDGE_COLLAPSE_MAX_EDGES_TO_COLLAPSE_COUNT;


        if(max_collapsed_edges_count == -1){
            PRINT_IF_NOT_SILENT(" --> -1 max collapsed edges count, using candidates list size"<<std::endl);
            max_collapsed_edges_count = to_collapse.size();
        }


        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" TRYING TO COLLAPSE "<<max_collapsed_edges_count<<" EDGES AMONG "<<to_collapse.size()<<" CANDIDATES WITH MODE "<<mode<<"..."<<std::endl);

        TopoHelper topo_helper(mesh_);

        int iteration_count(0);
        int collapse_count(0);
        int total_collapse_count(0);

        const int max_iteration_count(EDGE_COLLAPSE_MAX_ITERATION_COUNT);


        auto start_time = std::chrono::high_resolution_clock::now();

        auto valence_at_last_iteration_prop = mesh_.request_vertex_property<size_t>("", 0);
        auto moved_at_last_iteration_prop = mesh_.request_vertex_property<bool>("", true);

        for(auto from_to_pair: to_collapse){
            if(valence_at_last_iteration_prop[from_to_pair.first]){
                PRINT_IF_NOT_SILENT(" - valence of vertex "<<from_to_pair.first<<" is "<<valence_at_last_iteration_prop[from_to_pair.first]<<std::endl);
            }
        }
        bool timedout(false);
        do{
            collapse_count = 0;
            int unchanged_valence_count(0);
            int unchanged_neighborhood_count(0);
            int invalid_edges_count(0);
            int uncollapsible_count(0);
            int unexp_one_ring_count(0);
            int tried_count(0);
            float average_steps_count(0);
            int centroid_used_count(0);
            //PRINT_IF_NOT_SILENT(" ============================"<<std::endl);
            //PRINT_IF_NOT_SILENT(" running iteration "<<(iteration_count + 1)<<std::endl);
            auto iteration_start_time = std::chrono::high_resolution_clock::now();

            auto already_tried_during_this_iteration_prop = mesh_.request_vertex_property<bool>();

            for(auto from_to_pair: to_collapse){


                auto const from_v = from_to_pair.first;
                auto const to_v   = from_to_pair.second;

                //std::cout<<" - trying pair "<<from_v<<" (val "<<codomain_mesh_.valence(from_v)<<"), "<<to_v<<std::endl);

                moved_at_last_iteration_prop[from_v] = false;

                if(mesh_.is_boundary(from_v) ||
                        mesh_.is_deleted(from_v) ||
                        mesh_.is_deleted(to_v) ||
                        !expanded_prop_[from_v] ||
                        !expanded_prop_[to_v]){
                    /*std::cout<<" - skipping edge ("<<from_v<<", "<<to_v<<"): "<<std::endl);
                    std::cout<<"         - from-v is boundary: "<<codomain_mesh_.is_boundary(from_v)<<std::endl);
                    std::cout<<"         - deleted: "<<codomain_mesh_.is_deleted(from_v)<<", "<<codomain_mesh_.is_deleted(to_v)<<std::endl);
                    std::cout<<"         - expanded: "<<expanded_prop_[from_v]<<", "<<expanded_prop_[to_v]<<std::endl);*/

                    continue;
                }

                tried_count++;

                if(mesh_.valence(from_v) == valence_at_last_iteration_prop[from_v] &&
                        !already_tried_during_this_iteration_prop[from_v]){
                   //PRINT_IF_NOT_SILENT(" --> same valence (="<<codomain_mesh_.valence(from_v)<<") as last iteration, skipping"<<std::endl);
                    unchanged_valence_count++;
                    continue;
                }

                bool one_neighbor_moved_at_last_iteration(false);
                for(auto vv_it = mesh_.vv_iter(from_v); vv_it.valid(); vv_it++){
                    if(moved_at_last_iteration_prop[*vv_it]){
                        one_neighbor_moved_at_last_iteration = true;
                        break;
                    }
                }
                if(!one_neighbor_moved_at_last_iteration){
                    //PRINT_IF_NOT_SILENT(" -- no neighbors moved at last iteration, skipping"<<std::endl);
                    unchanged_neighborhood_count++;
                    continue;
                }

                valence_at_last_iteration_prop[from_v] = mesh_.valence(from_v);
                already_tried_during_this_iteration_prop[from_v] = true;

                try{
                    check_for_timeout();
                }catch(TimeOutException e){
                    PRINT_IF_NOT_SILENT(" --> global timeout"<<std::endl);
                    timedout = true;
                    break;
                }

                float current_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count() / 1000000;
                if(current_duration_s > EDGE_COLLAPSE_ITERATION_TIMEOUT_S){
                    PRINT_IF_NOT_SILENT(" --> local timeout"<<std::endl);
                    timedout = true;
                    break;
                }


                //PRINT_IF_NOT_SILENT(" - checking split "<<split_list_[i]<<std::endl);

                VertexHandle surviving_vertex(-1);

                auto he = mesh_.find_halfedge(from_v, to_v);
                bool is_collapsible = false;
                if(he.is_valid()){
                    is_collapsible = topo_helper.isCollapsible(mesh_.edge_handle(he));
                    uncollapsible_count += !is_collapsible;
                }
                invalid_edges_count += !he.is_valid();

                //Note for he.is_valid(): the split edges might not exist anymore
                //e.g. split edge (a-b) -> (a-c-b) and then -> (a-c-d-b),
                //then edge (c-b) doesn't exist anymore
                if(he.is_valid() && is_collapsible){
                    //PRINT_IF_NOT_SILENT(" -- edge "<<codomain_mesh_.find_halfedge(he)<<" is collapsible, checking star-shape"<<std::endl);

                    bool valid_in_domain(true);
                    //check that this is a valid collapse in the domain mesh
                    auto initial_domain_pos = domain_vertex_position_prop_[from_v];
                    domain_vertex_position_prop_[from_v] = domain_vertex_position_prop_[to_v];
                    for(auto vc_it = mesh_.vc_iter(from_v); vc_it.valid(); vc_it++){
                        bool found_to_vertex(false);
                        for(auto v: mesh_.get_cell_vertices(*vc_it)){
                            if(v == to_v){
                                found_to_vertex = true;
                                break;
                            }
                        }
                        if(!found_to_vertex){
                            auto is_deg = OVMtetToCGALtet(mesh_, domain_vertex_position_prop_, *vc_it).is_degenerate();

                            if(is_deg){
                                valid_in_domain = false;
                                //std::cout<<" --> collapse would create degenerate tets, skipping"<<std::endl);
                                break;
                            }
                        }
                        if(OVMtetToCGALtet(mesh_, domain_vertex_position_prop_, *vc_it).orientation() == CGAL::NEGATIVE){
                            valid_in_domain = false;
                            //std::cout<<" --> collapse would create degenerate tets, skipping"<<std::endl);
                            break;
                        }

                    }
                    domain_vertex_position_prop_[from_v] = initial_domain_pos;

                    if(valid_in_domain){
                        int displacement_result(1);
                        VertexPosition new_pos;

                        auto one_ring_setup_start_time = std::chrono::high_resolution_clock::now();
                        ExpansionCone cone;
                        ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(mesh_,
                                                                                    vertex_position_prop_,
                                                                                    {from_v, to_v},
                                                                                    cone);
                        auto one_ring_setup_end_time = std::chrono::high_resolution_clock::now();
                        float one_ring_setup_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(one_ring_setup_end_time - one_ring_setup_start_time).count() / 1000000;
                        //PRINT_IF_NOT_SILENT(" - 1-ring setup duration for edge "<<idx<<": "<<one_ring_setup_duration_s<<"s"<<std::endl);
                        //PRINT_IF_NOT_SILENT(" - 1-ring size = "<<cone.n_vertices()<<std::endl);
                        int max, total;
                        double avg;
                        VertexHandle max_precision_vh;
                        compute_precision_stats(cone, cone.vertex_position_prop(), max, total, avg, max_precision_vh);
                        //PRINT_IF_NOT_SILENT(" - 1-ring max precision = "<<max<<", total = "<<total<<", average = "<<avg<<std::endl);


                        if(mode){

                            //PRINT_IF_NOT_SILENT(" - new pos = "<<vec2vec(new_pos)<<std::endl);

                            auto cone_from_v = cone.mesh_to_cone_handle(from_v);
                            auto   cone_to_v = cone.mesh_to_cone_handle(to_v);

                            auto cone_he = cone.find_halfedge(cone_from_v, cone_to_v);
                            if(!cone_he.is_valid()){
                                PRINT_IF_NOT_SILENT(" ERROR - couldn't find cone halfedge "<<cone_from_v<<"-"<<cone_to_v<<std::endl);
                                return -1;
                            }
                            //PRINT_IF_NOT_SILENT(" cone he: "<<cone_he<<": "<<cone.find_halfedge(cone_he)<<std::endl);
                            TopoHelper cone_topo_helper(cone);
                            //PRINT_IF_NOT_SILENT(" cone collapsible: "<<cone_collapsible<<std::endl);
                            cone.collapse_edge(cone_he);
                            //PRINT_IF_NOT_SILENT(" new v = "<<new_v<<std::endl);

                            ExactBadTetFinder bad_tet_finder(cone, cone.vertex_position_prop());

                            const int max_step_count(EDGE_COLLAPSE_MAX_STEP_COUNT);
                            int step_count(0);

                            for(int i(0); i <= max_step_count; i++){

                                step_count++;
                                double alpha = (double)i / (double)step_count;
                                new_pos = (alpha * vertex(from_v) + (1 - alpha) * vertex(to_v));
                                cone.set_vertex(cone_to_v, new_pos);

                                displacement_result = 0;
                                //std::cout<<" - checking validity for step "<<step_count<<" for edge "<<codomain_mesh_.find_halfedge(cone_he)<<std::endl);
                                for(auto vc_it = cone.vc_iter(cone_to_v); vc_it.valid(); vc_it++){
                                    //PRINT_IF_NOT_SILENT(" - checking cell "<<cone.get_cell_vertices(*vc_it)<<std::endl);
                                    if(bad_tet_finder.isFlipped(*vc_it) || bad_tet_finder.isDegenerate(*vc_it)){
                                        displacement_result = 1;
                                        break;
                                    }
                                }

                                if(!displacement_result){
                                    average_steps_count += step_count;
                                    //PRINT_IF_NOT_SILENT(" --> found valid position after "<<i<<" step "<<std::endl);
                                    break;
                                }
                            }


                            //also try with centroid of 1-ring neighborhood
                            if(displacement_result){
                                VertexPosition centroid(0,0,0);
                                int n(0);
                                for(auto vv_it = cone.vv_iter(cone_to_v); vv_it.valid(); vv_it++){
                                    centroid += cone.vertex(*vv_it);
                                    n++;
                                }
                                centroid /= n;

                                cone.set_vertex(cone_to_v, centroid);

                                displacement_result = 0;
                                for(auto vc_it = cone.vc_iter(cone_to_v); vc_it.valid(); vc_it++){
                                    //PRINT_IF_NOT_SILENT(" - checking cell "<<cone.get_cell_vertices(*vc_it)<<std::endl);
                                    if(bad_tet_finder.isFlipped(*vc_it) || bad_tet_finder.isDegenerate(*vc_it)){
                                        //PRINT_IF_NOT_SILENT(" --> bad tet"<<std::endl);
                                        displacement_result = 1;
                                        break;
                                    }
                                }
                                if(!displacement_result){
                                    //PRINT_IF_NOT_SILENT(" --> worked with centroid"<<std::endl);
                                    centroid_used_count++;
                                }
                                new_pos = centroid;
                            }

                        }else{

                            auto exp_check_start_time = std::chrono::high_resolution_clock::now();
                            displacement_result = cone.is_geo_expandable(new_pos);
                            auto exp_check_end_time = std::chrono::high_resolution_clock::now();
                            float exp_check_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(exp_check_end_time - exp_check_start_time).count() / 1000000;
                            //PRINT_IF_NOT_SILENT(" - 1-ring expandability check duration for edge ("<<from_v<<"-"<<to_v<<"): "<<exp_check_duration_s<<"s"<<std::endl);
                        }

                        if(!displacement_result){
                            //PRINT_IF_NOT_SILENT(" --> union of 1-ring neighborhoods for he "<<codomain_mesh_.find_halfedge(he)<<" is star-shaped"<<std::endl);
                            mesh_.collapse_edge(he);
                            this->set_vertex(to_v, new_pos);
                            collapse_count++;
                            //std::cout<<" --> collapsed edge  "<<codomain_mesh_.find_halfedge(he)<<std::endl);

                            //split_list_[i].cone_collapsed_to_vertex = he_to_v;
                            surviving_vertex = to_v;
                            moved_at_last_iteration_prop[to_v] = true;

                            //temp warning disabled that for now
                            //break;
                        }else{
                            uncollapsible_count++;
                        }
                    }
                }

                if(surviving_vertex.is_valid()){
#if ENABLE_ALL_CHECKS
                    if(ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(mesh_, vertex_position_prop_, surviving_vertex)){
                        PRINT_IF_NOT_SILENT(" ERROR - mesh contains flipped tets after collapsing edge "<<from_v<<" -> "<<surviving_vertex<<std::endl);
                        return -1;
                    }
#endif

                    const auto initial_byte_size = byte_size(vertex(surviving_vertex));
                    //PRINT_IF_NOT_SILENT(" - vertex "<<v<<" position byte size = "<<initial_byte_size<<std::endl);
                    if(initial_byte_size >= DEFAULT_POSITION_BYTE_SIZE_THRESHOLD){

                        ExpansionCone cone;
                        ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(mesh_,
                                                                                    vertex_position_prop_,
                                                                                    surviving_vertex,
                                                                                    cone);

                        VertexPosition new_pos = binary_search_minimum_precision_position(surviving_vertex,
                                                                                          cone,
                                                                                          this->vertex(surviving_vertex),
                                                                                          DEFAULT_POSITION_BYTE_SIZE_THRESHOLD);

                        if(byte_size(new_pos) < initial_byte_size){
                            //PRINT_IF_NOT_SILENT(" --> reduced position size from "<<initial_byte_size<<" to "<<byte_size(new_pos)<<" for vertex "<<v<<std::endl);
                            this->set_vertex(surviving_vertex, new_pos);
                        }
                    }

                    if(collapse_count && !(collapse_count % 20)){

                        auto iteration_end_time = std::chrono::high_resolution_clock::now();
                        float iteration_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(iteration_end_time - iteration_start_time).count() / 1000000;

                        PRINT_IF_NOT_SILENT(" -- collapsed "<<collapse_count<<" in "<<iteration_duration_s<<" seconds"<<std::endl);

                    }
                }


                if(collapse_count >= max_collapsed_edges_count){
                    PRINT_IF_NOT_SILENT(" --> reached max number of edges to collapse for this iteration"<<std::endl);
                    break;
                }
                //std::cout<<" what"<<std::endl);
            }

            auto iteration_end_time = std::chrono::high_resolution_clock::now();
            float iteration_duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(iteration_end_time - iteration_start_time).count() / 1000000;

            average_steps_count /= (collapse_count - centroid_used_count);
            PRINT_IF_NOT_SILENT(" -----------------------------------------");
            PRINT_IF_NOT_SILENT(" - collapsed "<<collapse_count<<" edges at iteration "<<(iteration_count + 1)<<" in "<<iteration_duration_s<<"s"<<(timedout ? " (timed out)":"")<<std::endl);
            PRINT_IF_NOT_SILENT(" --> tried "<<tried_count<<" edges: "<<std::endl);
            PRINT_IF_NOT_SILENT("       unchanged valence: "<<unchanged_valence_count<<std::endl);
            PRINT_IF_NOT_SILENT("  unchanged neighborhood: "<<unchanged_neighborhood_count<<std::endl);
            PRINT_IF_NOT_SILENT("           invalid edges: "<<invalid_edges_count<<std::endl);
            PRINT_IF_NOT_SILENT("     uncollapsible edges: "<<uncollapsible_count<<std::endl);
            PRINT_IF_NOT_SILENT("    unexpandable 1-rings: "<<unexp_one_ring_count<<std::endl);
            PRINT_IF_NOT_SILENT("average steps on success: "<<average_steps_count<<std::endl);
            PRINT_IF_NOT_SILENT("          #centroid used: "<<centroid_used_count<<std::endl);

            /*if(!iteration_count && unchanged_valence_count){
                std::cout<<" ERROR - first iteration but unchanged valence vertices..."<<std::endl);
                return -1;
            }*/

            iteration_count++;
            total_collapse_count += collapse_count;

        }while(collapse_count &&
               iteration_count < max_iteration_count &&
               !timedout);

        total_collapsed_edges_count_ += total_collapse_count;

        auto end_time = std::chrono::high_resolution_clock::now();
        float duration_s = (float)std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000000;

        PRINT_IF_NOT_SILENT(" ...done! Collapsed "<<total_collapse_count<<"/"<<to_collapse.size()<<
                            " edges in total in "<<iteration_count<<
                            " iterations in "<<duration_s<<"s"<<(timedout ? " (timed out)":"")<<std::endl);
        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);

        total_edge_collapsing_time_s_ += duration_s;

        bool found_bad_cell(false);
        for(auto c: mesh_.cells()){
            if(mesh_.get_cell_vertices(c).size() != 4){
                PRINT_IF_NOT_SILENT(" - cell "<<c<<" does not have 4 vertices: "<<mesh_.get_cell_vertices(c)<<std::endl);
                found_bad_cell = true;
            }
        }
        if(found_bad_cell){
            PRINT_IF_NOT_SILENT(" ERROR - found at least one cell with more or less than 4 vertices"<<std::endl);
            return -1;
        }

        reset_EC_memory();

        return 0;

    }


    int Expander::split_edges_between_cone_vertices_but_not_part_of_the_cone(ExpansionCone& ss_cone){
        PRINT_IF_NOT_SILENT(" ------------------------------------------------------------------------------"<<std::endl);
        PRINT_IF_NOT_SILENT(" splitting non-cone edges connecting cone vertices..."<<std::endl);

        std::vector<EdgeHandle> to_split;
        auto visited_prop = ss_cone.request_vertex_property<bool>();

        //ss_cone.print_details();

        for(auto cone_v: ss_cone.vertices()){


            visited_prop[cone_v] = true;

            if(ss_cone.is_cone_tip(cone_v)){
                continue;
            }

            auto mesh_v = ss_cone.cone_to_mesh_handle(cone_v);

            //std::cout<<" - visiting vertex "<<cone_v<<" -> "<<mesh_v<<std::endl);

            if(!expanded_prop_[mesh_v]){
                std::cout<<" ERROR - mesh vertex "<<mesh_v<<" is not expanded..."<<std::endl;
                return -1;
            }

            for(auto out_he: mesh_.outgoing_halfedges(mesh_v)){
                auto mesh_neighbor = mesh_.to_vertex_handle(out_he);

                auto cone_neighbor = ss_cone.mesh_to_cone_handle(mesh_neighbor);

                //std::cout<<" -- checking neighbor "<<cone_neighbor<<" -> "<<mesh_neighbor<<std::endl);

                if(cone_neighbor.is_valid() &&
                        !ss_cone.is_deleted(cone_neighbor) &&
                        !visited_prop[cone_neighbor] &&
                        !ss_cone.is_cone_tip(cone_neighbor)){

                    auto cone_he = ss_cone.find_halfedge(cone_v, cone_neighbor);


                    if(!expanded_prop_[mesh_neighbor]){
                        std::cout<<" ERROR - mesh neighbor "<<mesh_neighbor<<" is not expanded..."<<std::endl;
                        return -1;
                    }

                    if(!cone_he.is_valid()){
                        PRINT_IF_NOT_SILENT(" - edge "<<mesh_.halfedge(out_he)<<" is part of the mesh but not the cone, adding for split"<<std::endl);
                        to_split.push_back(mesh_.edge_handle(out_he));
                    }
                }
            }
        }

        PRINT_IF_NOT_SILENT(" found "<<to_split.size()<<" edges to split"<<std::endl);
        for(auto e: to_split){
            auto mesh_from_v = mesh_.edge(e).from_vertex();
            auto mesh_to_v   = mesh_.edge(e).to_vertex();
            auto mid_vertex = split_edge(e);
            auto mid_pos = (vertex(mesh_from_v) + vertex(mesh_to_v))/2;
            //to update the exact vertex position prop
            this->set_vertex(mid_vertex, mid_pos);

            split_list_.push_back({mesh_from_v,
                                  mesh_to_v,
                                  mid_vertex,
                                  mid_pos});

            //expanded_prop_[mid_vertex] = false;
            set_expanded_prop(mid_vertex, true);
        }
        PRINT_IF_NOT_SILENT(" ...done"<<std::endl);


        PRINT_IF_NOT_SILENT(" ------------------------------------------------------------------------------"<<std::endl);


        return 0;
    }

    int Expander::split_edges_of_tets_with_three_faces_on_cone_base(ExpansionCone& ss_cone){

            PRINT_IF_NOT_SILENT(" ------------------------------------------------------------------------------");
            PRINT_IF_NOT_SILENT(" looking for tets that are incident to three faces on the cone's base...");

            std::vector<EdgeHandle> to_split;
            auto split_prop = mesh_.request_edge_property<bool>();

            for(auto cone_e: ss_cone.edges()){
                if(!ss_cone.is_boundary(cone_e)){
                    continue;
                }

                auto cone_from_v = ss_cone.edge(cone_e).from_vertex();
                auto cone_to_v = ss_cone.edge(cone_e).to_vertex();

                if(ss_cone.is_cone_tip(cone_from_v) || ss_cone.is_cone_tip(cone_to_v)){
                    continue;
                }

                auto mesh_from_v = ss_cone.cone_to_mesh_handle(cone_from_v);
                auto mesh_to_v   = ss_cone.cone_to_mesh_handle(cone_to_v);

                auto mesh_he = mesh_.find_halfedge(mesh_from_v, mesh_to_v);

                if(!mesh_he.is_valid()){
                    std::cout<<" ERROR - mesh edge "<<mesh_from_v<<" - "<<mesh_to_v<<" doesn't exist"<<std::endl;
                    return -1;
                }

                for(auto hehf_it = mesh_.hehf_iter(mesh_he); hehf_it.valid(); hehf_it++){
                    auto mesh_op_vertex = mesh_.to_vertex_handle(mesh_.next_halfedge_in_halfface(mesh_he, *hehf_it));

                    if(!mesh_op_vertex.is_valid()){
                        std::cout<<" ERROR - invalid mesh op vertex on face "<<mesh_.get_halfface_vertices(*hehf_it)<<std::endl;
                        return -1;
                    }

                    auto cone_op_vertex = ss_cone.mesh_to_cone_handle(mesh_op_vertex);

                    if(cone_op_vertex.is_valid() && !ss_cone.is_cone_tip(cone_op_vertex)){
                        //std::cout<<" --> found cone base vertex "<<cone_op_vertex<<", looking for hf"<<std::endl);
                        auto cone_hf = ss_cone.find_halfface({cone_from_v, cone_to_v, cone_op_vertex});

                        if(!cone_hf.is_valid()){
                            PRINT_IF_NOT_SILENT("  ---> found tricky tet"<<std::endl);
                            to_split.push_back({mesh_.edge_handle(mesh_.find_halfedge(mesh_from_v, mesh_to_v))});
                            to_split.push_back({mesh_.edge_handle(mesh_.find_halfedge(mesh_from_v, mesh_op_vertex))});
                        }

                    }
                }
            }

            PRINT_IF_NOT_SILENT(" found "<<to_split.size()<<" edges to split"<<std::endl);
            for(auto e: to_split){
                if(split_prop[e]){
                    continue;
                }
                auto mesh_from_v = mesh_.edge(e).from_vertex();
                auto mesh_to_v   = mesh_.edge(e).to_vertex();
                auto mid_vertex = split_edge(e);
                auto mid_pos = (vertex(mesh_from_v) + vertex(mesh_to_v))/2;
                //to update the exact vertex position prop
                this->set_vertex(mid_vertex, mid_pos);

                PRINT_IF_NOT_SILENT(" - split edge "<<mesh_.edge(e)<<" -> "<<mid_vertex<<std::endl);

                split_list_.push_back({mesh_from_v,
                                      mesh_to_v,
                                      mid_vertex,
                                      mid_pos});

                //expanded_prop_[mid_vertex] = false;
                set_expanded_prop(mid_vertex, true);

                auto cone_from_v = ss_cone.mesh_to_cone_handle(mesh_from_v);
                auto cone_to_v   = ss_cone.mesh_to_cone_handle(mesh_to_v);
                auto cone_e      = ss_cone.edge_handle(ss_cone.find_halfedge(cone_from_v, cone_to_v));

                auto cone_mid_v = ss_cone.split_edge(cone_e, mid_vertex);

                PRINT_IF_NOT_SILENT(" - cone mid-vertex is "<<cone_mid_v<<" -> "<<ss_cone.cone_to_mesh_handle(cone_mid_v)<<" at "<<vec2vec(ss_cone.vertex(cone_mid_v))<<std::endl);

                split_prop[e] = true;

            }
            PRINT_IF_NOT_SILENT(" ...done"<<std::endl);
            PRINT_IF_NOT_SILENT(" ------------------------------------------------------------------------------"<<std::endl);

            return 0;
        }



    int Expander::star_shapify_cluster(const ExpansionCone& cluster_ec){

#if !CLUSTER_SS_ENABLED
        PRINT_IF_NOT_SILENT(" CLUSTER STAR-SHAPIFICATION DISABLED -> STOPPING THERE");
        return 1;
#endif

        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);
        PRINT_IF_NOT_SILENT(" - STAR-SHAPIFYING CLUSTER "<<cluster_ec<<"..."<<std::endl);

       // PRINT_IF_NOT_SILENT(" cluster EC: "); if(!silent_mode_){cluster_ec.print_details();}


        /*VertexPosition blah;
        auto sadfasd_copy = cluster_ec;
        std::cout<<" cluster EC exp: "<<sadfasd_copy.is_expandable(blah, true)<<std::endl);*/

        //make a new "virtual" EC with a single tip vertex
        ExpansionCone virtual_cone = cluster_ec;

        std::vector<VertexHandle> mesh_tip_vertices;
        //stores all halfedges going from a base boundary vertex to a tip vertex
        auto tip_neighbors_prop = virtual_cone.request_vertex_property<std::vector<VertexHandle>>();
        for(auto cone_v: virtual_cone.vertices()){
            if(!virtual_cone.is_cone_tip(cone_v)/* && base_boundary_prop[cone_v]*/){
                for(auto out_he: virtual_cone.outgoing_halfedges(cone_v)){
                    //auto mesh_neighbor = virtual_cone.to_vertex_handle(out_he);
                    //auto cone_neighbor = cluster_ec.mesh_to_cone_handle(mesh_neighbor);
                    auto cone_neighbor = virtual_cone.to_vertex_handle(out_he);
                    if(cone_neighbor.is_valid() &&
                       virtual_cone.is_cone_tip(cone_neighbor)){
                        tip_neighbors_prop[cone_v].push_back(cone_neighbor);
                    }
                }
            }else{
                auto mesh_v = virtual_cone.cone_to_mesh_handle(cone_v);
                mesh_tip_vertices.push_back(mesh_v);
            }
        }

        /*PRINT_IF_NOT_SILENT(" tip neighbors: "<<std::endl);
        for(auto cone_v: virtual_cone.vertices()){
            PRINT_IF_NOT_SILENT(" - "<<cone_v<<" -> "<<tip_neighbors_prop[cone_v]<<std::endl);
        }*/


        auto new_tip = virtual_cone.merge_tips();



        //star-shapify the virtual cone
        StarShapifyableExpansionCone ss_cone(virtual_cone, remaining_seconds_before_timeout());

        const int original_cone_max_vertex_index(ss_cone.n_vertices() -1);

        PRINT_IF_NOT_SILENT(" - original cone max vertex index: "<<original_cone_max_vertex_index<<std::endl);

        auto ss_result = ss_cone.star_shapify(false, true);


        if(ss_result){
            PRINT_IF_NOT_SILENT(" CSS with precision reduction failed, trying without"<<std::endl);
            ss_cone =  StarShapifyableExpansionCone(virtual_cone, remaining_seconds_before_timeout());

            ss_result = ss_cone.star_shapify(false, false);
        }

        if(ss_result){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't star-shapify virtual cluster EC"<<std::endl);
            return -1;
        }



        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);
        PRINT_IF_NOT_SILENT("   DUPLICATING STAR-SHAPIFICATION SPLITS TO CLUSTER CONE..."<<std::endl);




        //auto ss_cone_copy = ss_cone;
        auto cluster_ec_copy = cluster_ec;
        //adding a dummy vertex
        auto dummy_vertex = cluster_ec_copy.add_vertex(VertexHandle(-1), {0,0,0});
        cluster_ec_copy.delete_vertex(dummy_vertex);

        //auto base_vertex_mid_vertex_index_prop = ss_cone.request_vertex_property<int>();


        /*auto ss_cone_tip_neighbors_prop = ss_cone.request_vertex_property<std::vector<VertexHandle>>();
        for(auto v: virtual_cone.vertices()){
            ss_cone_tip_neighbors_prop[v] = tip_neighbors_prop[v];
        }*/

        auto ss_cone_base_replacements_prop = ss_cone.request_vertex_property<std::vector<VertexHandle>>();
        //only setting for the cluster EC because we only want the original base vertices
        for(auto v: ss_cone.vertices()){
            if(!ss_cone.is_cone_tip(v)){
                ss_cone_base_replacements_prop[v] = {v};
            }
        }

        auto split_layer_idx_prop = cluster_ec_copy.request_vertex_property<int>();

        //dimensions 0 -> layer, 1 -> cluster, 2 -> cluster vertex
        std::vector<std::vector<std::vector<VertexHandle>>> cone_mid_vertices_clusters_layers;

        std::vector<VertexHandle> full_cone_mid_vertices_list;
        auto cone_mid_vertices_initial_positions_prop = cluster_ec_copy.request_vertex_property<VertexPosition>();

        //used to update the split list with the correct positions after expanding all mid-vertices
        //(which is done after the split list application)
        auto cone_mid_vertex_to_split_list_index_prop = cluster_ec_copy.request_vertex_property<int>();
        //std::vector<int> mid_vertices_indices_in_split_list;
        const int start_index_in_split_list(split_list_.size());
        PRINT_IF_NOT_SILENT(" - starting index in split list = "<<start_index_in_split_list<<std::endl);
        PRINT_IF_NOT_SILENT(" - virtual cone split list size: "<<ss_cone.get_split_list().size()<<std::endl);

        if(ss_cone.get_split_list().size() > 100){
            PRINT_IF_NOT_SILENT(" - too many splits, not printing full list"<<std::endl);
        }else{
            PRINT_IF_NOT_SILENT(" - splits: "<<std::endl);
            for(const auto& split: ss_cone.get_split_list()){
                /*if(split.cone_from_vertex == new_tip ||
                    split.cone_to_vertex == new_tip){*/
                PRINT_IF_NOT_SILENT(" -- "<<split<<std::endl)
                        //}
            }
        }

        //update the split list by duplicating the spoke edge splits
        SplitList cluster_split_list;
        int new_vertex_index(ss_cone.n_vertices());
        PRINT_IF_NOT_SILENT(" - next new vertex index is "<<new_vertex_index<<std::endl);
        //int first_spoke_split_index(-1);
        std::vector<VertexHandle> base_mid_vertices;
        int k(0);
        for(const auto& split: ss_cone.get_split_list()){
            //PRINT_IF_NOT_SILENT(" ----------------------------------"<<std::endl);
            //PRINT_IF_NOT_SILENT(" - processing split n°"<<k<<": "<<split<<std::endl);

            VertexHandle cone_base_vertex(-1);
            if(split.cone_from_vertex == new_tip){
                cone_base_vertex = split.cone_to_vertex;
            }else if(split.cone_to_vertex == new_tip){
                cone_base_vertex = split.cone_from_vertex;
            }


            if(cone_base_vertex.is_valid()){
                //std::cout<<" WARNING - HARD BYPASS OF SPOKE DUPLICATION"<<std::endl);
                //continue;

                const int base_vertex_layer = split_layer_idx_prop[cone_base_vertex];
                const int current_layer_idx = base_vertex_layer + 1;

                //PRINT_IF_NOT_SILENT(" --> base vertex is "<<cone_base_vertex<<
                //                    " from layer "<<base_vertex_layer<<
                //                    " -> spoke edge, duplicating it "<<std::endl);

                auto mid_vertex_position = split.new_vertex_position;
                //PRINT_IF_NOT_SILENT(" -- mid-vertex position = "<<vec2vec(mid_vertex_position)<<std::endl);

                std::vector<VertexHandle> mid_vertices_cluster;
                int i(0);
                for(auto cone_base_replacement: ss_cone_base_replacements_prop[cone_base_vertex]){
                    //PRINT_IF_NOT_SILENT(" -- handling cone base replacement n°"<<i<<": "<<cone_base_vertex<<" -> "<<cone_base_replacement<<std::endl);
                    //for(int i(0); i< (int)ss_cone_tip_neighbors_prop[base_vertex].size(); i++){
                    std::vector<HalfEdgeHandle> spokes_to_split;
                    for(auto out_he: cluster_ec_copy.outgoing_halfedges(cone_base_replacement)){
                        //auto cone_tip_v = ss_cone_tip_neighbors_prop[base_vertex][i];
                        auto cone_tip_v = cluster_ec_copy.to_vertex_handle(out_he);
                        if(!cluster_ec_copy.is_cone_tip(cone_tip_v)){
                            continue;
                        }

                        //PRINT_IF_NOT_SILENT(" --- found spoke edge "<<cluster_ec_copy.find_halfedge(out_he)<<", splitting it "<<std::endl);

                        spokes_to_split.push_back(out_he);
                    }

                    for(auto to_split: spokes_to_split){
                        auto cone_new_mid_vertex = cluster_ec_copy.split_edge(to_split);
                        cluster_ec_copy.set_vertex(cone_new_mid_vertex, mid_vertex_position);
                        cluster_split_list.push_back({cone_base_replacement,
                                                      cluster_ec_copy.to_vertex_handle(to_split),
                                                      cone_new_mid_vertex,
                                                      mid_vertex_position});
                        cone_mid_vertex_to_split_list_index_prop[cone_new_mid_vertex] = start_index_in_split_list + cluster_split_list.size() -1;



                        //if(i){
                        full_cone_mid_vertices_list.push_back(cone_new_mid_vertex);
                        cone_mid_vertices_initial_positions_prop[cone_new_mid_vertex] = mid_vertex_position;

                        mid_vertices_cluster.push_back(cone_new_mid_vertex);
                        split_layer_idx_prop[cone_new_mid_vertex] = current_layer_idx;
                        //PRINT_IF_NOT_SILENT(" -----> added new mid-vertex "<<cone_new_mid_vertex<<
                        //                    " for expansion, layer "<<split_layer_idx_prop[cone_new_mid_vertex]<<std::endl);


                    }
                    i++;
                }
                if(mid_vertices_cluster.size() > 1){
                //if(!mid_vertices_cluster.empty()){

                    while((int)cone_mid_vertices_clusters_layers.size() <= current_layer_idx){
                        cone_mid_vertices_clusters_layers.push_back({});
                    }
                    cone_mid_vertices_clusters_layers[current_layer_idx].push_back(mid_vertices_cluster);
                }

                //ss_cone_tip_neighbors_prop[split.cone_new_vertex] = ss_cone_tip_neighbors_prop[cone_base_vertex];
                ss_cone_base_replacements_prop[split.cone_new_vertex] = mid_vertices_cluster;
                //PRINT_IF_NOT_SILENT(" --- set replacement vertices for layer "<<current_layer_idx<<
                //                    " mid vertex "<<split.cone_new_vertex<<
                //                    " as "<<mid_vertices_cluster<<std::endl);
                if(mid_vertices_cluster.empty()){
                    std::cout<<" ERROR - no replacement vertices for vertex "<<split.cone_new_vertex<<std::endl;
                    return -1;
                }
                //base_vertex_mid_vertex_index_prop[original_base_vertex]++;
                //PRINT_IF_NOT_SILENT(" --> cone tip neighbors for mid-vertex "<<split.cone_new_vertex<<": "<<ss_cone_tip_neighbors_prop[split.cone_new_vertex]<<std::endl);

            }else{
                PRINT_IF_NOT_SILENT(" --> split "<<split<<" not involving the tip vertex, simply copying it"<<std::endl);
                //cluster_split_list.push_back(split);
                auto base_edge = cluster_ec_copy.find_halfedge(split.cone_from_vertex, split.cone_to_vertex);
                if(!base_edge.is_valid()){
                    std::cout<<" ERROR - couldn't find non-spoke edge to split ("<<split.cone_from_vertex<<", "<<split.cone_to_vertex<<")"<<std::endl;
                    return -1;
                }

                base_mid_vertices.push_back(split.cone_new_vertex);

                auto cone_new_mid_vertex = cluster_ec_copy.split_edge(base_edge);
                auto mid_vertex_position = split.new_vertex_position;

                VertexHandle cone_collapsed_to_vertex = split.cone_collapsed_to_vertex;
                if(cone_collapsed_to_vertex.is_valid()){
                    PRINT_IF_NOT_SILENT(" --> split was followed by a collapse, forwarding it to the cluster cone"<<std::endl);
                    auto to_collapse = cluster_ec_copy.find_halfedge(split.cone_new_vertex, split.cone_collapsed_to_vertex);
                    if(!to_collapse.is_valid()){
                        std::cout<<" ERROR - couldn't recover edge to collapse ("<<split.cone_new_vertex<<" - "<<split.cone_collapsed_to_vertex<<")"<<std::endl;
                        return -1;
                    }
                    cluster_ec_copy.collapse_edge(to_collapse);

                }

                cluster_ec_copy.set_vertex(cone_new_mid_vertex, mid_vertex_position);
                cluster_split_list.push_back({cluster_ec_copy.from_vertex_handle(base_edge),
                                              cluster_ec_copy.to_vertex_handle(base_edge),
                                              cone_new_mid_vertex,
                                              mid_vertex_position,
                                              split.cone_collapsed_to_vertex});
                cone_mid_vertex_to_split_list_index_prop[cone_new_mid_vertex] = start_index_in_split_list + cluster_split_list.size() -1;
                //PRINT_IF_NOT_SILENT(" --- copied base split "<<cluster_split_list.back()<<" to cluster split list at index "<<cone_mid_vertex_to_split_list_index_prop[cone_new_mid_vertex]<<std::endl);

                //add the new base mid-vertex as its own replacement because it wasn't part of the original
                //cone, and thus the property wasn't initialized for it before the duplication
                ss_cone_base_replacements_prop[split.cone_new_vertex] = {split.cone_new_vertex};
                //PRINT_IF_NOT_SILENT(" --- set base replacement for vertex "<<split.cone_new_vertex<<" as "<<ss_cone_base_replacements_prop[split.cone_new_vertex]<<std::endl);

            }

            k++;
        }

        PRINT_IF_NOT_SILENT(" - cluster split list size: "<<cluster_split_list.size()<<std::endl);

        PRINT_IF_NOT_SILENT(" ... DONE WITH SPLITS DUPLICATION");
        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);




        //auto cone_to_mesh_prop = codomain_mesh_.request_vertex_property<VertexHandle>();
        //the cone copy should have the right size
        auto cone_to_mesh_prop = cluster_ec_copy.request_vertex_property<VertexHandle>();
        for(int i(0); i<(int)ss_cone.n_vertices(); i++){
            auto v = VertexHandle(i);
            cone_to_mesh_prop[v] = cluster_ec_copy.cone_to_mesh_v_handle_prop()[v];
            //PRINT_IF_NOT_SILENT(" - set prop "<<v<<" -> "<<cone_to_mesh_prop[v]<<" and split-list index for mesh vertex "<<cone_to_mesh_prop[v]<<" as "<<mid_vertex_index_in_split_list_prop[cone_to_mesh_prop[v]]<<std::endl);
        }
        int ss_app_result = apply_split_list(cluster_split_list,
                                             cone_to_mesh_prop,
                                             true);

        if(ss_app_result){
            PRINT_IF_NOT_SILENT(" - ERROR: couldn't apply star-shapification to mesh"<<std::endl);
            if(ss_cone.n_cells() > 50){
                PRINT_IF_NOT_SILENT(" --> cone contains "<<ss_cone.n_cells()<<" cells. Change this line if you really want to print it"<<std::endl);
            }else{
                ss_cone.print_details();
            }
            return -1;
        }

#if ENABLE_ALL_CHECKS
        if(ExactBadTetFinder::meshContainsFlippedTets(mesh_, vertex_position_prop_)){
            std::cout<<" ERROR - mesh contains flipped tets after applying cluster star-shapification"<<std::endl;
            auto bad_tets = ExactBadTetFinder::findBadTets(mesh_, vertex_position_prop_);
            for(auto flipped_tet: bad_tets.second){
                std::cout<<" - "<<flipped_tet<<": "<<std::endl;
                for(auto v: mesh_.get_cell_vertices(flipped_tet)){
                    std::cout<<"  -- "<<v<<" at "<<vec2vec(vertex_position_prop_[v])<<std::endl;
                }
            }

            std::cout<<" - base mid-vertices: "<<base_mid_vertices<<std::endl;
            for(auto v: base_mid_vertices){
                std::cout<<" -- cells around "<<v<<" : "<<std::endl;
                for(auto vc_it = ss_cone.vc_iter(v); vc_it.valid(); vc_it++){
                    std::cout<<"  --- "<<(*vc_it)<<" (volume = "<<CGAL::to_double(OVMtetToCGALtet(ss_cone, *vc_it).volume())<<") : "<<std::endl;
                    for(auto cv: ss_cone.get_cell_vertices(*vc_it)){
                        std::cout<<"  ---- "<<cv<<" -> "<<cone_to_mesh_prop[cv]<<std::endl;
                    }
                }
            }
            return EXPANSION_ERROR;
        }
#endif
        PRINT_IF_NOT_SILENT(" --> split list application succesful, no flipped tets!"<<std::endl);


        //cluster expandability check + getting the new position
        std::vector<VertexHandle> unexpanded_vertices;
        std::vector<int> cluster_indices;
        for(auto cone_tip: cluster_ec.cone_tip_vertices()){
            unexpanded_vertices.push_back(cluster_ec.cone_to_mesh_handle(cone_tip));
            cluster_indices.push_back(cluster_indices.size());
        }
        ExpansionCone cluster;
        auto cluster_result = merge_expansion_cones(unexpanded_vertices,
                                                    cluster_indices,
                                                    cluster);

        if(cluster_result){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't merge cluster ");
            for(auto c: cluster_indices){
                PRINT_IF_NOT_SILENT(unexpanded_vertices[c]<<" ");
            }
            PRINT_IF_NOT_SILENT(std::endl);
            return -1;
        }

        VertexPosition new_position;
        auto cluster_expandability = cluster.is_expandable(new_position);

        if(cluster_expandability){
            PRINT_IF_NOT_SILENT(" ERROR - cluster is still not expandable after split-list application, code: "<<cluster_expandability<<std::endl);
            cluster.print_details();

            PRINT_IF_NOT_SILENT(" original cluster: "<<std::endl);
            cluster_ec.print_details();
            PRINT_IF_NOT_SILENT(" cone for each base vertex: "<<std::endl);
            for(auto cone_v:cluster_ec.vertices()){
                if(!cluster_ec.is_cone_tip(cone_v)){
                    PRINT_IF_NOT_SILENT(" - base vertex "<<cone_v<<" : "<<std::endl);
                    for(auto vv_it = cluster_ec.vv_iter(cone_v); vv_it.valid(); vv_it++){
                        if(cluster_ec.is_cone_tip(*vv_it)){
                            PRINT_IF_NOT_SILENT("      - "<<(*vv_it)<<std::endl);
                        }
                    }
                }
            }
            PRINT_IF_NOT_SILENT(" cone tip manifoldness: "<<std::endl);
            for(auto cone_v:cluster_ec.cone_tip_vertices()){
                PRINT_IF_NOT_SILENT(" ---------------------------------- "<<std::endl);
                TopoHelper::manifoldVertex(cluster_ec, cone_v);
            }
            //auto cluster_expandability = cluster.is_topo_expandable(true);
            return -1;
        }

        for(auto sub_v: cluster_ec.cone_tip_vertices()){
            auto mesh_v = cluster_ec.cone_to_mesh_handle(sub_v);
            this->set_vertex(mesh_v, new_position);
            PRINT_IF_NOT_SILENT(" - tip vertex "<<sub_v<<" -> "<<mesh_v<<" now at "<<vec2vec(this->vertex(mesh_v))<<" is expanded: "<<expanded_prop_[mesh_v]<<std::endl);
        }
        PRINT_IF_NOT_SILENT(" --> cluster is now expandable, performing mid-vertices expansion..."<<std::endl);

        //then, mid-vertices cluster expansion
        for(int i(0); i< (int)cone_mid_vertices_clusters_layers.size(); i++){
            const auto& cone_mid_vertices_clusters = cone_mid_vertices_clusters_layers[i];
            for(const auto& cone_cluster: cone_mid_vertices_clusters){
                for(auto cone_v: cone_cluster){
                    auto mesh_v = cone_to_mesh_prop[cone_v];
                    //expanded_prop_[mesh_v] = false;
                    set_expanded_prop(mesh_v, false);
                }
            }
        }
        if(cone_mid_vertices_clusters_layers.empty()){
            std::cout<<" ERROR: no mid-vertices for some reason..."<<std::endl;
            return -1;
        }

#define ENABLE_LAYER_STUFF true

        bool full_expansion_successful(ENABLE_LAYER_STUFF);

#if ENABLE_LAYER_STUFF
        PRINT_IF_NOT_SILENT(" - expanding "<<cone_mid_vertices_clusters_layers.size()<<" cone mid-vertices layers "<<std::endl);
        //std::cout<<" WARNING -  expanding them in reverse order (shouldn't make any difference though)"<<std::endl);
        //for(const auto& cone_mid_vertices_cluster: cone_mid_vertices_clusters){
        for(int i(cone_mid_vertices_clusters_layers.size()-1); i>=0; i--){
        //for(int i(0); i< (int)cone_mid_vertices_clusters_layers.size(); i++){

            if(!full_expansion_successful){
                break;
            }

            const auto& cone_mid_vertices_clusters = cone_mid_vertices_clusters_layers[i];
            int cluster_idx(0);

            PRINT_IF_NOT_SILENT(" -------------------------------- expanding layer n°"<<i<<std::endl);
            //for(int k(cone_mid_vertices_clusters.size()-1); k>=0; k--){
            for(auto k(0); k<(int)cone_mid_vertices_clusters.size(); k++){

                if(!full_expansion_successful){
                    break;
                }

                const auto& cone_mid_vertices_cluster = cone_mid_vertices_clusters[k];

                cluster_idx = k;
                PRINT_IF_NOT_SILENT(" -- expanding cone cluster n°"<<cluster_idx<<
                                    ": "<<cone_mid_vertices_cluster<<std::endl);


                int local_expanded_count(0);
                int i(0);
                while(local_expanded_count < (int)cone_mid_vertices_cluster.size() &&
                      i < (int)cone_mid_vertices_cluster.size() &&
                      full_expansion_successful){
                    //PRINT_IF_NOT_SILENT(" ----- "<<i<<"-th mid-vertices expansion iteration"<<std::endl);
                    bool expanded_something(false);
                    for(int j(0); j< (int)cone_mid_vertices_cluster.size(); j++){
                        auto cone_v = cone_mid_vertices_cluster[j];
                        auto mesh_v = cone_to_mesh_prop[cone_v];
                        if(!mesh_v.is_valid()){
                            PRINT_IF_NOT_SILENT(" ERROR - couldn't find mesh vertex for cone mid-vertex "<<cone_v<<std::endl);
                            return -1;
                        }

                        if(expanded_prop_[mesh_v]){
                            //PRINT_IF_NOT_SILENT(" -> already expanded, skipping"<<std::endl);
                            continue;
                        }

                        check_for_timeout();

                        //PRINT_IF_NOT_SILENT(" -- trying unexpanded vertex "<<cone_v<<" -> "<<mesh_v<<std::endl);

                        auto to_expand = mesh_v;

                        //std::cout<<" - setting up 1-ring of vertex "<<to_expand<<std::endl);
#warning TODO: replace with EC::set_up_1_ring...
                        ExpansionCone cone;
                        auto added_vertex_prop = mesh_.request_vertex_property<bool>();
                        cone.add_vertex(to_expand, vertex(to_expand));
                        added_vertex_prop[to_expand] = true;

                        for(auto vc_it = mesh_.vc_iter(to_expand); vc_it.valid(); vc_it++){
                            auto c_vertices = mesh_.get_cell_vertices(*vc_it);
                            for(auto v: c_vertices){
                                if(!added_vertex_prop[v]){
                                    cone.add_vertex(v, vertex(v));
                                    added_vertex_prop[v] = true;
                                }
                            }

                            auto cone_c = cone.add_cell(*vc_it, c_vertices);
                            if(cone_c.idx() == -1){
                                PRINT_IF_NOT_SILENT(" error while adding cluster interface cell "<<c_vertices<<" incident to mid-vertex "<<to_expand<<std::endl);
                                return -1;
                            }
                        }

                        //PRINT_IF_NOT_SILENT(" - 1-ring neighborhood of "<<to_expand<<": "); cone.print_details();

                        VertexPosition pos;
                        auto exp_result = cone.is_geo_expandable(pos);
                        if(exp_result == -1){
                            PRINT_IF_NOT_SILENT(" --> error while checking expandability of 1-ring neighborhood of mid-vertex "<<to_expand<<std::endl);
                            PRINT_IF_NOT_SILENT(" cone details: "); if(!silent_mode_){cone.print_details();}
                            return -1;
                        }else if(exp_result){
                            //PRINT_IF_NOT_SILENT(" --> couldn't expand 1-ring neighborhood of mid-vertex "<<to_expand<<", result: "<<exp_result<<std::endl);


                        }else{
                           //PRINT_IF_NOT_SILENT(" --> mid-vertex "<<to_expand<<" is expandable!"<<std::endl);

                            set_vertex(to_expand, pos);
                            //expanded_prop_[to_expand] = true;
                            set_expanded_prop(to_expand, true);
                            //split.new_vertex_position = pos;
                            //split_list_.push_back(split);
                            local_expanded_count++;

                            if(split_list_[cone_mid_vertex_to_split_list_index_prop[cone_v]].cone_new_vertex != to_expand){
                                PRINT_IF_NOT_SILENT(" ERROR - cone mid-vertex "<<cone_v<<" index in split list is "<<split_list_[cone_mid_vertex_to_split_list_index_prop[cone_v]].cone_new_vertex<<" but expanded vertex is "<<to_expand<<" -> "<<mesh_v<<std::endl);
                                return -1;
                            }
                            split_list_[cone_mid_vertex_to_split_list_index_prop[cone_v]].new_vertex_position = pos;
                            //RINT_IF_NOT_SILENT(" --> updated split n° "<<cone_mid_vertex_to_split_list_index_prop[cone_v]<<
                            //                    " mid-vertex "<<to_expand<<" position as "<<vec2vec(pos)<<std::endl);

                            expanded_something = true;
                            //break;

                            //PRINT_IF_NOT_SILENT(" --> added split "<<split_list_.back()<<" to expander split list"<<std::endl);
                            //PRINT_IF_NOT_SILENT(" --> expanded count: "<<local_expanded_count<<"/"<<local_split_list.size()<<std::endl);
                        }
                    }

                    if(!expanded_something){

                        PRINT_IF_NOT_SILENT(" couldn't expand a single vertex at iteration "<<i<<", resorting to random-ordered expansions"<<std::endl);

                        full_expansion_successful = false;
                        break;


                    }

                    i++;
                }


                if(local_expanded_count != (int)cone_mid_vertices_cluster.size()){
                    PRINT_IF_NOT_SILENT(" couly only expand "<<local_expanded_count<<" mid-vertices"<<std::endl);
                    full_expansion_successful = false;
                    //PRINT_IF_NOT_SILENT(" ERROR - could only expand "<<local_expanded_count<<" mid-vertices among "<<cone_mid_vertices_cluster.size()<<" in mid-vertices cluster"<<std::endl);
                    //return -1;
                }


                cluster_idx++;
            }
            PRINT_IF_NOT_SILENT(" --> succesfully expanded all layer "<<i<<" clusters."<<std::endl);

        }
#endif

        if(!full_expansion_successful){
            LightWeightStopWatch sw;
            const int max_attempts_count(std::numeric_limits<int>::max());

            //std::cout<<" making "<<max_attempts_count<<" attempts at expanding "<< full_cone_mid_vertices_list.size()<<" random-ordered mid-vertices"<<std::endl);

            int i(0);
            //for(int i(0); i<max_attempts_count; i++){
            while(!full_expansion_successful){
                /*if(full_expansion_successful){
                    break;
                }*/

                if(i && !(i%100)){
                    PRINT_IF_NOT_SILENT(" --------------------- trying randomized iteration "<<i<<std::endl);
                }

                for(auto cone_v: full_cone_mid_vertices_list){
                    auto mesh_v = cone_to_mesh_prop[cone_v];
                    set_vertex(mesh_v, cone_mid_vertices_initial_positions_prop[cone_v]);
                    set_expanded_prop(mesh_v, false);
                }


                //std::random_device rd;
                std::mt19937 g(RAND_SEED);
                std::shuffle(full_cone_mid_vertices_list.begin(), full_cone_mid_vertices_list.end(), g);


                int expanded_count(0);
                bool expanded_something(true);
                int k(0);
                while(expanded_count < (int)full_cone_mid_vertices_list.size() &&
                      expanded_something){

                    expanded_something = false;
                    for(auto cone_v: full_cone_mid_vertices_list){
                        auto mesh_v = cone_to_mesh_prop[cone_v];

                        if(expanded_prop_[mesh_v]){
                            //PRINT_IF_NOT_SILENT(" -> already expanded, skipping"<<std::endl);
                            continue;
                        }

                        //std::cout<<" - setting up 1-ring for mesh vertex "<<mesh_v<<" with valence "<<codomain_mesh_.valence(mesh_v)<<std::endl);
                        ExpansionCone cone;
                        ExpansionCone::set_up_1_ring_neighborhood_as_expansion_cone(mesh_,
                                                                                    vertex_position_prop_,
                                                                                    mesh_v,
                                                                                    cone);


                        VertexPosition pos;
                        auto exp_result = cone.is_geo_expandable(pos);
                        if(exp_result == -1){
                            std::cout<<" error while checking cone expandability"<<std::endl;
                            return -1;
                        }else if(!exp_result){
                            expanded_something = true;

                            set_vertex(mesh_v, pos);
                            set_expanded_prop(mesh_v, true);
                            expanded_count++;

                            if(split_list_[cone_mid_vertex_to_split_list_index_prop[cone_v]].cone_new_vertex != mesh_v){
                                std::cout<<" ERROR - cone mid-vertex "<<cone_v<<" index in split list is "<<split_list_[cone_mid_vertex_to_split_list_index_prop[cone_v]].cone_new_vertex<<" but expanded vertex is "<<mesh_v<<" -> "<<mesh_v<<std::endl;
                                return -1;
                            }
                            split_list_[cone_mid_vertex_to_split_list_index_prop[cone_v]].new_vertex_position = pos;

                            //std::cout<<" --> expanded vertex "<<mesh_v<<std::endl);

                           /* if(! (expanded_count % 20)){
                                std::cout<<" --> expanded "<<expanded_count<<" vertices"<<std::endl);
                            }*/
                        }
                        check_for_timeout();
                    }

                    k++;

                }
                if(expanded_count == (int)full_cone_mid_vertices_list.size()){
                    PRINT_IF_NOT_SILENT(" --> random expansions successful after "<<(i+1)<<"-th attempt in "<<k<<" iterations!"<<std::endl);
                    full_expansion_successful = true;
                }else{
                    //std::cout<<" --> failure, only "<<expanded_count<<" expanded vertices among "<<full_cone_mid_vertices_list.size()<<std::endl);
                }
                if(i && !(i % 100)){
                    PRINT_IF_NOT_SILENT(" - made "<<i<<" random-ordered attempts in "<<sw.lap_duration()<<" seconds"<<std::endl);
                }

                i++;
            }

            if(!full_expansion_successful){
                std::cout<<" ERROR - "<<max_attempts_count<<" weren't enough to fully expand the duplicate vertices..."<<std::endl;
                return -1;
            }
        }


        PRINT_IF_NOT_SILENT(" --> succesfully expanded all mid-vertices. Expanding cluster..."<<std::endl);


        //and finally, cluster expansion
        ExpansionCone cluster2;
        cluster_result = merge_expansion_cones(unexpanded_vertices,
                                               cluster_indices,
                                               cluster2);

        if(cluster_result){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't merge cluster after expanding mid-vertices");
            for(auto c: cluster_indices){
                PRINT_IF_NOT_SILENT(unexpanded_vertices[c]<<" ");
            }
            PRINT_IF_NOT_SILENT(std::endl);
            return -1;
        }
        auto cluster_expansion_result = expand_cluster(cluster2, new_position);
        if(cluster_expansion_result){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't expand cluster after star-shapification and split-list application... "<<std::endl);
            return -1;
        }


        for(auto sub_v: cluster_ec.cone_tip_vertices()){
            auto mesh_v = cluster_ec.cone_to_mesh_handle(sub_v);

            if(!expanded_prop_[mesh_v]){
                std::cout<<" ERROR - tip vertex "<<sub_v<<" -> "<<mesh_v<<" is not expanded after cluster expansion"<<std::endl;
                return -1;
            }
        }


        PRINT_IF_NOT_SILENT(" --> and cluster EC is now expanded!"<<std::endl);
        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);




        PRINT_IF_NOT_SILENT(" ...done! Star-shapified cluster EC"<<std::endl);
        PRINT_IF_NOT_SILENT(" ==============================================================="<<std::endl);

        return 0;
    }





    int Expander::split_expanded_spokes_out_of_expansion_cone(const VertexHandle& center_vertex,
                                                              const ExpansionCone& cone,
                                                              std::vector<VertexHandle>& steiner_vertices){

        steiner_vertices.clear();

        std::vector<HalfEdgeHandle> spokes_to_split;

        //PRINT_IF_NOT_SILENT(" splitting spoke edges that couldn't be added to the expansion cone"<<std::endl);
        //finally, we split all edges connecting vertices that we couldn't add the the expansion cone
        for(auto out_he_it: mesh_.outgoing_halfedges(center_vertex)){
            auto neighbor = mesh_.to_vertex_handle(out_he_it);

            if(expanded_prop_[neighbor]){
                //PRINT_IF_NOT_SILENT(" - checking expanded neighbor "<<neighbor<<std::endl);
                auto cone_handle = cone.mesh_to_cone_handle(neighbor);

                if(cone_handle.idx() == -1){
                    //PRINT_IF_NOT_SILENT(" -- neighbor "<<neighbor<<" is not part of cone, adding edge "<<codomain_mesh_.find_halfedge(out_he_it)<<" to split list"<<std::endl);

                    spokes_to_split.push_back(out_he_it);
                }
            }
        }

        for(auto to_split: spokes_to_split){

            auto mid_vertex = split_edge(to_split);
            //PRINT_IF_NOT_SILENT(" --> mid vertex: "<<mid_vertex<<std::endl);
            if(mid_vertex.idx() == -1){
                PRINT_IF_NOT_SILENT(" ERROR - couldn't split edge "<<mesh_.halfedge(to_split)<<std::endl);
                return -1;
            }

            //expanded_prop_[mid_vertex] = false;
            set_expanded_prop(mid_vertex, false);
            this->set_vertex(mid_vertex, {0,0,0});
            to_expand_count_++;
            steiner_vertices.push_back(mid_vertex);
        }

        return 0;
    }




    VertexExpanse Expander::finalize_expanse(const VertexHandle& center_vertex,
                                             const ExpansionCone& cone,
                                             const VertexPosition& new_vertex_position,
                                             int to_expand_count_increase){

#warning cannot use this one with clusters as it is since the condition is on the presence of degenerates.

        auto cone_copy = cone;
        auto min_precision_position = binary_search_minimum_precision_position(center_vertex,
                                                                               cone_copy,
                                                                               new_vertex_position,
                                                                               DEFAULT_POSITION_BYTE_SIZE_THRESHOLD);


        // min_precision_position = new_vertex_position;

        /*PRINT_IF_NOT_SILENT(" --> moved vertex "<<center_vertex<<
                 ", initial solution byte size = "<<max_size(new_vertex_position)<<
                 ", reduced solution byte size = "<<max_size(min_precision_position)<<std::endl);*/


        /*if(max_size(min_precision_position) < 50){
            PRINT_IF_NOT_SILENT(" -----> position = "<<min_precision_position<<std::endl);
        }*/

        move_vertex_and_mark_as_expanded(center_vertex,
                                         min_precision_position);



        data_logger_.add_evolution_point({(to_expand_count() - expanded_count()),
                                          expanded_count(),
                                          (int)mesh_.n_logical_vertices(),
                                          (int)cone.n_logical_vertices(),
                                          to_expand_count_increase,
                                          last_single_sv_ec_count_,
                                          max_byte_size(min_precision_position),
                                          last_unexp_valence_,
                                          PEHelpers::DBCI_vertices_count(mesh_),
                                          center_vertex.idx()});

        last_single_sv_ec_count_ = 0;
        last_unexp_valence_ = 0;

        //reset the property
        expanse_at_last_iteration_prop_[center_vertex] = VertexExpanse();


        /*PRINT_IF_NOT_SILENT(" -------- expanded vertex "<<center_vertex<<" with neighbors: ";
        for(auto out_he: codomain_mesh_.outgoing_halfedges(center_vertex)){
            PRINT_IF_NOT_SILENT(codomain_mesh_.to_vertex_handle(out_he)<<" ";
        }
        PRINT_IF_NOT_SILENT(std::endl;*/

        return {center_vertex,
                cone,
                min_precision_position};
    }



    void Expander::move_vertex_and_mark_as_expanded(const VertexHandle& tip_vertex,
                                                    const VertexPosition& new_vertex_position){

        //temp remove
        auto current_pos = this->vertex(tip_vertex);

        expanded_count_++;
        //expanded_prop_[tip_vertex] = true;
        set_expanded_prop(tip_vertex, true);
        this->set_vertex(tip_vertex, new_vertex_position);

        if(byte_size(new_vertex_position) > 200){
            PRINT_IF_NOT_SILENT(" -- moved vertex "<<tip_vertex<<" to position of size "<<byte_size(new_vertex_position)<<" (B)"<<std::endl);
        }

        //std::cout<<" -- moved vertex "<<tip_vertex<<" from "<<vec2vec(current_pos)<<" to "<<vec2vec(new_vertex_position)<<std::endl);
    }


#if 0
    VertexPosition Expander::find_minimum_precision_position(const VertexHandle& center_vertex,
                                                             const VertexPosition& new_vertex_position){

        //PRINT_IF_NOT_SILENT(" ======================================================="<<std::endl);
        //PRINT_IF_NOT_SILENT("  looking for minimum precision position with dichotomy..."<<std::endl);

        const auto initial_pos = this->vertex(center_vertex);
        auto current_pos = new_vertex_position;

        int byte_size = max_min_size(new_vertex_position);

        //PRINT_IF_NOT_SILENT(" - max-min size = "<<byte_size<<std::endl);
        //PRINT_IF_NOT_SILENT(" - initial max size = "<<max_size(new_vertex_position)<<std::endl);

        int i(0);

        //times 7 so it's at most below 1000 bytes and not more
        int shift_factor = 7 * byte_size;

        do{
            if(max_byte_size(current_pos) < 200 || shift_factor <= 0){
                break;
            }


            PEHelpers::lower_precision(shift_factor, current_pos);

            //PRINT_IF_NOT_SILENT(" ------ "<<std::endl);
            //PRINT_IF_NOT_SILENT(" ---- max size at iteration "<<i<<": "<<max_size(current_pos)<<std::endl);
            //PRINT_IF_NOT_SILENT(" - shift factor = "<<shift_factor<<std::endl);
            //PRINT_IF_NOT_SILENT(" - pos for vertex "<<center_vertex<<" = "<<current_pos<<std::endl);


            this->set_vertex(center_vertex, current_pos);


            shift_factor -= byte_size;
            i++;
        }while(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_));

        this->set_vertex(center_vertex, initial_pos);

        //PRINT_IF_NOT_SILENT(" ...done, final max size = "<<max_size(current_pos)<<std::endl);
        //PRINT_IF_NOT_SILENT(" ======================================================="<<std::endl);

        return current_pos;

    }
#endif


    VertexPosition Expander::binary_search_minimum_precision_position(const VertexHandle& center_vertex,
                                                                      ExpansionCone& cone,
                                                                      const VertexPosition& new_vertex_position,
                                                                      int position_byte_size_threshold){

        //TEMPORARY
        const auto& tip_vertex = center_vertex;

        //return new_vertex_position;

        const int byte_size_lower_bound(position_byte_size_threshold);
        const int min_shift_factor_diff(32);
        const int min_byte_size_diff(byte_size_lower_bound/2);
        const int max_iterations(10);

        auto initial_byte_size = byte_size(new_vertex_position);

        if(initial_byte_size < byte_size_lower_bound){
            return new_vertex_position;
        }
        const auto initial_pos = this->vertex(center_vertex);

        this->set_vertex(center_vertex, new_vertex_position);

        //check which cells are degenerate with the Chebyshev position
        //std::cout<<" - incident cells:"<<std::endl);
        std::vector<CellHandle> incident_cells;
        std::vector<bool> incident_cells_degeneracy;
        for(auto vc_it = mesh_.vc_iter(center_vertex); vc_it.valid(); vc_it++){
            incident_cells.push_back(*vc_it);
            incident_cells_degeneracy.push_back(OVMtetToCGALtet(mesh_, vertex_position_prop_, *vc_it).is_degenerate());
            //std::cout<<" -- "<<(*vc_it)<<" is deg: "<<incident_cells_degeneracy.back()<<std::endl);
        }

        //first, try by simply using double precision
        auto current_pos = vec2vec(vec2vec(new_vertex_position));
        this->set_vertex(center_vertex, current_pos);
        cone.set_vertex(cone.mesh_to_cone_handle(center_vertex), current_pos);


        bool found_change(false);
        for(int i(0); i<(int)incident_cells.size(); i++){
            auto c = incident_cells[i];
            if(OVMtetToCGALtet(mesh_, vertex_position_prop_, c).is_degenerate() != incident_cells_degeneracy[i]){
                //std::cout<<" - cell "<<c<<" was "<<(incident_cells_degeneracy[i] ? "": "non-")<<" degenerate and is now "<<OVMtetToCGALtet(codomain_mesh_, vertex_position_prop_, c).is_degenerate()<<std::endl);
                found_change = true;
                break;
            }
        }

        //PRINT_IF_NOT_SILENT(" ======================================================="<<std::endl);
        //PRINT_IF_NOT_SILENT(" byte size > "<<byte_size_lower_bound<<", looking for minimum precision position with binary search..."<<std::endl);

        if(!found_change &&
                !ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(mesh_, vertex_position_prop_, tip_vertex)){

            //PRINT_IF_NOT_SILENT(" --> double precision enough, using that"<<std::endl);
            return current_pos;

        }

        current_pos = new_vertex_position;

        //PRINT_IF_NOT_SILENT(" - max-min size = "<<byte_size<<std::endl);
        //PRINT_IF_NOT_SILENT(" - initial max size = "<<max_size(new_vertex_position)<<std::endl);



        int max_shift_factor = PEHelpers::find_maximum_shift_factor(new_vertex_position);
        int min_shift_factor = 0;

        int shift_factor = max_shift_factor;

        bool found_best_factor(false);
        int last_byte_size(initial_byte_size);
        int last_shift_factor(shift_factor);

        int i(0);
        while(!found_best_factor && i < max_iterations){
            check_for_timeout();

            current_pos = new_vertex_position;
            PEHelpers::lower_precision(shift_factor, current_pos);
            int current_byte_size = byte_size(current_pos);

            /*PRINT_IF_NOT_SILENT(" ------ "<<std::endl);
            PRINT_IF_NOT_SILENT(" ---- max size at iteration "<<i<<": "<<current_byte_size<<std::endl);
            PRINT_IF_NOT_SILENT(" - min shift factor = "<<min_shift_factor<<std::endl);
            PRINT_IF_NOT_SILENT(" - max shift factor = "<<max_shift_factor<<std::endl);
            PRINT_IF_NOT_SILENT(" - shift factor = "<<shift_factor<<std::endl);*/


            //PRINT_IF_NOT_SILENT(" - pos for vertex "<<center_vertex<<" = "<<current_pos<<std::endl);

            this->set_vertex(center_vertex, current_pos);
            cone.set_vertex(cone.mesh_to_cone_handle(center_vertex), current_pos);


            bool found_change(false);
            for(int i(0); i<(int)incident_cells.size(); i++){
                auto c = incident_cells[i];
                if(OVMtetToCGALtet(mesh_, vertex_position_prop_, c).is_degenerate() != incident_cells_degeneracy[i]){
                    //std::cout<<" - cell "<<c<<" was "<<(incident_cells_degeneracy[i] ? "": "non-")<<" degenerate and is now "<<OVMtetToCGALtet(codomain_mesh_, vertex_position_prop_, c).is_degenerate()<<std::endl);
                    found_change = true;
                    break;
                }
            }
            //if the precision is too low, we use it as lower bound
            if(found_change || ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(mesh_, vertex_position_prop_, tip_vertex)){
                //     cone.contains_degenerate_tets()){

                max_shift_factor = shift_factor;

                //if it's too high, we use it as UPPER bound
            }else{
                min_shift_factor = shift_factor;
            }

            if((shift_factor <= min_shift_factor && current_byte_size < byte_size_lower_bound) ||
                    std::abs(current_byte_size - last_byte_size) < min_byte_size_diff){

                /*PRINT_IF_NOT_SILENT(" ------------- done: "<<std::endl);
                PRINT_IF_NOT_SILENT(" - min shift factor = "<<min_shift_factor<<std::endl);
                PRINT_IF_NOT_SILENT(" - max shift factor = "<<max_shift_factor<<std::endl);
                PRINT_IF_NOT_SILENT(" - shift factor = "<<shift_factor<<std::endl);
                PRINT_IF_NOT_SILENT(" - current byte size = "<<current_byte_size<<std::endl);
                PRINT_IF_NOT_SILENT(" - current byte size diff = "<<std::abs(current_byte_size - last_byte_size)<<std::endl);
                PRINT_IF_NOT_SILENT(" --> ok"<<std::endl);*/
                found_best_factor = true;
            }

            last_byte_size = current_byte_size;
            last_shift_factor = shift_factor;
            shift_factor = (min_shift_factor + max_shift_factor) / 2;
            if(std::abs(last_shift_factor - shift_factor) < min_shift_factor_diff){
                /*PRINT_IF_NOT_SILENT(" last_shift_factor = "<<last_shift_factor<<", shift factor = "<<shift_factor<<std::endl);
                PRINT_IF_NOT_SILENT(" --> shift factor change = "<<std::abs(last_shift_factor - shift_factor)<<" < "<<min_shift_factor_diff<<" -> ok"<<std::endl);*/
                found_best_factor = true;
            }

            //PRINT_IF_NOT_SILENT(" - updated shift factor = "<<shift_factor<<std::endl);

            i++;
        }

        //re-compute best valid position
        current_pos = new_vertex_position;
        PEHelpers::lower_precision(min_shift_factor, current_pos);

        this->set_vertex(center_vertex, initial_pos);


        //PRINT_IF_NOT_SILENT(" ...done, reduced precision from "<<initial_byte_size<<" bytes to "<<byte_size(current_pos)<<std::endl);
        //PRINT_IF_NOT_SILENT(" ======================================================="<<std::endl);

        return current_pos;


#if 0

        const int byte_size_lower_bound(100);
        const int min_shift_factor_diff(16);
        const int max_iterations(3);

        auto new_vertex_initial_position = vertex(tip_vertex);

        auto initial_byte_size = byte_size(new_vertex_initial_position);
        if(initial_byte_size < byte_size_lower_bound){
            return new_vertex_initial_position;
        }

        //PRINT_IF_NOT_SILENT(" ======================================================="<<std::endl);
        //PRINT_IF_NOT_SILENT("  looking for minimum precision position for new tip vertex "<<tip_vertex<<" initially at "<<vec2vec(vertex(tip_vertex))<<" with binary search..."<<std::endl);

        auto v = new_vertex_initial_position;
        /*if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_) ||
                cone.contains_degenerate_tets()){
            PRINT_IF_NOT_SILENT(" ERROR - integrity already not maintained with initial position"<<std::endl);
            return {0,0,0};
        }*/

        //first, try by simply using double precision
        /*const auto initial_pos = this->vertex(tip_vertex);
        auto current_pos = vec2vec(vec2vec(new_vertex_initial_position));
        this->set_vertex(tip_vertex, current_pos);
        cone.set_vertex(cone.mesh_to_cone_handle(center_vertex), current_pos);

        if(!(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_) ||
                cone.contains_degenerate_tets())){
            PRINT_IF_NOT_SILENT(" --> double precision enough, using that"<<std::endl);
            PRINT_IF_NOT_SILENT(" ======================================================="<<std::endl);
            return current_pos;
        }*/

        auto current_pos = new_vertex_initial_position;

        //PRINT_IF_NOT_SILENT(" - initial byte size = "<<initial_byte_size<<std::endl);
        int max_shift_factor = PEHelpers::find_maximum_shift_factor(new_vertex_initial_position);
        int min_shift_factor = 0;
        int shift_factor = max_shift_factor;

        bool found_best_factor(false);

        int i(0);
        while(!found_best_factor && i < max_iterations){

            current_pos = new_vertex_initial_position;
            PEHelpers::lower_precision(shift_factor, current_pos);

            /*PRINT_IF_NOT_SILENT(" ------ "<<std::endl);
            PRINT_IF_NOT_SILENT(" ---- byte size at iteration "<<i<<"/"<<max_iterations<<" : "<<byte_size(current_pos)<<std::endl);
            PRINT_IF_NOT_SILENT("          min shift factor = "<<min_shift_factor<<std::endl);
            PRINT_IF_NOT_SILENT("          max shift factor = "<<max_shift_factor<<std::endl);
            PRINT_IF_NOT_SILENT("              shift factor = "<<shift_factor<<std::endl);
            */
            this->set_vertex(tip_vertex, current_pos);
            cone.set_vertex(cone.mesh_to_cone_handle(center_vertex), current_pos);

            //if the precision is too low, we use it as lower bound
            if(ExactBadTetFinder::meshContainsFlippedTets(codomain_mesh_, vertex_position_prop_) ||
                    cone.contains_degenerate_tets()){
                max_shift_factor = shift_factor;

                //if it's too high, we use it as UPPER bound
            }else{
                min_shift_factor = shift_factor;
            }

            shift_factor = (min_shift_factor + max_shift_factor) / 2;
            //PRINT_IF_NOT_SILENT(" - updated shift factor = "<<shift_factor<<std::endl);

            if(shift_factor <= min_shift_factor &&
                    ((shift_factor - min_shift_factor) < min_shift_factor_diff ||
                    byte_size(current_pos) < byte_size_lower_bound)){
                //PRINT_IF_NOT_SILENT(" -> ok"<<std::endl);
                found_best_factor = true;
            }
            i++;
        }

        //PRINT_IF_NOT_SILENT(" final min shift factor = "<<min_shift_factor<<std::endl);

        //re-compute best valid position
        current_pos = new_vertex_initial_position;
        PEHelpers::lower_precision(min_shift_factor, current_pos);

        //PRINT_IF_NOT_SILENT(" final byte size = "<<byte_size(current_pos)<<std::endl);

        this->set_vertex(tip_vertex, new_vertex_initial_position);

        //PRINT_IF_NOT_SILENT(" ...done, reduced precision from "<<initial_byte_size<<" bytes to "<<byte_size(current_pos)<<" after "<<i<<" iterations"<<std::endl);
        //PRINT_IF_NOT_SILENT(" ======================================================="<<std::endl);
        return current_pos;
#endif
    }







    CellHandle Expander::add_cell_to_cone(const CellHandle& ch,
                                          ExpansionCone& cone) const{

        auto cell_vertices = mesh_.get_cell_vertices(ch);
        int already_added_vertices_count(0);
        //PRINT_IF_NOT_SILENT(" adding cell "<<ch<<" : "<<cell_vertices<<" to cone: "<<cone<<std::endl);
        for(auto v: cell_vertices){

            //PRINT_IF_NOT_SILENT(" - checking vertex "<<v<<std::endl);
            auto cone_vertex = cone.mesh_to_cone_handle(v);
            if(cone_vertex.idx() == -1){


                cone_vertex = cone.add_vertex(v, this->vertex(v));
                if(cone_vertex.idx() == -1){
                    PRINT_IF_NOT_SILENT(" ERROR - couldn't add vertex "<<v<<" to expansion cone"<<std::endl);
                    return CellHandle(-1);
                }
                //PRINT_IF_NOT_SILENT(" -- added vertex "<<v<<" to the expansion cone"<<std::endl);
            }else{
                already_added_vertices_count++;
                //PRINT_IF_NOT_SILENT(" -- vertex "<<v<<" was already added to the cone"<<std::endl);
            }
        }

        CellHandle cone_cell(-1);

        //if(already_added_vertices_count != 4) {

            cone_cell = cone.add_cell(ch, cell_vertices);
            if (cone_cell.idx() == -1) {
                std::cout << " ERROR - couldn't add cell " << ch << " : " << cell_vertices << " to expansion cone"
                          << std::endl;
                return CellHandle(-1);
            }
        /*}else{
            cone_cell = cone.mesh_to_cone_handle(ch);
        }*/
        return cone_cell;
    }


    FaceHandle Expander::add_face_to_cone(const FaceHandle& fh,
                                          ExpansionCone& cone) const{

        auto face_vertices = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh, 0));
        //PRINT_IF_NOT_SILENT(" adding face "<<fh<<" : "<<face_vertices<<" to the expansion cone"<<std::endl);
        for(auto v: face_vertices){

            auto cone_vertex = cone.mesh_to_cone_handle(v);
            if(cone_vertex.idx() == -1){

                cone_vertex = cone.add_vertex(v, this->vertex(v));

                if(cone_vertex.idx() == -1){
                    PRINT_IF_NOT_SILENT(" ERROR - couldn't add vertex "<<v<<" to expansion cone"<<std::endl);
                    return FaceHandle(-1);
                }
                //PRINT_IF_NOT_SILENT(" -- added vertex "<<v<<" to the expansion cone"<<std::endl);
            }else{
                //PRINT_IF_NOT_SILENT(" -- vertex "<<v<<" was already added to the cone"<<std::endl);
            }
        }

        auto cone_face = cone.add_face(face_vertices);

        if(cone_face.idx() == -1){
            PRINT_IF_NOT_SILENT(" ERROR - couldn't add face "<<fh<<": "<<face_vertices<<" to the expansion cone"<<std::endl);
            return FaceHandle(-1);
        }

        return cone_face;
    }



    int Expander::expanded_vertices_count(const CellHandle& ch) const{
        int count(0);
        for(auto cv_it = mesh_.cv_iter(ch); cv_it.valid(); cv_it++){
            count += expanded_prop_[*cv_it];
        }
        return count;
    }

    int Expander::expanded_vertices_count(const FaceHandle& fh) const{
        int count(0);
        for(auto fv_it = mesh_.fv_iter(fh); fv_it.valid(); fv_it++){
            count += expanded_prop_[*fv_it];
        }
        return count;
    }


    bool Expander::all_unexpanded_vertices_are_at_the_center(const VertexHandle& center_vertex) const{

        for(auto v: mesh_.vertices()){
            if(v != center_vertex &&
               !expanded_prop_[v] &&
               this->vertex(v) != VertexPosition(0,0,0)){

                PRINT_IF_NOT_SILENT(" --> vertex "<<v<<" is unexpanded but located at "<<this->vertex(v)<<std::endl);
                return false;

            }
        }
        return true;
    }



    int Expander::expansion_valence(const VertexHandle& center_vertex) const{

        int expanded_valence(0);
        for(auto out_he_it: mesh_.outgoing_halfedges(center_vertex)){
            expanded_valence += expanded_prop_[mesh_.to_vertex_handle(out_he_it)];
        }
        return expanded_valence;
    }

    int Expander::unexpansion_valence(const VertexHandle& center_vertex) const{

        return (int)mesh_.valence(center_vertex) - expansion_valence(center_vertex);
    }



    VertexPosition Expander::vertex(const VertexHandle& v) const{

        return vertex_position_prop_[v];
    }


    void Expander::set_vertex(const VertexHandle& v,
                              const VertexPosition& pos){
        vertex_position_prop_[v] = pos;
        mesh_.set_vertex(v, vec2vec(pos));

        max_vertex_precision_B_ = std::max(max_vertex_precision_B_,
                                         max_byte_size(pos));

        data_logger_.add_precision_point(max_byte_size(pos));
    }


    void Expander::set_expanded_prop(const VertexHandle& v,
                                     const bool status){
        expanded_prop_[v] = status;
    }

    bool Expander::expanded_prop(const VertexHandle& v) const {
        return expanded_prop_[v];
    }

    VertexHandle Expander::split_edge(const EdgeHandle& e){

        auto from_v_pos = vertex(mesh_.edge(e).from_vertex());
        auto to_v_pos   = vertex(mesh_.edge(e).to_vertex());
        auto mid_v = mesh_.split_edge(e);

        domain_vertex_position_prop_[mid_v] = 0.5 * (domain_vertex_position_prop_[mesh_.edge(e).from_vertex()] + domain_vertex_position_prop_[mesh_.edge(e).to_vertex()]);

        return mid_v;
    }

    VertexHandle Expander::split_edge(const HalfEdgeHandle& e){
        return split_edge(mesh_.edge_handle(e));
    }


    void Expander::check_for_timeout() const{

        auto current_time = std::chrono::high_resolution_clock::now();
        int current_duration_s = (float)std::chrono::duration_cast<std::chrono::milliseconds>(current_time - global_start_time_).count() /1000;

        if(current_duration_s > max_allocated_time_s_){
            throw TimeOutException();
        }

    }

    int Expander::remaining_seconds_before_timeout() const{
        auto current_time = std::chrono::high_resolution_clock::now();
        int current_duration_s = (float)std::chrono::duration_cast<std::chrono::milliseconds>(current_time - global_start_time_).count() /1000;
        return max_allocated_time_s_ - current_duration_s;
    }


    int Expander::iteration_count() const{
        return iteration_count_;
    }

}
