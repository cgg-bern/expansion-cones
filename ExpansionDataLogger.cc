#include "ExpansionDataLogger.hh"


namespace OpenVolumeMesh{


ExpansionDataLogger::ExpansionDataLogger(const std::string& mesh_name)
    : mesh_name_(mesh_name),
      expansion_valence_histo_(max_histogram_expansion_valence + 1),
      start_time_(std::chrono::high_resolution_clock::now())

{}

void ExpansionDataLogger::add_histogram_point(int expansion_valence){
    expansion_valence_histo_[std::min(std::max(expansion_valence, 0),
                                      max_histogram_expansion_valence)]++;
}


/*
void ExpansionDataLogger::add_total_timing(const int){
    auto start_time = std::chrono::high_resolution_clock::now();

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout<<std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()<<" SECONDS"<<std::endl;
}
*/
void ExpansionDataLogger::add_evolution_point(const ExpansionStatus status){

    /*std::cout<<" added evolution point "<<
               status.left_to_expand_count<<", "<<
               status.expanded_count<<", "<<
               status.vertices_count<<", "<<
               status.ec_size<<", "<<
               status.to_expand_count_increase<<std::endl;*/
    expansion_evolution_.push_back(status);

}

void ExpansionDataLogger::add_precision_point(const int position_precision){

    auto current_time = std::chrono::high_resolution_clock::now();
    int duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time_).count();
    position_precision_record_.push_back({duration_ms, position_precision});
}


void ExpansionDataLogger::add_cluster_computation_point(const int time_mus,
                                                        const int max_prec,
                                                        const int total_prec,
                                                        const int n_vertices){
    cluster_computation_record_.emplace_back(std::tuple<int, int, int, int>(time_mus, max_prec, total_prec, n_vertices));
}

void ExpansionDataLogger::add_smoothing_computation_point(const int time_mus,
                                                          const int max_val,
                                                          const int total_val,
                                                          const int n_smoothed_vertices){

    smoothing_computation_record_.emplace_back(std::tuple<int, int, int, int>(time_mus, max_val, total_val, n_smoothed_vertices));
}


void ExpansionDataLogger::add_iteration_timing_point(const int time_mus){
    iteration_time_record_.push_back(time_mus);
}


bool ExpansionDataLogger::write_data_to_file(const std::string& output_file_path,
                                             const int result,
                                             const int expansion_time_s,
                                             const int max_precision_B,
                                             const int total_precision_B,
                                             const double average_precision_B,
                                             const int total_saved_precision_B) const{

    //std::cout<<" expansion  result in log file: "<<result<<std::endl;
    //std::cout<<" expansion time(s) in log file: "<<expansion_time_s<<std::endl;

    /*if(output_file_path.empty()){
        return false;
    }*/

    if(output_file_path.empty()){
        std::cout<<" can't write data to empty file..."<<std::endl;
        return false;
    }

    std::ofstream output_data_file;
    output_data_file.open(output_file_path);

    if(!output_data_file.is_open()){
        std::cout<<" ERROR - couldn't open file "<<output_file_path<<std::endl;
        return false;
    }

    output_data_file<<mesh_name_<<std::endl;
    output_data_file<<result<<std::endl;
    output_data_file<<expansion_time_s<<std::endl;
    for(auto status: expansion_evolution_){
        output_data_file<<status.left_to_expand_count<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.expanded_count<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.vertices_count<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.steiner_vertices_created<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.available_single_steiner_vertex_ecs<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.position_precision<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.unexp_valence<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.DBCI_vertices_count<<" ";
    }
    output_data_file<<std::endl;

    for(auto status: expansion_evolution_){
        output_data_file<<status.expanded_vertex_index<<" ";
    }
    output_data_file<<std::endl;

    for(auto val: expansion_valence_histo_){
        output_data_file<<val<<" ";
    }
    output_data_file<<std::endl;

    output_data_file<<max_precision_B<<std::endl;
    output_data_file<<total_precision_B<<std::endl;
    output_data_file<<average_precision_B<<std::endl;
    output_data_file<<total_saved_precision_B<<std::endl;

    output_data_file.close();



    std::cout<<" WARNING: disabled additional diagnostics exportation"<<std::endl;

    /*std::ofstream precision_output_data_file;
    precision_output_data_file.open(output_file_path+".prec.txt");
    if(!precision_output_data_file.is_open()){
        std::cout<<" ERROR - couldn't open file "<<(output_file_path+".prec.txt")<<std::endl;
        return false;
    }
    for(auto prec_point: position_precision_record_){
        precision_output_data_file<<prec_point.first<<" ";
    }
    precision_output_data_file<<std::endl;
    for(auto prec_point: position_precision_record_){
        precision_output_data_file<<prec_point.second<<" ";
    }
    precision_output_data_file<<std::endl;
    precision_output_data_file.close();*/




    /*std::ofstream cluster_output_data_file;
    cluster_output_data_file.open(output_file_path+".cluster.txt");
    if(!cluster_output_data_file.is_open()){
        std::cout<<" ERROR - couldn't open file "<<(output_file_path+".cluster.txt")<<std::endl;
        return false;
    }
    for(auto cluster_point: cluster_computation_record_){
        cluster_output_data_file<<std::get<0>(cluster_point)<<
                                  " "<<std::get<1>(cluster_point)<<
                                  " "<<std::get<2>(cluster_point)<<
                                  " "<<std::get<3>(cluster_point)<<std::endl;
    }
    cluster_output_data_file.close();*/




    /*std::ofstream smoothing_output_data_file;
    smoothing_output_data_file.open(output_file_path+".smoothing.txt");
    if(!smoothing_output_data_file.is_open()){
        std::cout<<" ERROR - couldn't open file "<<(output_file_path+".smoothing.txt")<<std::endl;
        return false;
    }
    for(auto smoothing_point: smoothing_computation_record_){
        smoothing_output_data_file<<std::get<0>(smoothing_point)<<
                                  " "<<std::get<1>(smoothing_point)<<
                                  " "<<std::get<2>(smoothing_point)<<
                                  " "<<std::get<3>(smoothing_point)<<std::endl;
    }

    smoothing_output_data_file.close();*/




    std::ofstream iteration_time_output_file;
    iteration_time_output_file.open(output_file_path+".iterations.txt");
    if(!iteration_time_output_file.is_open()){
        std::cout<<" ERROR - couldn't open file "<<(output_file_path+".iterations.txt")<<std::endl;
        return false;
    }
    for(auto it_time: iteration_time_record_){
        iteration_time_output_file<<it_time<<" ";
    }
    iteration_time_output_file<<std::endl;
    iteration_time_output_file.close();

    return true;
}

void ExpansionDataLogger::print_out_expansion_histogram() const{

   for(size_t i(0); i<expansion_valence_histo_.size(); i++){
       std::cout<<" "<<std::setw(3)<<i<<" : "<<expansion_valence_histo_[i]<<std::endl;
   }
}

}
