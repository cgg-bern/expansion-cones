#pragma once

#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace OpenVolumeMesh{

struct ExpansionStatus{
    int left_to_expand_count;
    int expanded_count;
    int vertices_count;
    int ec_size;
    int steiner_vertices_created;
    int available_single_steiner_vertex_ecs;
    int position_precision;
    int unexp_valence;
    int DBCI_vertices_count;
    int expanded_vertex_index;
};


class ExpansionDataLogger{

public:

    ExpansionDataLogger(const std::string& mesh_name);

    static constexpr int max_histogram_expansion_valence = 20;

    void add_histogram_point(int expansion_valence);

    void add_evolution_point(const ExpansionStatus status);

    void add_precision_point(const int position_precision);

    void add_cluster_computation_point(const int time_mus,
                                       const int max_prec,
                                       const int total_prec,
                                       const int n_vertices);

    void add_smoothing_computation_point(const int time_mus,
                                         const int max_val,
                                         const int total_val,
                                         const int n_smoothed_vertices);

    void add_iteration_timing_point(const int time_mus);


    /* FILE DATA FORMAT (one data element per line)
     * mesh name (given in ctor)
     * expansion result (0 for succes, -1 for error, >0 for failure
     * total time to expand (in seconds)
     * vertices to expand evolution (one value per iteration, separeted by a ws)
     * vertices expanded evolution (same as above)
     * total vertices count evolution (same as above)
     * histogram values, separated by a ws */
    bool write_data_to_file(const std::string& output_file_name,
                            const int result,
                            const int expansion_time_s,
                            const int max_precision_B,
                            const int total_precision_B,
                            const double average_precision_B,
                            const int total_saved_precision_B) const;


    void print_out_expansion_histogram() const;

    std::string mesh_name() const{
        return mesh_name_;
    }

private:

    std::string mesh_name_;

    std::vector<int> expansion_valence_histo_;

    std::vector<ExpansionStatus> expansion_evolution_;

    //first is time in microseconds, second is precision
    std::vector<std::pair<int, int>> position_precision_record_;

    //for each cluster: [computation time (mus), max precision, total precision, #vertices]
    std::vector<std::tuple<int, int, int, int>> cluster_computation_record_;

    //for each smoothing iteration: [computation time (mus), max valence, total valence, #smoothed vertices]
    std::vector<std::tuple<int, int, int, int>> smoothing_computation_record_;

    //microseconds for each iteration
    std::vector<int> iteration_time_record_;

    const std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;

};

}
