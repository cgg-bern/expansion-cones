#include "PrecisionDiagnostics.hh"


namespace OpenVolumeMesh{


void compute_precision_stats(const TetrahedralMesh& mesh,
                             const VertexPropertyT<ExactVertexPosition>& exact_vertex_positions,
                             int& max_precision_B,
                             int& total_precision_B,
                             double& average_precision_B,
                             VertexHandle& max_precision_vertex){

    max_precision_B     = 0;
    total_precision_B   = 0;
    average_precision_B = 0;

    double n(0);
    for(auto v: mesh.vertices()){
        auto precision_B = byte_size(exact_vertex_positions[v]);
        if(precision_B > max_precision_B){
            max_precision_B = precision_B;
            max_precision_vertex = v;
        }
        average_precision_B += precision_B;
        total_precision_B   += precision_B;
        n++;
    }
    average_precision_B /= n;
}


int byte_size(const ExactVertexPosition& pos){
    return pos[0].size() + pos[1].size() + pos[2].size();
}

int max_byte_size(const ExactVertexPosition& pos){
    return std::max(std::max(pos[0].size(), pos[1].size()), pos[2].size());
}


int max_min_size(const ExactVertexPosition& pos){

    auto x_min_size = std::min(pos[0].numerator().bit_size(), pos[0].denominator().bit_size());
    auto y_min_size = std::min(pos[1].numerator().bit_size(), pos[1].denominator().bit_size());
    auto z_min_size = std::min(pos[2].numerator().bit_size(), pos[2].denominator().bit_size());

    return std::max(std::max(x_min_size, y_min_size), z_min_size) / 8;
}
};
