#pragma once

#include "BadTetFinder.hh"

namespace OpenVolumeMesh{


/** NOTE: precision is the sum of all three components' precision */
void compute_precision_stats(const TetrahedralMesh& mesh,
                             const VertexPropertyT<ExactVertexPosition>& exact_vertex_positions,
                             int& max_precision_B,
                             int& total_precision_B,
                             double& average_precision_B,
                             VertexHandle& max_precision_vertex);

    int byte_size(const ExactVertexPosition& pos);

    /* Returns the maximum of the minimum between the
     * numerator and denominator size in bytes of each element*/
    int max_min_size(const ExactVertexPosition& pos);

    int max_byte_size(const ExactVertexPosition& pos);
};


