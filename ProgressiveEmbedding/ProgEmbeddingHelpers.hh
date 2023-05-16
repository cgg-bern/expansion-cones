#pragma once

#include "TopoHelper.hh"
#include "OtherHelpers.hh"
#include "ExpansionCone.hh"
#include <Eigen/Dense>



namespace OpenVolumeMesh{



/* Progressive Embedding Helpers */
class PEHelpers{
public:

    /* This assumes that the mesh has two persistent properties called
     * str_exact_position and str_domain_exact_position
     * These should be string vertex props containing three rational numbers, separated by a ; */
    static void compute_and_output_quality_statistics(TetrahedralMesh& mesh, const std::string& output_file_path);

    static void compute_and_output_quality_statistics(const TetrahedralMesh& domain_mesh,
                                                      const TetrahedralMesh& codomain_mesh,
                                                      const std::string& output_file_path);

    template<typename PositionType, typename ScalarType>
    static void conformal_and_areal_components(const std::vector<PositionType>& xyz_coords,
                                               const std::vector<PositionType>& uvw_coords,
                                               const ScalarType& domain_volume,
                                               const ScalarType& codomain_volume,
                                               ScalarType& conformal,
                                               ScalarType& areal){

        auto& xyz = xyz_coords;
        auto& uvw = uvw_coords;

        Eigen::Matrix<ScalarType, 3, 3> U;
        Eigen::Matrix<ScalarType, 3, 3> X;

        auto u01 = uvw[1]-uvw[0];
        auto u02 = uvw[2]-uvw[0];
        auto u03 = uvw[3]-uvw[0];
        U<<u01[0], u02[0], u03[0],
           u01[1], u02[1], u03[1],
           u01[2], u02[2], u03[2];


        auto x01 = xyz[1]-xyz[0];
        auto x02 = xyz[2]-xyz[0];
        auto x03 = xyz[3]-xyz[0];
        X<<x01[0], x02[0], x03[0],
           x01[1], x02[1], x03[1],
           x01[2], x02[2], x03[2];

        /*std::cout<<" x01 = "<<x01<<std::endl;
        std::cout<<" x02 = "<<x02<<std::endl;
        std::cout<<" x03 = "<<x03<<std::endl;

        std::cout<<" x01 = "<<vec2vec(x01)<<std::endl;
        std::cout<<" x02 = "<<vec2vec(x02)<<std::endl;
        std::cout<<" x03 = "<<vec2vec(x03)<<std::endl;*/

        //std::cout<<" U = "<<U<<std::endl;
        //std::cout<<" X = "<<X<<std::endl;
        //auto X_inv = X.colPivHouseholderQr().solve(I);
        //J=U*X.inverse();

        auto X_inv = X.inverse();

        //std::cout<<" X-1 = "<<X_inv<<std::endl;

        auto J = U * X_inv;

        //std::cout<<" J = "<<J<<std::endl;

        auto det_J = J.determinant();

        //std::cout<<" |J| = "<<det_J<<std::endl;

        //auto volume_scaling = codomain_volume / domain_volume;

        double d_det_J = CGAL::to_double(det_J);
        if(d_det_J != 0){
            conformal = (J.transpose()*J).trace()/pow(d_det_J, 2.0/3.0);
        }else{
            conformal = std::numeric_limits<double>::max();
        }
        areal = (codomain_volume * det_J)/domain_volume + domain_volume/(codomain_volume*det_J);

        //std::cout<<" ---------- "<<std::endl;
    }




    static VertexPosition string_to_position(const std::string& str_exact_pos);



    /** DBCI: "Double-Boundary-Connected Interior"
     * NOTE: it can be connceted to more than two boundary vertices*/
    static int DBCI_vertices_count(const TetrahedralMesh& mesh);

    /** "Boundary-Connecting Interior" edges: edges whose endpoints are
     *  both boundary vertices but the edge itself is not a boundary edge
     **/
    static int BCI_edges_count(const TetrahedralMesh& mesh);

    static int degenerate_boundary_faces_count(const TetrahedralMesh& mesh);


    static size_t find_maximum_shift_factor(const ExactVertexPosition& new_vertex_position);


    /** uniformly makes a integer division on each of
     *  pos' elements' numerator and denominator
     *  using positional shift
     *
     *  i.e returns x'= (x0, x1, x2) where
     *      xi' = (p'/q') where
     *      p = p' >> factor + r_q,
     *      q = q' >> factor + r_p and
     *      xi = p/q and
     *      pos = (x0, x1, x2)   */
    static void lower_precision(unsigned long shift_factor,
                                ExactVertexPosition& pos);

    static int count_entangled_edges(const TetrahedralMesh& mesh);

};

}

