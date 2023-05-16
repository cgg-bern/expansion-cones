#pragma once

#include <CGAL/Tetrahedron_3.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>


#include "CommonMeshDefinitions.hh"


#define ENABLE_EXACT_REPRESENTATION 1

#if ENABLE_EXACT_REPRESENTATION
using ExactType = CGAL::Gmpq;
#else
using ExactType = double;
#endif

using ExactVertexPosition = OpenVolumeMesh::VectorT<ExactType, 3>;


//using CGAL_Kernel      = CGAL::Exact_predicates_exact_constructions_kernel;
using CGAL_Kernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGAL_Point3      = CGAL::Point_3<CGAL_Kernel>;
using CGAL_Tetrahedron = CGAL::Tetrahedron_3<CGAL_Kernel>;
using CGAL_Triangle3   = CGAL::Triangle_3<CGAL_Kernel>;


using CGAL_ExactKernel      = CGAL::Cartesian<ExactType>;
//using CGAL_ExactKernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGAL_ExactPoint3      = CGAL::Point_3<CGAL_ExactKernel>;
using CGAL_ExactTetrahedron = CGAL::Tetrahedron_3<CGAL_ExactKernel>;
using CGAL_ExactTriangle    = CGAL::Triangle_3<CGAL_ExactKernel>;



namespace OpenVolumeMesh{



CGAL_Point3 OVMvec3ToCGALPoint3(const TetrahedralMesh::PointT& vec);

CGAL_Tetrahedron OVMtetToCGALtet(const TetrahedralMesh& mesh,
                                 const CellHandle& ch);

ExactVertexPosition OVMvec3ToCGALPoint3(const CGAL_ExactPoint3& vec,
                                        const int scale_factor = 1);

/* scale_factor applied to all components */
CGAL_ExactPoint3 OVMvec3ToCGALPoint3(const ExactVertexPosition& vec,
                                     const int scale_factor = 1);

CGAL_ExactTetrahedron OVMtetToCGALtet(const TetrahedralMesh& mesh,
                                      const VertexPropertyT<ExactVertexPosition>& exact_positions,
                                      const CellHandle& ch);

class BadTetFinder
{
public:

    //.first = degenerate tets, .second = flipped tets
    using BadTetList = std::pair<
    std::vector<CellHandle>,
    std::vector<CellHandle>>;

    static BadTetList findBadTets(const TetrahedralMesh& mesh);

    static bool meshContainsFlippedTets(const TetrahedralMesh& mesh);

    //sum of all flipped tets' volume
    static double total_negative_volume(const TetrahedralMesh& mesh,
                                        const CellPropertyT<bool>* skip_cells_prop = nullptr);

    //sum of all tets' volume
    static double total_signed_volume(const TetrahedralMesh& mesh);

    static double total_unsigned_volume(const TetrahedralMesh& mesh);

    static double total_negative_volume_in_1_ring(const TetrahedralMesh& mesh,
                                                  const std::set<VertexHandle>& vs,
                                                  const CellPropertyT<bool>* skip_cells_prop = nullptr);


protected:

    BadTetFinder(const TetrahedralMesh& mesh);

    bool isFlipped(OpenVolumeMesh::CellHandle cell);

    bool isDegenerate(OpenVolumeMesh::CellHandle cell) const;

    const TetrahedralMesh& mesh_;
};





class ExactBadTetFinder : public BadTetFinder{

public:

    static BadTetList findBadTets(const TetrahedralMesh& mesh,
                                  const VertexPropertyT<ExactVertexPosition>& exact_positions);

    static bool meshContainsFlippedTets(const TetrahedralMesh& mesh,
                                        const VertexPropertyT<ExactVertexPosition>& exact_positions);


    static bool meshContainsDegenerateTets(const TetrahedralMesh& mesh,
                                           const VertexPropertyT<ExactVertexPosition>& exact_positions);


    static bool meshContainsBadTets(const TetrahedralMesh& mesh,
                                    const VertexPropertyT<ExactVertexPosition>& exact_positions);




    static bool meshContainsFlippedTetsIn1Ring(const TetrahedralMesh& mesh,
                                               const VertexPropertyT<ExactVertexPosition>& exact_positions,
                                               const VertexHandle& vh);


    static bool meshContainsDegenerateTetsIn1Ring(const TetrahedralMesh& mesh,
                                                  const VertexPropertyT<ExactVertexPosition>& exact_positions,
                                                  const VertexHandle& vh);




    ExactBadTetFinder(const TetrahedralMesh& mesh,
                      const VertexPropertyT<ExactVertexPosition>& exact_positions);

    bool isFlipped(OpenVolumeMesh::CellHandle cell);

    bool isDegenerate(OpenVolumeMesh::CellHandle cell) const;

private:

    const VertexPropertyT<ExactVertexPosition>& exact_positions_;
};


std::ostream& operator<<(std::ostream& os,
                         const std::vector<OpenVolumeMesh::VertexHandle>& vertices);


std::ostream& operator<<(std::ostream& os,
                         const std::set<OpenVolumeMesh::VertexHandle>& vertices);

std::ostream& operator<<(std::ostream& os,
                         const std::list<OpenVolumeMesh::VertexHandle>& vertices);

}
