#include "BadTetFinder.hh"

namespace OpenVolumeMesh{




/** Simple BadTetFinder */
BadTetFinder::BadTetFinder(const TetrahedralMesh& mesh) : mesh_(mesh){


}



CGAL_Point3 OVMvec3ToCGALPoint3(const TetrahedralMesh::PointT& vec){
    return {vec[0], vec[1], vec[2]};
}



CGAL_Tetrahedron OVMtetToCGALtet(const TetrahedralMesh& mesh,
                                 const CellHandle& ch){

    auto cell_vertices = mesh.get_cell_vertices(ch);
    return CGAL_Tetrahedron(
                OVMvec3ToCGALPoint3(mesh.vertex(cell_vertices[0])),
            OVMvec3ToCGALPoint3(mesh.vertex(cell_vertices[1])),
            OVMvec3ToCGALPoint3(mesh.vertex(cell_vertices[2])),
            OVMvec3ToCGALPoint3(mesh.vertex(cell_vertices[3])));
}



bool BadTetFinder::isFlipped(OpenVolumeMesh::CellHandle cell){

    std::vector<VertexHandle> cell_vertices = mesh_.get_cell_vertices(cell);

    if(cell_vertices.size() != 4){
        std::cout<<" ERROR: Cell "<<cell<<": ("<<cell_vertices<<") is not a tet, it has "<<cell_vertices.size()<<" vertices."<<std::endl;
        return true;
    }

    auto tet = OVMtetToCGALtet(mesh_,
                               cell);

    return tet.orientation() == CGAL::NEGATIVE;
    return tet.volume() < 0;

}



bool BadTetFinder::isDegenerate(OpenVolumeMesh::CellHandle cell) const{
    std::vector<VertexHandle> cell_vertices = mesh_.get_cell_vertices(cell);

    if(cell_vertices.size() != 4){
        std::cout<<" ERROR: Cell "<<cell<<": ("<<cell_vertices<<") is not a tet, it has "<<cell_vertices.size()<<" vertices."<<std::endl;
        return true;
    }

    return OVMtetToCGALtet(mesh_, cell).is_degenerate();
}



BadTetFinder::BadTetList BadTetFinder::findBadTets(const TetrahedralMesh& mesh){


    BadTetFinder bad_tet_finder(mesh);

    BadTetFinder::BadTetList bad_tets;

    for(auto c_it: mesh.cells()){

        if(bad_tet_finder.isDegenerate(c_it)){
            bad_tets.first.push_back(c_it);
        }else if(bad_tet_finder.isFlipped(c_it)){
            bad_tets.second.push_back(c_it);
        }
    }
    return bad_tets;
}



bool BadTetFinder::meshContainsFlippedTets(const TetrahedralMesh& mesh){

    BadTetFinder bad_tet_finder(mesh);

    for(auto c_it: mesh.cells()){
        if(bad_tet_finder.isFlipped(c_it)){
            return true;
        }
    }
    return false;
}


double BadTetFinder::total_signed_volume(const TetrahedralMesh& mesh){

    BadTetFinder bad_tet_finder(mesh);

    double vol_sum(0);
    for(auto c_it: mesh.cells()){
        auto vol = OVMtetToCGALtet(mesh, c_it).volume();
        vol_sum += vol;
    }

    return vol_sum;
}


double BadTetFinder::total_unsigned_volume(const TetrahedralMesh& mesh){

    BadTetFinder bad_tet_finder(mesh);

    double vol_sum(0);
    for(auto c_it: mesh.cells()){
        auto vol = OVMtetToCGALtet(mesh, c_it).volume();
        vol_sum += abs(vol);
    }

    return vol_sum;
}



double BadTetFinder::total_negative_volume(const TetrahedralMesh& mesh,
                                           const CellPropertyT<bool>* skip_cells_prop){

    BadTetFinder bad_tet_finder(mesh);

    double vol_sum(0);
    for(auto c_it: mesh.cells()){
        if(skip_cells_prop && (*skip_cells_prop)[c_it]){
            continue;
        }
        auto vol = OVMtetToCGALtet(mesh, c_it).volume();
        vol_sum += (vol < 0 ? vol : 0);
    }

    return vol_sum;
}



double BadTetFinder::total_negative_volume_in_1_ring(const TetrahedralMesh& mesh,
                                                     const std::set<VertexHandle>& vs,
                                                     const CellPropertyT<bool>* skip_cells_prop){


    std::vector<bool> visited_prop(mesh.n_cells());

    double vol_sum(0);

    for(auto v: vs){
        for(auto vc_it = mesh.vc_iter(v); vc_it.valid(); vc_it++){
            if(skip_cells_prop && (*skip_cells_prop)[*vc_it]){
                continue;
            }

            if(!visited_prop[vc_it->idx()]){

                auto vol = OVMtetToCGALtet(mesh, *vc_it).volume();
                vol_sum += (vol < 0 ? vol : 0);
                visited_prop[vc_it->idx()] = true;
            }

        }
    }

    return vol_sum;

}





/* Exact BadTetFinder */


#define DISABLE_ALL_CHECKS 0

#if DISABLE_ALL_CHECKS
#define RETURN_FALSE_IF_CHECKS_DISABLED return false;
#define RETURN_EMPTY_STRUCT_IF_CHECKS_DISABLED return BadTetList();
#else
#define RETURN_FALSE_IF_CHECKS_DISABLED
#define RETURN_EMPTY_STRUCT_IF_CHECKS_DISABLED
#endif



ExactBadTetFinder::ExactBadTetFinder(const TetrahedralMesh& mesh,
                                     const VertexPropertyT<ExactVertexPosition>& exact_positions)
    : BadTetFinder(mesh),
      exact_positions_(exact_positions){
}

ExactVertexPosition OVMvec3ToCGALPoint3(const CGAL_ExactPoint3& vec,
                                        const int scale_factor){
    return {vec[0] * scale_factor, vec[1] * scale_factor, vec[2] * scale_factor};
}

CGAL_ExactPoint3 OVMvec3ToCGALPoint3(const ExactVertexPosition& vec,
                                     const int scale_factor){
    return {vec[0] * scale_factor, vec[1] * scale_factor, vec[2] * scale_factor};
}


CGAL_ExactTetrahedron OVMtetToCGALtet(const TetrahedralMesh& mesh,
                                      const VertexPropertyT<ExactVertexPosition>& exact_positions,
                                      const CellHandle& ch){

    auto cell_vertices = mesh.get_cell_vertices(ch);
    return CGAL_ExactTetrahedron(
                OVMvec3ToCGALPoint3(exact_positions[cell_vertices[0]]),
            OVMvec3ToCGALPoint3(exact_positions[cell_vertices[1]]),
            OVMvec3ToCGALPoint3(exact_positions[cell_vertices[2]]),
            OVMvec3ToCGALPoint3(exact_positions[cell_vertices[3]]));
}



bool ExactBadTetFinder::isFlipped(OpenVolumeMesh::CellHandle cell){

    RETURN_FALSE_IF_CHECKS_DISABLED;

    std::vector<VertexHandle> cell_vertices = mesh_.get_cell_vertices(cell);

    if(cell_vertices.size() != 4){
        std::cout<<" ERROR: Cell "<<cell<<": ("<<cell_vertices<<") is not a tet, it has "<<cell_vertices.size()<<" vertices."<<std::endl;
        return true;
    }

    auto tet = OVMtetToCGALtet(mesh_,
                               exact_positions_,
                               cell);

    return tet.orientation() == CGAL::NEGATIVE;
    return tet.volume() < 0;

}



bool ExactBadTetFinder::isDegenerate(OpenVolumeMesh::CellHandle cell) const{

    RETURN_FALSE_IF_CHECKS_DISABLED;

    std::vector<VertexHandle> cell_vertices = mesh_.get_cell_vertices(cell);

    if(cell_vertices.size() != 4){
        std::cout<<" ERROR: Cell "<<cell<<": ("<<cell_vertices<<") is not a tet, it has "<<cell_vertices.size()<<" vertices."<<std::endl;
        return true;
    }

    return OVMtetToCGALtet(mesh_,
                           exact_positions_,
                           cell).is_degenerate();
}



BadTetFinder::BadTetList ExactBadTetFinder::findBadTets(const TetrahedralMesh& mesh,
                                                        const VertexPropertyT<ExactVertexPosition>& exact_positions){

    RETURN_EMPTY_STRUCT_IF_CHECKS_DISABLED

    ExactBadTetFinder bad_tet_finder(mesh,
                                     exact_positions);

    BadTetFinder::BadTetList bad_tets;

    for(auto c_it: mesh.cells()){

        if(bad_tet_finder.isDegenerate(c_it)){
            bad_tets.first.push_back(c_it);
        }else if(bad_tet_finder.isFlipped(c_it)){
            bad_tets.second.push_back(c_it);
        }
    }
    return bad_tets;
}



bool ExactBadTetFinder::meshContainsFlippedTets(const TetrahedralMesh& mesh,
                                                const VertexPropertyT<ExactVertexPosition>& exact_positions){


    RETURN_FALSE_IF_CHECKS_DISABLED;

    ExactBadTetFinder bad_tet_finder(mesh,
                                     exact_positions);

    for(auto c_it: mesh.cells()){
        if(bad_tet_finder.isFlipped(c_it)){
            /*std::cout<<" --> found tet with volume = "<<OVMtetToCGALtet(mesh,
                                                                        c_it,
                                                                        exact_positions).volume()<<std::endl;*/

            //std::cout<<" --> in double precision: "<<OVMtetToCGALtet(mesh,c_it).volume()<<std::endl;
            return true;
        }
    }
    return false;
}


bool ExactBadTetFinder::meshContainsBadTets(const TetrahedralMesh& mesh,
                                            const VertexPropertyT<ExactVertexPosition>& exact_positions){


    RETURN_FALSE_IF_CHECKS_DISABLED;

    ExactBadTetFinder bad_tet_finder(mesh,
                                     exact_positions);

    std::cout<<" lakjsdflkajs"<<std::endl;

    for(auto c_it: mesh.cells()){
        std::cout<<" - checking cell "<<c_it<<", deleted = "<<mesh.is_deleted(c_it)<<std::endl;
        if(mesh.is_deleted(c_it)){
            continue;
        }
        std::cout<<" - checking cell "<<c_it<<": "<<mesh.get_cell_vertices(c_it)<<std::endl;
        if(bad_tet_finder.isFlipped(c_it) || bad_tet_finder.isDegenerate(c_it)){
            std::cout<<" --> found tet with volume = "<<OVMtetToCGALtet(mesh,
                                                                        exact_positions,
                                                                        c_it).volume()<<std::endl;

            /*std::cout<<" --> flipped! In double precision: "<<OVMtetToCGALtet(mesh,
                                                                     c_it).volume()<<std::endl;*/
            return true;
        }
    }
    return false;
}


bool ExactBadTetFinder::meshContainsDegenerateTets(const TetrahedralMesh& mesh,
                                                   const VertexPropertyT<ExactVertexPosition>& exact_positions){


    RETURN_FALSE_IF_CHECKS_DISABLED;

    ExactBadTetFinder bad_tet_finder(mesh,
                                     exact_positions);

    for(auto c_it: mesh.cells()){
        if(bad_tet_finder.isDegenerate(c_it)){
            /*std::cout<<" --> found tet "<<mesh.get_cell_vertices(c_it)<<
                       " with volume = "<<OVMtetToCGALtet(mesh,
                                                          c_it,
                                                          exact_positions).volume()<<std::endl;

            std::cout<<" --> in double precision: "<<BadTetFinder::OVMtetToCGALtet(mesh,
                                                                                   c_it).volume()<<std::endl;
            */
            return true;
        }
    }
    return false;
}




bool ExactBadTetFinder::meshContainsFlippedTetsIn1Ring(const TetrahedralMesh& mesh,
                                                       const VertexPropertyT<ExactVertexPosition>& exact_positions,
                                                       const VertexHandle& vh){

    ExactBadTetFinder bad_tet_finder(mesh,
                                     exact_positions);

    for(auto vc_it = mesh.vc_iter(vh); vc_it.valid(); vc_it++){
        if(bad_tet_finder.isFlipped(*vc_it)){
            return true;
        }
    }
    return false;
}


bool ExactBadTetFinder::meshContainsDegenerateTetsIn1Ring(const TetrahedralMesh& mesh,
                                                          const VertexPropertyT<ExactVertexPosition>& exact_positions,
                                                          const VertexHandle& vh){
    ExactBadTetFinder bad_tet_finder(mesh,
                                     exact_positions);

    for(auto vc_it = mesh.vc_iter(vh); vc_it.valid(); vc_it++){
        if(bad_tet_finder.isDegenerate(*vc_it)){
            return true;
        }
    }
    return false;
}



std::ostream& operator<<(std::ostream& os,
                         const std::vector<VertexHandle>& vertices){

    for(auto v: vertices){
        os<<" "<<v;
    }
    return os;
}


std::ostream& operator<<(std::ostream& os,
                         const std::set<VertexHandle>& vertices){

    for(auto v: vertices){
        os<<" "<<v;
    }
    return os;
}


std::ostream& operator<<(std::ostream& os,
                         const std::list<OpenVolumeMesh::VertexHandle>& vertices){

    for(auto v: vertices){
        os<<" "<<v;
    }
    return os;
}
}
