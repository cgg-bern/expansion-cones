#ifndef MESHTOPOLOGICALCOMPARATOR_HH
#define MESHTOPOLOGICALCOMPARATOR_HH

#include <vector>
#include <algorithm>
#include "SubTriangleMap.hh"
#include "ProgEmbeddingHelpers.hh"


namespace OpenVolumeMesh{


bool is_subdivided_topologically_equivalent(const TetrahedralMesh& reference_mesh,
                                            const TetrahedralMesh& subdivided_mesh,
                                            const SubTriangleMap& sub_triangle_map);



/** NOTE: assumes that both meshes' vertices are indexed the same way.
 * i.e. if two meshes are topologically equivalent BUT their vertices are not the same
 * then it will return false */
template<class TetrahedralMesh>
bool is_topologically_equivalent(const TetrahedralMesh& mesh_one, const TetrahedralMesh& mesh_two){

    if(mesh_one.n_logical_vertices() != mesh_two.n_logical_vertices() ||
            mesh_one.n_logical_edges() != mesh_two.n_logical_edges() ||
            mesh_one.n_logical_faces() != mesh_two.n_logical_faces() ||
            mesh_one.n_logical_cells() != mesh_two.n_logical_cells()){

        std::cout<<" MESHES ARE NOT TOPOLOGICALLY EQUIVALENT: "<<std::endl;
        std::cout<<" -- "<<mesh_one.n_logical_vertices()<<" vs. "<<mesh_two.n_logical_vertices()<<" vertices"<<std::endl;
        std::cout<<" -- "<<mesh_one.n_logical_edges()<<" vs. "<<mesh_two.n_logical_edges()<<" vertices"<<std::endl;
        std::cout<<" -- "<<mesh_one.n_logical_faces()<<" vs. "<<mesh_two.n_logical_faces()<<" vertices"<<std::endl;
        std::cout<<" -- "<<mesh_one.n_logical_cells()<<" vs. "<<mesh_two.n_logical_cells()<<" vertices"<<std::endl;

        return false;
    }

    //check that all cells exist in both meshes
    for(auto c: mesh_one.cells()){
        std::vector<OpenVolumeMesh::VertexHandle> cell_vertices;
        for(auto cv_it = mesh_one.cv_iter(c); cv_it.valid(); cv_it++){
            cell_vertices.push_back(*cv_it);
        }

        bool found_cell_in_mesh_two(false);
        //find the cell in mesh two
        for(auto vc_it = mesh_two.vc_iter(cell_vertices[0]); vc_it.valid(); vc_it++){
            if(found_cell_in_mesh_two){
                break;
            }
            //check all vertices
            int found_vertices(0);
            for(auto cv_it = mesh_two.cv_iter(*vc_it); cv_it.valid(); cv_it++){
                found_vertices += (std::find(cell_vertices.begin(), cell_vertices.end(), *cv_it) != cell_vertices.end());
            }
            if(found_vertices == 4){
                found_cell_in_mesh_two = true;
            }
        }

        if(!found_cell_in_mesh_two){
            return false;
        }
    }

    return true;
}

}

#endif // MESHTOPOLOGICALCOMPARATOR_HH
