#include "BadTriangleFinder.hh"


using namespace OpenMesh;


std::vector<OpenMesh::FaceHandle> BadTriangleFinder::find_flipped_triangles(TriMesh& mesh){

    std::vector<OpenMesh::FaceHandle> flipped_triangles;

    bool first_non_zero_tri(true);
    bool first_tri_has_positive_area(false);


    std::cout<<" checking area of triangles..."<<std::endl;
    for(auto tri: mesh.faces()){
        auto area = signed_triangle_area(mesh, tri);

        if(!area){
            continue;
        }
        bool positive_area = (area >= 0);

        if(area != 0 &&
                first_non_zero_tri){
            first_non_zero_tri = false;
            first_tri_has_positive_area = positive_area;
            std::cout<<" first non-zero tri "<<tri<<
                       " area = "<<area<<
                       ", first is positive: "<<first_tri_has_positive_area<<std::endl;
        }

        if(first_tri_has_positive_area != positive_area){
            flipped_triangles.push_back(tri);
        }

        std::cout<<" tri "<<tri<<": ";
        for(auto fv_it = mesh.fv_iter(tri); fv_it.is_valid(); fv_it++){
            std::cout<<*fv_it<<" ";
        }
        std::cout<<", area = "<<area<<", is positive: "<<positive_area<<std::endl;
    }

    return flipped_triangles;
}



double BadTriangleFinder::signed_triangle_area(TriMesh& mesh,
                                               const OpenMesh::FaceHandle& tri){

    std::vector<ACG::Vec3d> positions;
    for(auto fv_it = mesh.fv_iter(tri); fv_it.is_valid(); fv_it++){
        positions.push_back(mesh.point(*fv_it));
    }

    const auto& x = positions[0];
    const auto& y = positions[1];
    const auto& z = positions[2];

    return 0.5 * ((y[0]-x[0]) * (z[1]-y[1]) - (z[0]-y[0]) * (y[1]-x[1])) ;

}

