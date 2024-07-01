#include "SphereMapper.hh"


namespace OpenVolumeMesh{




void SphereMapper::perturbate_vertices_positions(TetrahedralMesh& mesh){

    const double lambda(0.05);

    for(auto v: mesh.vertices()){
        if(mesh.is_boundary(v)){

            auto vertex_pos = mesh.vertex(v);

            Vec3d perturbation = {0,0,0};

            VertexHandle prev_neighbor(-1);

            for(auto out_he_it: mesh.outgoing_halfedges(v)){
                if(mesh.is_boundary(out_he_it)){
                    auto current_neighbor = mesh.to_vertex_handle(out_he_it);

                    if(prev_neighbor.idx() != -1 &&
                            prev_neighbor != current_neighbor){
                        auto prev_neighbor_pos = mesh.vertex(prev_neighbor);
                        auto cur_neighbor_pos = mesh.vertex(current_neighbor);

                        /*std::cout<<" prev v : "<<prev_neighbor<<std::endl;
                    std::cout<<" cur v: "<<current_neighbor<<std::endl;
                    std::cout<<" center v: "<<v<<std::endl;
                    std::cout<<" prev pos: "<<prev_neighbor_pos<<std::endl;
                    std::cout<<" cur pos : "<<cur_neighbor_pos<<std::endl;
                    std::cout<<" center pos: "<<vertex_pos<<std::endl;*/

                        //add small perturbation based on cross prod of surrounding vertices

                        auto perturbation_direction = (prev_neighbor_pos-vertex_pos).cross(cur_neighbor_pos-vertex_pos);

                        if(perturbation_direction.norm() == perturbation_direction.norm()){
                            perturbation += lambda * perturbation_direction.normalize();
                        }

                        if(perturbation.norm() != perturbation.norm()){
                            std::cout<<" nan perturbation... direction: "<<perturbation_direction<<std::endl;
                            return;
                        }
                    }

                    prev_neighbor = current_neighbor;
                }
            }
            /*std::cout<<" vertex "<<v<<" position : "<<vertex_pos<<std::endl;
        std::cout<<" perturbation : "<<perturbation<<std::endl;
        std::cout<<" updated position : "<<perturbation + vertex_pos<<std::endl;
        std::cout<<" ------"<<std::endl;*/
            mesh.set_vertex(v, perturbation + vertex_pos);
        }
    }
}


bool SphereMapper::project_boundary_to_sphere(TetrahedralMesh& mesh,
                                              double sphere_radius){

    //std::cout<<"-------------------------- "<<std::endl;
    //std::cout<<" projecting mesh boundary onto a sphere of radius "<<sphere_radius<<std::endl;

    for(auto v: mesh.vertices()){
        if(mesh.is_boundary(v)){
            auto v_pos = mesh.vertex(v);
            //std::cout<<" moved vertex "<<v<<" at "<<std::setprecision(20)<<mesh.vertex(v);
            mesh.set_vertex(v, v_pos * (sphere_radius / v_pos.norm()));
            //std::cout<<"  to "<<std::setprecision(20)<<mesh.vertex(v)<<std::endl;

        }
    }


    //std::cout<<"...done"<<std::endl;
    //std::cout<<"-------------------------- "<<std::endl;

    return true;
}


bool SphereMapper::is_boundary_star_shaped_and_non_degenerate(TetrahedralMesh mesh){

    auto bad_tets = BadTetFinder::findBadTets(mesh);
    std::cout<<" - deg tets: "<<bad_tets.first.size()<<std::endl;
    std::cout<<" - flipped tets: "<<bad_tets.second.size()<<std::endl;

    auto boundary_deg_triangles_count = PEHelpers::degenerate_boundary_faces_count(mesh);

    if(boundary_deg_triangles_count){
        std::cout<<" WARNING - degenerate boundary triangles found"<<std::endl;
    }

    return bad_tets.second.empty() && !boundary_deg_triangles_count;

}


Vec3d SphereMapper::boundary_centroid(const TetrahedralMesh& mesh){
    Vec3d centroid(0,0,0);
    double v_count(0);

    for(auto v: mesh.vertices()){
        if(mesh.is_boundary(v)){
            centroid += mesh.vertex(v);
            v_count++;
        }
    }

    centroid /= v_count;

    return centroid;
}



bool SphereMapper::center_mesh(TetrahedralMesh& mesh,
                               const TetrahedralMesh::PointT& target_centroid){

    //std::cout<<"-------------------------- "<<std::endl;
    //std::cout<<" centering mesh..."<<std::endl;

    auto centroid = boundary_centroid(mesh);

    //std::cout<<" - centroid at "<<centroid<<std::endl;

    for(auto v: mesh.vertices()){
        //std::cout<<" moved vertex "<<v<<" at "<<std::setprecision(20)<<mesh.vertex(v);
        mesh.set_vertex(v, mesh.vertex(v) - centroid + target_centroid);
        //std::cout<<"  to "<<std::setprecision(20)<<mesh.vertex(v)<<std::endl;
    }

    //std::cout<<"...done"<<std::endl;
    //std::cout<<"-------------------------- "<<std::endl;

    return true;
}


double SphereMapper::encompassing_sphere_radius(TetrahedralMesh& mesh){

    double max_distance(0);

    for(auto v: mesh.vertices()){
        double distance = mesh.vertex(v).norm();
        max_distance = std::max(distance, max_distance);
    }
    return max_distance;

}


void SphereMapper::scale_mesh_to(TetrahedralMesh& mesh,
                                 double new_scale){

    double max_distance(0);

    for(auto v: mesh.vertices()){
        double distance = mesh.vertex(v).norm();
        max_distance = std::max(distance, max_distance);
    }

    for(auto v: mesh.vertices()){
        mesh.set_vertex(v, mesh.vertex(v) * (max_distance / new_scale));
    }


}

void SphereMapper::move_centroid_and_scale_mesh_boundary_surface(TetrahedralMesh& mesh,
                                                                 const TetrahedralMesh::PointT& target_centroid,
                                                                 double target_surface_area,
                                                                 bool move_and_scale_interior){

    auto centroid = boundary_centroid(mesh);


    //std::cout<<" centroid = "<<centroid<<std::endl;

    // restore original surface area
    double area = boundary_area(mesh);

    //std::cout<<" boundary area = "<<area<<std::endl;
    double scale = sqrt(target_surface_area / area);
    for (auto v : mesh.vertices()){
        if(move_and_scale_interior || mesh.is_boundary(v)){

            auto new_pos = (mesh.vertex(v) - centroid) * scale + target_centroid ;

            mesh.set_vertex(v, new_pos);
        }
    }

    //std::cout<<" new centroid = "<<boundary_centroid(mesh)<<std::endl;
    //std::cout<<" new area = "<<boundary_area(mesh)<<std::endl;

}



void SphereMapper::map_interior_to_centroid(TetrahedralMesh& mesh,
                                            const TetrahedralMesh::PointT& target_centroid){

    for(auto v: mesh.vertices()){
        if(!mesh.is_boundary(v)){
            //std::cout<<" - moved "<<v<<" from "<<mesh.vertex(v)<<" with norm = "<<mesh.vertex(v).norm()<<" to "<<target_centroid<<std::endl;
            mesh.set_vertex(v, target_centroid);
        }
    }
}


double SphereMapper::boundary_area(const TetrahedralMesh& mesh){

    double area(0);

    for(auto f: mesh.faces()){
        if(mesh.is_boundary(f)){
            area += triangle_face_area(mesh, f);
        }
    }

    return area;
}



double SphereMapper::triangle_face_area(const TetrahedralMesh& mesh,
                                        const FaceHandle& face){

    auto vertices = mesh.get_halfface_vertices(mesh.halfface_handle(face, 0));

    if(vertices.size() != 3){
        return 0.;
    }

    auto A = mesh.vertex(vertices[0]);
    auto B = mesh.vertex(vertices[1]);
    auto C = mesh.vertex(vertices[2]);


    return 0.5 * ((B-A).cross(C-A)).norm();
}




double clamp_cotan(const double x){
    const double bound = 19.1; // 3 degrees
    return (x < -bound ? -bound : (x > bound ? bound : x));
}




}
