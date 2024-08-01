#include "BatchProcessor.hh"

using namespace OpenVolumeMesh;

using Mesh = TetrahedralMesh;



int BatchProcessor::generate_boundary_conditions(MeshProcessingFunction function,
                                                 std::string filename,
                                                 bool pre_process,
                                                 TetrahedralMesh& input_mesh,
                                                 TetrahedralMesh& result_mesh){

    Mesh mesh;
    //int mesh;
    OpenVolumeMesh::IO::FileManager fileManager;
    //std::cout<<" READING MESH "<<filename.toStdString()<<"..."<<std::endl;
    auto file_read_status = fileManager.readFile(filename, mesh);
    if(!file_read_status){
        std::cout << " couldn't load mesh " << filename << std::endl;
        return -1;
    }

    for(auto v: mesh.vertices()){
        if(std::isnan(mesh.vertex(v).norm())){
            std::cout<<" ERROR - found nan-norm vertex position in input mesh loaded from file"<<std::endl;
            return -1;
        }
    }
    for(auto v: mesh.vertices()){
        if(mesh.vertex(v).norm() > 1e20){
            std::cout<<" ERROR - found vertex position with norm > 1e20 in input mesh loaded from file"<<std::endl;
            return -1;
        }
    }

    for(auto e: mesh.edges()){
        int b_face(0);
        for(auto ef_it = mesh.ef_iter(e); ef_it.valid(); ef_it++){
            b_face += mesh.is_boundary(*ef_it);
        }
        if(b_face > 2){
            std::cout<<" ERROR - found non-manifold edge "<<mesh.edge(e)<<std::endl;
            return -1;
        }
    }

    if(pre_process){
        TetMeshBoundarySplitter::preProcessProblematicRegions(mesh);

        auto initial_bad_tets = BadTetFinder::findBadTets(mesh);

        if(initial_bad_tets.first.size()){
            std::cerr<<" WARNING - pre-processed mesh contains "<<initial_bad_tets.first.size()<<" degenerate tets"<<std::endl;
            //return -1;
        }
    }
    input_mesh = mesh;

    SphereMapper::center_mesh(mesh);

    auto result = function(mesh);
    result_mesh = mesh;

    return result;
}


#if 0
int BatchProcessor::shrink_and_expand(std::string mesh_file_path,
                                      bool pre_process,
                                      int option1,
                                      int option2,
                                      TetrahedralMesh& input_mesh,
                                      TetrahedralMesh& result_mesh,
                                      const std::string& mesh_name,
                                      const std::string& output_file_path){

    Mesh mesh;
    //int mesh;
    OpenVolumeMesh::IO::FileManager fileManager;

    //std::cout<<" READING MESH "<<filename.toStdString()<<"..."<<std::endl;
    auto file_read_status = fileManager.readFile(mesh_file_path, mesh);
    if(!file_read_status){
        std::cout<<" couldn't load mesh "<<mesh_file_path<<std::endl;
        return -1;
    }

    /*for(auto v: mesh.vertices()){
        if(std::isnan(mesh.vertex(v).norm())){
            std::cout<<" - "<<v<<" at "<<mesh.vertex(v)<<std::endl;
        }
    }
    std::cout<<" ---------------"<<std::endl;*/


    if(pre_process){

        TetMeshBoundarySplitter::preProcessProblematicRegions(mesh);

        auto initial_bad_tets = BadTetFinder::findBadTets(mesh);

        if(initial_bad_tets.first.size()){
            std::cerr<<" WARNING - pre-processed mesh contains "<<initial_bad_tets.first.size()<<" degenerate tets"<<std::endl;
            //return -1;
        }

    }
    input_mesh = mesh;

    auto result = function(mesh, mesh_name, output_file_path, option1, option2);
    result_mesh = mesh;

    return result;
}
#endif





