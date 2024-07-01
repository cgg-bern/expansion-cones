#include "BatchProcessor.hh"
#include "ConnectedVertexSubsetIterator.hh"

#include "LightWeightStopWatch.hh"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#define DEFAULT_FUNCTION_INDEX -2

using namespace OpenVolumeMesh;

bool extract_mesh_name_and_directory(std::string location,
                                     std::string& filename,
                                     std::string& directory);


void writeResultMeshInSubdirectory(std::string location,
                                   std::string sub_directory,
                                   TetrahedralMesh& result_mesh,
                                   std::string suffix = "output");


void export_input_and_output_for_FOF_and_TLC_comparison(std::string location,
                                           std::string sub_directory,
                                           TetrahedralMesh& input,
                                           TetrahedralMesh& output);



void print_usage(){
    std::cerr<<" EXPECTED USAGE: ShrinkAndExpand [input_mesh] [output_location] [function_index] [boundary_mapping] [option]"<<std::endl;
    std::cout<<"   ---------------------------------------------------- "<<std::endl;
    std::cout<<"   [input_mesh] the input .ovm mesh file"<<std::endl;
    std::cout<<"   ---------------------------------------------------- "<<std::endl;
    std::cout<<"   [output_location] where the output .ovm mesh will go, along with a .json file containing some expansion stats"<<std::endl;
    std::cout<<"   ---------------------------------------------------- "<<std::endl;
    std::cout<<"   with [function_index] in:"<<std::endl;
    std::cout<<"        0 - shrink-and-expand method "<<std::endl;
    std::cout<<"        1 - map mesh's boundary to be star-shaped"<<std::endl;
    std::cout<<"   ---------------------------------------------------- "<<std::endl;
    std::cout<<"   [boundary_mapping]: "<<std::endl;
    std::cout<<"        1 - tetrahedral boundary"<<std::endl;
    std::cout<<"        2 - stiff tetrahedral boundary"<<std::endl;
    std::cout<<"        3 - tetrahedral boundary then spherical projection"<<std::endl;
    std::cout<<"        4 - random star-shaped boundary"<<std::endl;
    std::cout<<"   ---------------------------------------------------- "<<std::endl;
    std::cout<<"   [option] "<<std::endl;
    std::cout<<"   for function_index = 0:"<<std::endl;
    std::cout<<"        0 - run expansion in silent mode (default)"<<std::endl;
    std::cout<<"        1 - print out expansion details. "<<std::endl;
    std::cout<<"            Use this if you want to see what's happening with the algorithm."<<std::endl;
    std::cout<<"            Note that this doesn't represent TOO MUCH output and doesn't really impeed on performance."<<std::endl;
    std::cout<<"            (This is NOT a debugging mode)"<<std::endl;
    std::cout<<"   for function_index = 1:"<<std::endl;
    std::cout<<"        0 -    output only export "<<std::endl;
    std::cout<<"        1 - input + output export (including boundary cleanup)"<<std::endl;
    std::cout<<"        2 - input + output export + TLC formatted"<<std::endl;
    std::cout<<"   ---------------------------------------------------- "<<std::endl;
}

/**
 * NOTE ABOUT VALIDATION FILES: To make sure the process actually terminated,
 * the trick is to create a file at the beginning and then to delete it.
 * That way if the file remains at the end it means the process did not terminated.

**/
int main(int argc, char** argv) {


#warning TODO: rename executable for something like "BatchProcessor"
    if(argc < 5 || argc > 6){
        print_usage();
        return -1;
    }

    std::istringstream ss;


    ss = std::istringstream(argv[1]);
    std::string input_mesh_path;
    if (!(ss >> input_mesh_path)) {
        std::cerr << "Invalid input path argument: " << argv[1] << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after location argument: " << argv[1] << '\n';
    }


    ss = std::istringstream(argv[2]);
    std::string output_path;
    if (!(ss >> output_path)) {
        std::cerr << "Invalid output path argument: " << argv[2] << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after location argument: " << argv[2] << '\n';
    }


    int function_index(DEFAULT_FUNCTION_INDEX);

    ss = std::istringstream(argv[3]);
    if (!(ss >> function_index)) {
        std::cerr << "Invalid function index: " << argv[3] << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv[3] << '\n';
    }


    int boundary_map(0);
    if(argc >= 5){
        ss = std::istringstream(argv[4]);
        ss >> boundary_map;
        std::cout<<" map: "<<boundary_map<<std::endl;
    }

    int option(0);
    if(argc >= 6){
        ss = std::istringstream(argv[5]);
        ss >> option;
        std::cout<<" OPTION: "<<option<<std::endl;
    }


    TetrahedralMesh input_mesh;
    TetrahedralMesh result_mesh;

    switch(function_index){

    case 0:{

        std::cout<<"  SHRINKING AND EXPANDING MESH "<<input_mesh_path<<std::endl;

        std::string filename, directory;
        if(!extract_mesh_name_and_directory(input_mesh_path, filename, directory)){
                std::cerr<<" --> couldn't extract directory and filename from location "<<input_mesh_path<<std::endl;
                exit(EXIT_FAILURE);
        }
        std::string boundary = "";
        switch(boundary_map){
        case 1:{
            boundary = "tet_mapped";
            break;
        }
        case 2:{
            boundary = "stiff_tet_mapped";
            break;
        }
        case 3:{
            boundary = "ball_mapped";
            break;
        }
        case 4:{
            boundary = "random_star_shape_mapped";
            break;
        }
        default:{
            std::cout<<" ERROR - unhandled boundary "<<boundary_map<<std::endl;
            print_usage();
            return -1;
        }
        }

        int result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::shrinkAndExpand,
                                                        std::string(input_mesh_path.c_str()),
                                                        true,
                                                        boundary_map,
                                                        option,
                                                        input_mesh,
                                                        result_mesh,
                                                        filename,
                                                        output_path+"/"+filename+"_"+boundary+".json");

        switch(result){
        case 0: {
            std::cout<<" ============================================================================================"<<std::endl;
            std::cout<<" ==== SUCCESFULLY MAPPED MESH TO A STAR-SHAPED DOMAIN WITH STRICTLY POSITIVE-VOLUME TETS ===="<<std::endl;
            std::cout<<" ============================================================================================"<<std::endl;

            writeResultMeshInSubdirectory(input_mesh_path,
                                          output_path,
                                          result_mesh,
                                          "output_"+boundary);
            break;
        }

        case 4: {
            std::cout<<" ====================================="<<std::endl;
            std::cout<<" ==== EXPANSION PROCESS TIMED OUT ===="<<std::endl;
            std::cout<<" ====================================="<<std::endl;
            break;
        }

        case -1: {
            std::cout<<" ===================================================="<<std::endl;
            std::cout<<" ==== FAILED TO MAP MESH TO A STAR-SHAPED DOMAIN ===="<<std::endl;
            std::cout<<" ===================================================="<<std::endl;
            break;
        }
        case -2: {
            std::cout<<" ============================================"<<std::endl;
            std::cout<<" ==== AN ERROR OCCURRED DURING EXPANSION ===="<<std::endl;
            std::cout<<" ============================================"<<std::endl;          
            break;
        }
        default: {
            std::cout<<" ================================================= "<<std::endl;
            std::cout<<" ==== PROCESS TERMINATED WITH UNHANDLED RESULT CODE "<<result<<" ===="<<std::endl;
            std::cout<<" ================================================ "<<std::endl;
            break;
        }
        }

        break;
    }

    case 1:{

        std::cout<<"  MAPPING TO STAR-SHAPED BOUNDARY MESH "<<input_mesh_path<<std::endl;

        int result(1);

        std::string output_sub_directory = "";
        switch(boundary_map){
        case 1:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_unit_tet,
                                                        std::string(input_mesh_path.c_str()),
                                                        true,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "tet_mapped";
            break;
        }
        case 2:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_stiff_unit_tet,
                                                        std::string(input_mesh_path.c_str()),
                                                        true,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "stiff_tet_mapped";
            break;
        }
        case 3:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_unit_ball_using_tet,
                                                        std::string(input_mesh_path.c_str()),
                                                        true,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "ball_mapped";
            break;
        }
        case 4:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_random_star_shape_using_tet,
                                                        std::string(input_mesh_path.c_str()),
                                                        true,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "random_star_shape_mapped";
            break;
        }
        default:{
            std::cout<<" unhandled boundary map "<<boundary_map<<std::endl;
            print_usage();
        }
        }

        result_mesh.collect_garbage();

        switch(result){
        case 0: {
            std::cout<<" ========================================================================"<<std::endl;
            std::cout<<" ======= SUCCESFULLY MAPPED MESH BOUNDARY TO A STAR-SHAPED DOMAIN ======="<<input_mesh_path<<std::endl;
            std::cout<<" ========================================================================"<<std::endl;

            std::cout<<" option = "<<option<<std::endl;
            if(option == 2){


                ProgressiveEmbedder::interior_uniform_smoothing(result_mesh);

                export_input_and_output_for_FOF_and_TLC_comparison(input_mesh_path,
                                                                  output_sub_directory,
                                                                  input_mesh,
                                                                  result_mesh);

            }else{
                if(option == 1){
                    std::cout<<" copying input "<<std::endl;

                    writeResultMeshInSubdirectory(input_mesh_path,
                                                  output_sub_directory,
                                                  input_mesh,
                                                  output_sub_directory+"_input");
                }

                writeResultMeshInSubdirectory(input_mesh_path,
                                              output_sub_directory,
                                              result_mesh,
                                              output_sub_directory);
            }
            break;
        }
        case 1: {
            std::cout<<" =============================================="<<std::endl;
            std::cout<<" ==== FAILED TO MAP MESH TO BE STAR-SHAPED ===="<<input_mesh_path<<std::endl;
            std::cout<<" =============================================="<<std::endl;
            break;
        }

        default: {

            std::cout<<" ====================================================== "<<std::endl;
            std::cout<<" ==== PROCESS TERMINATED WITH UNHANDLED RESULT CODE "<<result<<" ===="<<input_mesh_path<<std::endl;
            std::cout<<" ====================================================== "<<std::endl;
        }
        }

        break;

    }


    default:{
        std::cout<<" ERROR - UNKNOWN FUNCTION INDEX"<<std::endl;
        print_usage();
        return -1;
    }


    }

    return 0;
}




bool extract_mesh_name_and_directory(std::string location,
                                     std::string& filename,
                                     std::string& directory){

    char sep1 = '/';

    size_t directory_end_index = location.rfind(sep1, location.length());

    if (directory_end_index != std::string::npos) {
        //-4 to get rid of .ovm
        filename = location.substr(directory_end_index+1, location.length()-directory_end_index-4);
        directory = location.substr(0, directory_end_index);
    }else{
        filename = location.substr(0, location.length()-4);
        directory = "";
    }

    return true;
}




void writeResultMeshInSubdirectory(std::string location,
                                   std::string sub_directory,
                                   TetrahedralMesh& result_mesh,
                                   std::string suffix){

    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<" --> WRITING RESULT MESH "<<location<<std::endl;
    std::cout<<" result mesh stats: "<<std::endl;
    std::cout<<"  - vertices........"<<result_mesh.n_vertices()<<std::endl;
    std::cout<<"  -    edges........"<<result_mesh.n_edges()<<std::endl;
    std::cout<<"  -    faces........"<<result_mesh.n_faces()<<std::endl;
    std::cout<<"  -    cells........"<<result_mesh.n_cells()<<std::endl;

    OpenVolumeMesh::IO::FileManager fileManager;

    std::string filename, directory;

    if (extract_mesh_name_and_directory(location, filename, directory)) {

        std::string output_location = directory + "/" + sub_directory + "/" + filename + "_"+suffix+".ovm";

        std::cout<<" location = "<<location<<std::endl;
        std::cout<<" filename = "<<filename<<std::endl;
        std::cout<<" directory = "<<directory<<std::endl;
        std::cout<<" sub directory = "<<sub_directory<<std::endl;
        std::cout<<" output file: "<<output_location<<std::endl;

        result_mesh.collect_garbage();
        if(!fileManager.writeFile(output_location,
                                  result_mesh)){
            std::cerr<<" ERROR - COULD NOT WRITE FILE "<<directory + "/"+ sub_directory + "/" + filename + ".ovm"<<std::endl;
        }else{
            std::cout<<" wrote mesh as "<<output_location<<std::endl;
        }
    }else{
        std::cerr<<" ERROR - parsing directory path"<<std::endl;


        //std::cerr<<" ERROR - parsing directory path, performing emergency save at "<<(location + "_" + sub_directory + ".ovm")<<std::endl;
        //fileManager.readFile(location, result_mesh);
        //fileManager.writeFile(location + "_" + sub_directory + ".ovm", mesh);
    }
    std::cout<<" --> DONE"<<std::endl;
    std::cout<<"-------------------------------"<<std::endl;

}


void copyInputMeshInSubdirectory(std::string location, std::string sub_directory){
    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<" --> COPYING MESH "<<location<<std::endl;
    OpenVolumeMesh::IO::FileManager fileManager;

    TetrahedralMesh mesh;

    std::string filename, directory;

    if (extract_mesh_name_and_directory(location, filename, directory)) {


        fileManager.readFile(location, mesh);
        if(!fileManager.writeFile(directory + "/" + sub_directory + "/" + filename + ".ovm",
                                  mesh)){
            std::cerr<<" ERROR - COULD NOT WRITE FILE "<<directory + "/"+ sub_directory + "/" + filename + ".ovm"<<std::endl;
        }
    }else{
        std::cerr<<" ERROR - parsing directory path, performing emergency save at "<<(location + "_" + sub_directory + ".ovm")<<std::endl;
        fileManager.readFile(location, mesh);
        fileManager.writeFile(location + "_" + sub_directory + ".ovm", mesh);
    }
    std::cout<<" --> DONE"<<std::endl;
    std::cout<<"-------------------------------"<<std::endl;
}




void export_input_and_output_for_FOF_and_TLC_comparison(std::string location,
                                                        std::string sub_directory,
                                                        TetrahedralMesh& input,
                                                        TetrahedralMesh& output){

    std::cout<<" ---------------------------------------------"<<std::endl;
    std::cout<<" EXPORTING INPUT (REST) AND OUTPUT (INITIAL) MESHES FOR FOF AND TLC COMPARISONS"<<std::endl;

    input.collect_garbage();
    output.collect_garbage();

    if(input.n_vertices() != output.n_vertices()){
        std::cout<<" ERROR - input mesh has "<<input.n_vertices()<<" but output mesh has "<<output.n_vertices()<<std::endl;
        return;
    }


    std::string filename, directory;

    if (extract_mesh_name_and_directory(location, filename, directory)) {

        std::string tlc_output_location = directory + "/" + sub_directory + "/" + filename + ".ovm/TLC_input";

        //std::cout<<" location = "<<location<<std::endl;
        //std::cout<<" filename = "<<filename<<std::endl;
        //std::cout<<" directory = "<<directory<<std::endl;

        std::cout<<" exporting TLC-formatted mesh at "<<tlc_output_location<<std::endl;
        std::ofstream outfile;
        outfile.open(tlc_output_location);
        if(!outfile.is_open()){
            std::cout<<" ERROR - couldn't output TLC file "<<tlc_output_location<<std::endl;
        }

        int n_input_boundary_vertices(0);
        outfile<<input.n_vertices()<<" "<<3<<std::endl;
        for(auto v: input.vertices()){
            outfile<<input.vertex(v)<<std::endl;
            n_input_boundary_vertices += input.is_boundary(v);
        }

        int n_output_boundary_vertices(0);
        outfile<<output.n_vertices()<<" "<<3<<std::endl;
        for(auto v: output.vertices()){
            outfile<<output.vertex(v)<<std::endl;
            n_output_boundary_vertices += output.is_boundary(v);
        }

        if(n_input_boundary_vertices != n_output_boundary_vertices){
            std::cout<<" ERROR - "<<n_input_boundary_vertices<<" in input but "<<n_output_boundary_vertices<<" boundary vertices in output"<<std::endl;
            return;
        }

        outfile<<output.n_cells()<<" "<<4<<std::endl;
        for(auto c: output.cells()){
            outfile<<output.get_cell_vertices(c)<<std::endl;
        }

        outfile<<n_input_boundary_vertices<<std::endl;
        for(auto v: input.vertices()){
            if(input.is_boundary(v)){
                outfile<<v<<std::endl;
            }
        }
        outfile.close();
        //input mesh copy
        std::string ovm_location = directory + "/" + sub_directory + "/" + filename + ".ovm/";

        std::string ovm_input_mesh_location = ovm_location + "input.ovm";
        std::cout<<" copying .ovm input mesh at "<<ovm_input_mesh_location<<std::endl;
        OpenVolumeMesh::IO::FileManager fileManager;
        if(!fileManager.writeFile(ovm_input_mesh_location, input)){
            std::cerr<<" ERROR - COULD NOT WRITE FILE "<<ovm_input_mesh_location<<std::endl;
        }

        //also copy the result mesh as .ovm for comparison
        std::string ovm_boundary_mapped_mesh_location = ovm_location + sub_directory + ".ovm";
        std::cout<<" copying .ovm boundary-mapped mesh at "<<ovm_boundary_mapped_mesh_location<<std::endl;
        if(!fileManager.writeFile(ovm_boundary_mapped_mesh_location, output)){
            std::cerr<<" ERROR - COULD NOT WRITE FILE "<<ovm_boundary_mapped_mesh_location<<std::endl;
        }

    }else{
        std::cerr<<" ERROR - parsing directory path"<<std::endl;
    }
    std::cout<<" --> DONE"<<std::endl;

    std::cout<<" ---------------------------------------------"<<std::endl;

}





