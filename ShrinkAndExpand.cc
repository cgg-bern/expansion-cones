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


void copyInputMeshInSubdirectory(std::string location,
                                 std::string sub_directory);


void export_input_and_output_for_FOF_and_TLC_comparison(std::string location,
                                           std::string sub_directory,
                                           TetrahedralMesh& input,
                                           TetrahedralMesh& output);


void writeValidationFile(std::string location);


void removeValidationFile(std::string location);

void dispatchUncollapsingResultCode(int result_code,
                                    std::string location);

void short_test_function(const std::string& mesh_file_path);

void export_mesh_as_OBJ_with_scaled_tets(TetrahedralMesh& mesh,
                                         const std::string OBJ_mesh_file_path,
                                         float tet_scaling);

void smooth_interior(TetrahedralMesh& mesh,
                     int n_iters = 10,
                     double eps = 1e-7);


/**
 * NOTE ABOUT VALIDATION FILES: To make sure the process actually terminated,
 * the trick is to create a file at the beginning and then to delete it.
 * That way if the file remains at the end it means the process did not terminated.

**/
int main(int argc, char** argv) {


#warning TODO: rename executable for something like "BatchProcessor"
    if(argc < 3 || argc > 6){
        std::cerr<<" EXPECTED USAGE: ShrinkAndExpand [location] [function_index] [pre-process] [boundary-mapping] [shrinkage method/export option/domain mesh path]"<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        std::cout<<"   with [function_index] in:"<<std::endl;
        std::cout<<"  -3 - output quality statistics to a .txt file for the given mesh. Assumes it contains two string props containing the domain and codomain positions"<<std::endl;
        std::cout<<"  -2 - prints statistics on the mesh given in argument"<<std::endl;
        std::cout<<"  -1 - whatever's the current short test function (see StandaloneBatchCollaper_main::short_test_function()) "<<std::endl;
        std::cout<<"   7 - map mesh's boundary to be star-shaped"<<std::endl;
        std::cout<<"   8 - shrink-and-expand method. "<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        std::cout<<"   [location] can either be a single .ovm mesh file (if [function_index] is 0) or a directory"<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        std::cout<<"   [pre-process] 1 or 0 (default = 1)"<<std::endl;
        std::cout<<"                 if function_index is -2 or -3, used to pick which method's results are being checed (0:SAE, 1:TLC, 2:FOF)"<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        std::cout<<"   [boundary-mapping]: "<<std::endl;
        std::cout<<"        0 - smoothing until star-shapeness"<<std::endl;
        std::cout<<"        1 - tetrahedral boundary"<<std::endl;
        std::cout<<"        2 - stiff tetrahedral boundary"<<std::endl;
        std::cout<<"        3 - tetrahedral boundary then spherical projection"<<std::endl;
        std::cout<<"        4 - random star-shaped boundary"<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        std::cout<<"   [export option] (for function_index = 7):"<<std::endl;
        std::cout<<"        0 - normal export (output only)"<<std::endl;
        std::cout<<"        1 - input + output export (including boundary cleanup)"<<std::endl;
        std::cout<<"        2 - input + output export + TLC formatted"<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        std::cout<<"   [shrinkage method] (for function_index = 8):"<<std::endl;
        std::cout<<"        0 - full shrinkage"<<std::endl;
        std::cout<<"        1 - edge-based minimal shrinkage"<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        std::cout<<"   [domain mesh path] (for function_index = -3): path to the mesh in the original domain. [pre-process] should be either 1 or 2 in that case "<<std::endl;
        std::cout<<"   ---------------------------------------------------- "<<std::endl;
        return -1;
    }

    std::istringstream ss;


    ss = std::istringstream(argv[1]);
    std::string location;

    if (!(ss >> location)) {
        std::cerr << "Invalid location argument: " << argv[1] << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after location argument: " << argv[1] << '\n';
    }


    int function_index(DEFAULT_FUNCTION_INDEX);
    if(argc >= 3){

        std::cout<<" argv 2: "<<argv[2]<<std::endl;
        ss = std::istringstream(argv[2]);
        if (!(ss >> function_index)) {
            std::cerr << "Invalid function index: " << argv[2] << '\n';
        } else if (!ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[2] << '\n';
        }
    }

    int pre_process(0);
    if(argc >= 4){
        ss = std::istringstream(argv[3]);
        ss >> pre_process;
    }

    std::cout<<"       FUNCTION: "<<function_index<<std::endl;
    std::cout<<" PRE-PROCESSING: "<<pre_process<<std::endl;

    int option1(0);
    if(argc >= 5){
        ss = std::istringstream(argv[4]);
        ss >> option1;
        std::cout<<" OPTION1: "<<option1<<std::endl;
    }

    int option2(0);
    if(argc >= 6){
        ss = std::istringstream(argv[5]);
        ss >> option2;
        std::cout<<" OPTION2: "<<option2<<std::endl;
    }



    writeValidationFile(location);

    TetrahedralMesh input_mesh;
    TetrahedralMesh result_mesh;

    //only used to cases -2 and -3
    std::string method_name = "";
    auto method = pre_process;
    if(method == 0){
        method_name = "SAE";
    }else if(method == 1){
        method_name = "TLC";
    }else if(method == 2){
        method_name = "FOF";
    }

    switch(function_index){

    case -3:{

        std::cout<<" =================================================================================== "<<std::endl;
        std::cout<<" COMPUTING AND OUTPUTTING QUALITY STATISTICS FOR MESH AT "<<location<<std::endl;
        std::cout<<" =================================================================================== "<<std::endl;



        if(method){
            std::cout<<" --> using domain and codomain meshes for computations (method = "<<method_name<<")"<<std::endl;
        }else{
            std::cout<<" --> using single mesh with string properties for computations"<<std::endl;
        }

        TetrahedralMesh codomain_mesh;
        OpenVolumeMesh::IO::FileManager fileManager;
        //std::cout<<" READING MESH "<<filename.toStdString()<<"..."<<std::endl;
        auto file_read_status = fileManager.readFile(location, codomain_mesh);
        if(!file_read_status){
            std::cout<<" couldn't load mesh "<<location<<std::endl;
            return -1;
        }

        std::cout<<" -----------------------------"<<std::endl;
        std::cout<<" - codomain mesh: "<<std::endl;
        std::cout<<" - #vertices: "<<codomain_mesh.n_vertices()<<std::endl;
        std::cout<<" -    #edges: "<<codomain_mesh.n_edges()<<std::endl;
        std::cout<<" -    #faces: "<<codomain_mesh.n_faces()<<std::endl;
        std::cout<<" -    #cells: "<<codomain_mesh.n_cells()<<std::endl;
        std::cout<<" -----------------------------"<<std::endl;
        std::cout<<" ========================================= "<<std::endl;

        std::string filename, super_directory, directory;

        if(!extract_mesh_name_and_directory(location,
                                            filename,
                                            directory)){
                std::cout<<" ERROR - couldn't extract filename and location"<<std::endl;
                return -1;
        }

        auto mesh_name = filename;

        if(!extract_mesh_name_and_directory(directory,
                                            filename,
                                            super_directory)){
            std::cout<<" ERROR - couldn't extract super filename and location"<<std::endl;
            return -1;
        }


        if(method){

            //std::string output_path = super_directory+"/"+method_name+"_quality_stats.txt";
            std::string output_path = directory+"/"+method_name+"_quality_stats.txt";
            std::cout<<" -> writing stats to "<<output_path<<std::endl;

            auto bad_tets = BadTetFinder::findBadTets(codomain_mesh);

            if(!bad_tets.first.empty() || !bad_tets.second.empty()){
                std::cout<<" --> codomain mesh contains degenerate or flipped tets, skipping"<<std::endl;
                return 0;
            }

            TetrahedralMesh domain_mesh;
            ss = std::istringstream(argv[5]);
            std::string domain_mesh_location;
            if(!(ss>>domain_mesh_location)){
                std::cout<<" ERROR - couldn't read domain mesh path "<<argv[5]<<std::endl;
            }

            OpenVolumeMesh::IO::FileManager fileManager;
            //std::cout<<" READING MESH "<<filename.toStdString()<<"..."<<std::endl;
            auto file_read_status = fileManager.readFile(domain_mesh_location, domain_mesh);
            if(!file_read_status){
                std::cout<<" couldn't load mesh "<<domain_mesh_location<<std::endl;
                return -1;
            }

            std::cout<<" -----------------------------"<<std::endl;
            std::cout<<" - domain mesh: "<<std::endl;
            std::cout<<" - #vertices: "<<domain_mesh.n_vertices()<<std::endl;
            std::cout<<" -    #edges: "<<domain_mesh.n_edges()<<std::endl;
            std::cout<<" -    #faces: "<<domain_mesh.n_faces()<<std::endl;
            std::cout<<" -    #cells: "<<domain_mesh.n_cells()<<std::endl;
            std::cout<<" -----------------------------"<<std::endl;

            if(domain_mesh.n_vertices() != codomain_mesh.n_vertices()){
                std::cout<<" ERROR - domain and codomain meshes don't have the same number of vertices"<<std::endl;
                return -1;
            }

            PEHelpers::compute_and_output_quality_statistics(domain_mesh, codomain_mesh, output_path);

        }else{

            //std::string output_path = super_directory+"/"+method_name+"_quality_stats.txt";
            std::string output_path = directory+"/"+mesh_name+"_"+method_name+"_quality_stats.txt";
            std::cout<<" -> writing stats to "<<output_path<<std::endl;

            PEHelpers::compute_and_output_quality_statistics(codomain_mesh, output_path);
        }


        break;
    }

    case -2:{

        std::cout<<" ========================================= "<<std::endl;
        std::cout<<" PRINTING INFORMATION FOR MESH AT "<<location<<std::endl;
        std::cout<<" ========================================= "<<std::endl;
        TetrahedralMesh mesh;
        OpenVolumeMesh::IO::FileManager fileManager;
        //std::cout<<" READING MESH "<<filename.toStdString()<<"..."<<std::endl;
        auto file_read_status = fileManager.readFile(location, mesh);
        if(!file_read_status){
            std::cout<<" couldn't load mesh "<<location<<std::endl;
            return -1;
        }


        auto bad_tets = BadTetFinder::findBadTets(mesh);
        auto entangled_edges_count = PEHelpers::count_entangled_edges(mesh);

        std::cout<<" -----------------------------"<<std::endl;
        std::cout<<" - #vertices: "<<mesh.n_vertices()<<std::endl;
        std::cout<<" -    #edges: "<<mesh.n_edges()<<std::endl;
        std::cout<<" -    #faces: "<<mesh.n_faces()<<std::endl;
        std::cout<<" -    #cells: "<<mesh.n_cells()<<std::endl;
        std::cout<<" -----------------------------"<<std::endl;
        std::cout<<" - degenerate cells: "<<bad_tets.first.size()<<std::endl;
        std::cout<<" -    flipped cells: "<<bad_tets.second.size()<<std::endl;
        std::cout<<" -----------------------------"<<std::endl;
        std::cout<<"    entangled edges: "<<entangled_edges_count<<std::endl;
        std::cout<<" ========================================= "<<std::endl;

        std::string filename, super_directory, directory;

        if(!extract_mesh_name_and_directory(location,
                                            filename,
                                            directory)){
                std::cout<<" ERROR - couldn't extract filename and location"<<std::endl;
                return -1;
        }
        if(!extract_mesh_name_and_directory(directory,
                                            filename,
                                            super_directory)){
            std::cout<<" ERROR - couldn't extract super filename and location"<<std::endl;
            return -1;
        }


        std::string output_path = directory+"/"+method_name+"_stats.txt";
        std::cout<<" -> writing stats to "<<output_path<<std::endl;

        std::ofstream output_file(output_path);
        if(!output_file.is_open()){
            std::cout<<" ERROR - couldn't open stat file "<<output_path<<std::endl;
            return -1;
        }

        std::cout<<" location: "<<location<<std::endl;
        std::cout<<" filename: "<<filename<<std::endl;
        std::cout<<" directory: "<<directory<<std::endl;
        std::cout<<" super_directory: "<<super_directory<<std::endl;

        output_file<<filename<<std::endl;
        output_file<<mesh.n_vertices()<<std::endl;
        output_file<<mesh.n_edges()<<std::endl;
        output_file<<mesh.n_faces()<<std::endl;
        output_file<<mesh.n_cells()<<std::endl;
        output_file<<bad_tets.first.size()<<std::endl;
        output_file<<bad_tets.second.size()<<std::endl;
        output_file<<entangled_edges_count<<std::endl;

        output_file.close();

        break;
    }

    case 7:{


        std::cout<<"  MAPPING TO STAR-SHAPED BOUNDARY MESH "<<location<<std::endl;


        int result(1);

        std::string output_sub_directory = "";
        switch(option1){
        case 1:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_unit_tet,
                                                        std::string(location.c_str()),
                                                        pre_process,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "tet_mapped";
            break;
        }
        case 2:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_stiff_unit_tet,
                                                        std::string(location.c_str()),
                                                        pre_process,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "stiff_tet_mapped";
            break;
        }
        case 3:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_unit_ball_using_tet,
                                                        std::string(location.c_str()),
                                                        pre_process,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "ball_mapped";
            break;
        }
        case 4:{
            result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::map_to_random_star_shape_using_tet,
                                                        std::string(location.c_str()),
                                                        pre_process,
                                                        input_mesh,
                                                        result_mesh);
            output_sub_directory = "random_star_shape_mapped";
            break;
        }
        default:{
            std::cout<<" unhandled option1 "<<option1<<std::endl;
            return -1;
        }
        }


        //OpenVolumeMesh::IO::FileManager fileManager;
        //fileManager.writeFile("aa_output.ovm", result_mesh);

        result_mesh.collect_garbage();

        switch(result){
        case 0: {
            std::cout<<" ========================================================================"<<std::endl;
            std::cout<<" ==== SUCCESFULLY MAPPED MESH TO BE STAR-SHAPED WITH NO FLIPPED TETS ===="<<location<<std::endl;
            std::cout<<" ========================================================================"<<std::endl;
            //copyInputMeshInSubdirectory(location, "ball_mapped");

            std::cout<<" option 2 = "<<option2<<std::endl;
            if(option2 == 2){


                ProgressiveEmbedder::interior_uniform_smoothing(result_mesh);

                export_input_and_output_for_FOF_and_TLC_comparison(location,
                                                                  output_sub_directory,
                                                                  input_mesh,
                                                                  result_mesh);

            }else{
                if(option2 == 1){
                    std::cout<<" copying input "<<std::endl;

                    writeResultMeshInSubdirectory(location,
                                                  output_sub_directory,
                                                  input_mesh,
                                                  output_sub_directory+"_input");
                }

                writeResultMeshInSubdirectory(location,
                                              output_sub_directory,
                                              result_mesh,
                                              output_sub_directory);
            }
            removeValidationFile(location);
            break;
        }
        case 1: {
            std::cout<<" =============================================="<<std::endl;
            std::cout<<" ==== FAILED TO MAP MESH TO BE STAR-SHAPED ===="<<location<<std::endl;
            std::cout<<" =============================================="<<std::endl;
            copyInputMeshInSubdirectory(location, "boundary_mapping_failures");
            break;
        }

        default: {

            std::cout<<" ====================================================== "<<std::endl;
            std::cout<<" ==== PROCESS TERMINATED WITH UNHANDLED RESULT CODE "<<result<<" ===="<<location<<std::endl;
            std::cout<<" ====================================================== "<<std::endl;
            removeValidationFile(location);
        }
        }

        break;

    }

    case 8:{

        std::cout<<"  SHRINKING AND EXPANDING MESH "<<location<<std::endl;

        std::string filename, directory;
        if(!extract_mesh_name_and_directory(location, filename, directory)){
                std::cerr<<" --> couldn't extract directory and filename from location "<<location<<std::endl;
                exit(EXIT_FAILURE);
        }
        std::string boundary = "";
        switch(option1){
        case 1:{
            boundary = "tet_mapped/";
            break;
        }
        case 2:{
            boundary = "stiff_tet_mapped/";
            break;
        }
        case 3:{
            boundary = "ball_mapped/";
            break;
        }
        case 4:{
            boundary = "random_star_shape_mapped/";
            break;
        }
        default:{
            std::cout<<" ERROR - unhandled boundary option "<<option1<<std::endl;
            break;
        }
        }

        int result = BatchProcessor::collapseSingleMesh(ProgressiveEmbedder::shrinkAndExpand,
                                                        std::string(location.c_str()),
                                                        pre_process,
                                                        option1,
                                                        option2,
                                                        input_mesh,
                                                        result_mesh,
                                                        filename,
                                                        directory+"/expansion_data/"+boundary+filename+".json");

        switch(result){
        case 0: {
            std::cout<<" ============================================================================================"<<std::endl;
            std::cout<<" ==== SUCCESFULLY MAPPED MESH TO A STAR-SHAPED DOMAIN WITH STRICTLY POSITIVE-VOLUME TETS ===="<<std::endl;
            std::cout<<" ============================================================================================"<<std::endl;

            writeResultMeshInSubdirectory(location,
                                          "full_injective_map/"+boundary,
                                          result_mesh);
            break;
        }
        case 1: {
            std::cout<<" =============================================================================="<<std::endl;
            std::cout<<" ==== SUCCESFULLY MAPPED MESH TO A STAR-SHAPED DOMAIN BUT SHRINKAGE FAILED ===="<<std::endl;
            std::cout<<" =============================================================================="<<std::endl;
            //copyInputMeshInSubdirectory(location, "shrink_and_expand_failures");
            writeResultMeshInSubdirectory(location,
                                          "shrink_and_expand_failures",
                                          result_mesh);
            break;
        }
        case 2: {
            std::cout<<" ============================================================================================"<<std::endl;
            std::cout<<" ==== SUCCESFULLY MAPPED MESH TO A STAR-SHAPED DOMAIN AND SHRUNK IT BUT EXPANSION FAILED ===="<<std::endl;
            std::cout<<" ============================================================================================"<<std::endl;
            //copyInputMeshInSubdirectory(location, "shrink_and_expand_failures");
            writeResultMeshInSubdirectory(location,
                                          "shrink_and_expand_failures",
                                          result_mesh);
            break;
        }
            //TODO: handle this stuff better
        case 4: {
            std::cout<<" ====================================="<<std::endl;
            std::cout<<" ==== EXPANSION PROCESS TIMED OUT ===="<<std::endl;
            std::cout<<" ====================================="<<std::endl;
            //copyInputMeshInSubdirectory(location, "shrink_and_expand_failures");
            /*writeResultMeshInSubdirectory(location,
                                          "shrink_and_expand_failures",
                                          result_mesh);*/
            break;
        }
        case 5: {
            std::cout<<" ============================================"<<std::endl;
            std::cout<<" ==== EXPANSION REACHED MAX GROWTH RATIO ===="<<std::endl;
            std::cout<<" ============================================"<<std::endl;
            //copyInputMeshInSubdirectory(location, "shrink_and_expand_failures");
            /*writeResultMeshInSubdirectory(location,
                                          "shrink_and_expand_failures",
                                          result_mesh);*/
            break;
        }
        case 6: {
            std::cout<<" ====================================================================================="<<std::endl;
            std::cout<<" ==== SUCCESFULLY MAPPED MESH TO A STAR-SHAPED DOMAIN BUT REQUIRED FULL SHRINKING ===="<<std::endl;
            std::cout<<" ====================================================================================="<<std::endl;
            //copyInputMeshInSubdirectory(location, "shrink_and_expand_failures");
            removeValidationFile(location);
            writeResultMeshInSubdirectory(location,
                                          "full_shrinking_necessary",
                                          result_mesh);
            break;
        }
        case -1: {
            std::cout<<" ===================================================="<<std::endl;
            std::cout<<" ==== FAILED TO MAP MESH TO A STAR-SHAPED DOMAIN ===="<<std::endl;
            std::cout<<" ===================================================="<<std::endl;
            //removeValidationFile(location);
            copyInputMeshInSubdirectory(location, "boundary_mapping_failures");

            break;
        }
        case -2: {
            std::cout<<" ============================================"<<std::endl;
            std::cout<<" ==== AN ERROR OCCURRED DURING EXPANSION ===="<<std::endl;
            std::cout<<" ============================================"<<std::endl;
            //copyInputMeshInSubdirectory(location, "shrink_and_expand_failures");
            writeResultMeshInSubdirectory(location,
                                          "shrink_and_expand_failures",
                                          result_mesh);
            break;
        }
        case -3:{
            std::cout<<" =========================================================================================================================="<<std::endl;
            std::cout<<" ==== SUCCESFULLY MAPPED MESH TO A STAR-SHAPED DOMAIN WITH STRICTLY POSITIVE-VOLUME TETS WITH INITIAL INTERIOR MAPPING ===="<<std::endl;
            std::cout<<" =========================================================================================================================="<<std::endl;
            copyInputMeshInSubdirectory(location, "initial_mapping_ok");
            writeResultMeshInSubdirectory(location,
                                          "full_injective_map",
                                          result_mesh);
            break;
        }
        default: {

            std::cout<<" ================================================= "<<std::endl;
            std::cout<<" ==== PROCESS TERMINATED WITH UNHANDLED RESULT CODE "<<result<<" ===="<<std::endl;
            std::cout<<" ================================================ "<<std::endl;
            //copyInputMeshInSubdirectory(location, "shrink_and_expand_failures");
            writeResultMeshInSubdirectory(location,
                                          "shrink_and_expand_failures",
                                          result_mesh);
            break;
        }
        }

        removeValidationFile(location);



        break;
    }

    default:{
        std::cout<<" ERROR - UNKNOWN FUNCTION INDEX"<<std::endl;
    }


    }

    return 0;
}




void dispatchUncollapsingResultCode(int result_code,
                                    std::string location){

    switch(result_code){
    case-3:{
        std::cout<<" ==================================================== "<<std::endl;
        std::cout<<" ==== AN ERROR OCCURRED WHILE UN-COLLAPSING MESH ===="<<location<<std::endl;
        std::cout<<" ==================================================== "<<std::endl;
        copyInputMeshInSubdirectory(location, "uncollapsing_failures");
        break;

    }

    case -2:{
        std::cout<<" ============================================ "<<std::endl;
        std::cout<<" ==== FAILED TO THREE-STEP COLLAPSE MESH ===="<<location<<std::endl;
        std::cout<<" ============================================ "<<std::endl;
        removeValidationFile(location);
        break;
    }

    case -1:{
        std::cout<<" ================================ "<<std::endl;
        std::cout<<" ==== PRE-PROCESSING FAILURE ===="<<location<<std::endl;
        std::cout<<" ================================ "<<std::endl;
        removeValidationFile(location);
        break;
    }

    case 0:{
        std::cout<<" ========================================"<<std::endl;
        std::cout<<" ==== SUCCESSFULLY UN-COLLAPSED MESH ===="<<location<<std::endl;
        std::cout<<" ========================================"<<std::endl;
        removeValidationFile(location);
        break;
    }

    case 1:{

        std::cout<<" ============================================================== "<<std::endl;
        std::cout<<" ==== UN-COLLAPSED MESH IS NOT EQUIVALENT TO ORIGINAL MESH ===="<<std::endl;
        std::cout<<" ============================================================== "<<std::endl;
        copyInputMeshInSubdirectory(location, "uncollapsing_failures");
        break;
    }
    default:{

        std::cout<<" ============================= "<<std::endl;
        std::cout<<" ==== UNHANDLED RESULT VALUE "<<location<<std::endl;
        std::cout<<" ============================= "<<std::endl;
        copyInputMeshInSubdirectory(location, "uncollapsing_failures");
        break;
    }
    }
}



bool extract_mesh_name_and_directory(std::string location,
                                     std::string& filename,
                                     std::string& directory){

    char sep1 = '/';

    size_t directory_end_index = location.rfind(sep1, location.length());

    if (directory_end_index != std::string::npos) {
        //-4 to get rid of .ovm
        filename = location.substr(directory_end_index+1, location.length()-directory_end_index-5);
        directory = location.substr(0, directory_end_index);

        return true;
    }

    return false;
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


void writeValidationFile(std::string location){
    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<" --> WRITING VALIDATION FILE FOR MESH "<<location<<std::endl;

    std::string filename, directory;
    if (extract_mesh_name_and_directory(location, filename, directory)) {

        std::ofstream validation_file;
        validation_file.open(directory + "/validation/" + filename + ".txt");
        validation_file<<"ok"<<std::endl;

        if(!validation_file.is_open()){
            std::cout<<" ERROR - COULD NOT WRITE FILE "<<directory + "/validation/" + filename<<std::endl;
            return;
        }

        validation_file.close();

    }else{
        std::cout<<" --> ERROR, couldn't extract validation file name from location "<<location<<std::endl;
    }
    std::cout<<" --> DONE"<<std::endl;
    std::cout<<"-------------------------------"<<std::endl;
}


void removeValidationFile(std::string location){
    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<" --> REMOVING VALIDATION FILE FOR MESH "<<location<<std::endl;

    std::string filename, directory;
    if (extract_mesh_name_and_directory(location, filename, directory)) {

        if(!remove(std::string(directory + "/validation/" + filename + ".txt").c_str())){
            std::cout<<" --> SUCCESFULLY REMOVED VALIDATION FILE"<<std::endl;
        }else{

        }
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





