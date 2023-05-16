#pragma once

#include "ProgressiveEmbedder.hh"
#include "TetMeshBoundarySplitter.hh"

#include <OpenVolumeMesh/FileManager/FileManager.hh>


class BatchProcessor
{
public:


    /** return value should follow convention
     * -3 un-collapse failure
     * -2 collapse failure
     * -1 pre-processing/input error
     *  0 complete success
     *  >0 partial success, function-dependent */

    using MeshProcessingFunction =  int (*)(TetrahedralMesh&);

    using MeshProcessingFunctionWithFileOutput =  int (*)(TetrahedralMesh&, const std::string& mesh_name, const std::string& output_file_path, int option1, int option2);

    static int collapseSingleMesh(MeshProcessingFunction function,
                                  std::string filename,
                                  bool pre_process,
                                  TetrahedralMesh& input_mesh,
                                  TetrahedralMesh& result_mesh);

    static int collapseSingleMesh(MeshProcessingFunctionWithFileOutput function,
                                  std::string filename,
                                  bool pre_process,
                                  int option1,
                                  int option2,
                                  TetrahedralMesh& input_mesh,
                                  TetrahedralMesh& result_mesh,
                                  const std::string& mesh_name,
                                  const std::string& output_file_path);

};

