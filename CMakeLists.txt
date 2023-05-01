

project(ShrinkAndExpand)
cmake_minimum_required(VERSION 3.20)

    
#enable those only when debugging
#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
set (CMAKE_BUILD_TYPE Release)

find_package (OpenVolumeMesh REQUIRED _NO_MODULE)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (CGAL REQUIRED NO_MODULE)

include_directories(ProgressiveEmbedding/)
include_directories(.)

set(SOURCE_FILES
    ShrinkAndExpand.cc
    BatchProcessor.cc
    TetMeshBoundarySplitter.cc
    BadTetFinder.cc
    TopoHelper.cc
    Expander.cc
    ExpansionDataLogger.cc
    ExpansionCone.cc
    ConnectedVertexSubsetIterator.cc
    StarShapifyableExpansionCone.cc
    PrecisionDiagnostics.cc
    TetMapper.cc
    SphereMapper.cc
    LightWeightStopWatch.cc
    OtherHelpers.cc

    ProgressiveEmbedding/TopologicalLink.cc
    ProgressiveEmbedding/TopologicalFaceSet.cc
    ProgressiveEmbedding/ProgressiveEmbedder.cc
    ProgressiveEmbedding/ProgEmbeddingHelpers.cc)

add_executable(ShrinkAndExpand ${SOURCE_FILES})



    
#add_subdirectory(OpenVolumeMesh)


#target_compile_features(ShrinkAndExpand PUBLIC cxx_std_17)
    
#message("gurobi libraries: " ${GUROBI_LIBRARIES})

target_link_libraries(ShrinkAndExpand
    OpenVolumeMesh::OpenVolumeMesh
    Eigen3::Eigen
    CGAL::CGAL
    #OpenMeshCore
    #CGAL
    #Qt5::Core
    #${GUROBI_LIBRARIES}
    #CoMISo
    )


