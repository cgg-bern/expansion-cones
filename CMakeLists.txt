

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

FILE(GLOB SOURCE_FILES src/*.cc)

include_directories(src/)


add_executable(ShrinkAndExpand ${SOURCE_FILES})

    
    
target_link_libraries(ShrinkAndExpand
    OpenVolumeMesh::OpenVolumeMesh
    Eigen3::Eigen
    CGAL::CGAL
    )


