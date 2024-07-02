include(FetchContent)

set(FETCHCONTENT_QUIET OFF)
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)

project(ShrinkAndExpand)
cmake_minimum_required(VERSION 3.20)

    
#enable those only when debugging
#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
set (CMAKE_BUILD_TYPE Release)

#old version
#find_package (OpenVolumeMesh REQUIRED NO_MODULE)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (CGAL REQUIRED NO_MODULE)

if(NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(eigen
            #GIT_REPOSITORY https://gitlab.com/libeigen/eigen
            #GIT_TAG 3.4.0
            # temporary fix, 2025-01-16: Martin Heistermann's MR is not merged yet, but is required to fix a build error:
            GIT_REPOSITORY https://gitlab.com/mheistermann/eigen
            GIT_TAG fix/spqr-index-vs-storageindex
            SOURCE_DIR "${EXTERNAL_DIR}/eigen"
    )
    FetchContent_Populate(eigen)
    message("Downloaded Eigen3 to ${eigen_SOURCE_DIR}")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE "$<BUILD_INTERFACE:${eigen_SOURCE_DIR}>")
    #target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
endif()

if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    FetchContent_Declare(openvolumemesh
            GIT_REPOSITORY  https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git
            GIT_TAG v3.3.0
            SOURCE_DIR "${EXTERNAL_DIR}/OpenVolumeMesh"
    )
    FetchContent_MakeAvailable(openvolumemesh)
endif()


FILE(GLOB SOURCE_FILES src/*.cc)

include_directories(src/)

add_executable(ShrinkAndExpand ${SOURCE_FILES})

target_link_libraries(ShrinkAndExpand
    OpenVolumeMesh::OpenVolumeMesh
    Eigen3::Eigen
    CGAL::CGAL
    )


