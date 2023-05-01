#pragma once

#pragma clang diagnostic ignored "-W#warnings"
#pragma clang diagnostic ignored "-Woverloaded-virtual"

#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>


#include "LightWeightStopWatch.hh"

#include <list>


typedef OpenVolumeMesh::Vec3d Vec3d;
typedef OpenVolumeMesh::TetrahedralGeometryKernel<Vec3d, OpenVolumeMesh::TetrahedralMeshTopologyKernel> TetrahedralMesh;


typedef std::vector<OpenVolumeMesh::VertexHandle> FaceVertices;
typedef std::vector<FaceVertices> FaceVerticesVector;

