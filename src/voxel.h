#ifndef VOXEL_H
#define VOXEL_H

#include "CGL/CGL.h"
#include "CGL/vector2D.h"

using namespace CGL;

struct Voxel {
Voxel(Vector3D u, float T) : u(u), T(T) {}
    Vector3D u;
    float T;
};
#endif
