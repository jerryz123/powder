#ifndef CELL_H
#define CELL_H

#include "CGL/CGL.h"
#include "CGL/vector2D.h"

using namespace CGL;

struct Cell {
Cell(Vector3D u, float T, float phi, Vector3D pos)
: u(u), T(T), phi(phi), pos(pos){}
    Vector3D u;
    float T;
    float phi;
    Vector3D pos;
};
#endif
