#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "CGL/CGL.h"
#include "CGL/vector3D.h"

#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}
#define ID(x, y) (x + y*nx_cells)

using namespace std;

namespace CGL {

class Environment {
 public:

    Environment(size_t nx_cells, size_t ny_cells,
                float cell_width, float cell_height);
    void simulate(float delta_t, Vector3D gravity);
    
    /* Keep track of some variable shere */
    size_t nx_cells;
    size_t ny_cells;
    float cell_width;
    float cell_height;

    const float rho = 1.0;
    const float visc = 1.0;
    float* ux;
    float* uy;
    float* ux_p;
    float* uy_p;
    float* T;
    float* T_p;
 private:
    void simulate_vel(float delta_t);
}; 
}
#endif 
