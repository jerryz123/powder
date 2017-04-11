#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "CGL/CGL.h"
#include "CGL/vector3D.h"

#define ID(x, y) (x + y*nx_cells)
using namespace std;

namespace CGL {

class Environment {
public:
  /* Rope(vector<Mass *> &masses, vector<Spring *> &springs) */
  /*     : masses(masses), springs(springs) {} */
  /* Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, */
  /*      vector<int> pinned_nodes); */

    Environment(size_t nx_cells, size_t ny_cells,
                float cell_width, float cell_height);
    void simulate(float delta_t, Vector3D gravity);

    /* Keep track of some variable shere */
    size_t nx_cells;
    size_t ny_cells;
    float cell_width;
    float cell_height;
    const float rho = 1.0;
    float* ux;
    float* uy;
    float* ux_p;
    float* uy_p;

    float* T;
    float* T_p;
}; 
}
#endif 
