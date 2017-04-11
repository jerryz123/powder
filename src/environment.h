#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "CGL/CGL.h"
#include "CGL/vector3D.h"
#include "cell.h"
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

    // Each call to simulate involves
    // 1. sim_particle()
    // 2. sim_vel();
    // 3. sim_t();
    // 4. sim_smoke();

    void simulate(float delta_t, Vector3D gravity);

    void sim_vel();

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    // 4. project
    void sim_temp();

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    void sim_particle();

    void sim_smoke();

    void add_source(int N, int M, float * x, float * s, float dt);
    void diffuse(int N, int M, int b, float * x, float * x0, float diff, float dt);
    void project();
    void advect(int N, int M, int b, float * d, float * d0, float * u, float * v, float dt );
    void set_bnd(int N, int M, int b, float * x );

    /* Keep track of some variable shere */
    size_t nx_cells;
    size_t ny_cells;
    float cell_width;
    float cell_height;
    const float rho = 1.0;

    vector<float> ux_field;
    vector<float> uy_field;
    vector<Vector3D> velo_field;
    vector<float> temp_field;
    vector<float> density_field;

    /* vector<Mass *> masses; */
    /* vector<Spring *> springs; */
}; 
}
#endif 
