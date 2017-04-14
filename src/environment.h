#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "CGL/CGL.h"
#include "CGL/vector3D.h"

#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}
#define ID(x, y) ((x) + (y)*nx_cells)

using namespace std;

namespace CGL {

    enum InputMode {
        temperature
    };
    struct InputItem {
        Vector2D pos;
        InputMode input_mode;
    };

class Environment {
 public:

    Environment(size_t nx_cells, size_t ny_cells,
                float cell_width, float cell_height);

    // Each call to simulate involves
    // 1. sim_particle()
    // 2. sim_vel();
    // 3. sim_t();
    // 4. sim_smoke();

        size_t nx_cells;
    size_t ny_cells;
    float cell_width;
    float cell_height;

    const float rho = 1.0;

    const float visc = 1.0;
    const float T_diff = 0.005;
    const float u_diff = 0;
    float* ux;
    float* uy;
    float* ux_p;
    float* uy_p;
    float* T;
    float* T_p;



    void simulate(float delta_t, Vector3D gravity, vector<InputItem> inputs);

    void simulate_vel(float delta_t);

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    // 4. project
    void simulate_temp(float delta_t);
    void temp_decay(float * T, float delta_t);

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    void simulate_particle();
    void thermal_buoyancy(float* f, float delta_t);
    void simulate_smoke();
    void get_from_UI(float delta_t, Vector3D gravity, vector<InputItem> inputs);
    void add_source(float * curr, float * prev, float delta_t);
    void diffuse(int b, float * x, float * x0, float diff, float dt, bool isTemp);
    void project();
    void advect(int b, float * d, float * d0, float* u, float* v, float dt, bool isTemp);
    void set_bnd(int b, float * x );


}; 
}
#endif 
