#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"


#include "environment.h"
#include "particle.h"
#include <algorithm>
//#include <omp.h>

namespace CGL {

    Environment::Environment(size_t nx_cells, size_t ny_cells,
                             float cell_width, float cell_height) {
        this->nx_cells = nx_cells;
        this->ny_cells = ny_cells;
        this->cell_width = cell_width;
        this->cell_height = cell_height;

        int n_cells = nx_cells*ny_cells;
        this->ux = (float*)malloc(n_cells*(sizeof(float)));
        this->uy = (float*)malloc(n_cells*(sizeof(float)));

        this->ux_p = (float*)malloc(n_cells*(sizeof(float)));
        this->uy_p = (float*)malloc(n_cells*(sizeof(float)));

        this->T = (float*)malloc(n_cells*(sizeof(float)));
        this->T_p = (float*)malloc(n_cells*(sizeof(float)));

        this->vort = (float*)malloc(n_cells*(sizeof(float)));
        this->vort_f_y = (float*)malloc(n_cells*(sizeof(float)));
        this->vort_f_x = (float*)malloc(n_cells*(sizeof(float)));

        // Initialize fluid field (u_field)
        #pragma omp parallel for
        for (int y = 0; y < ny_cells; y++) {
            for (int x = 0; x < nx_cells; x++) {
                ux[ID(x, y)] = 10. * (rand() / (double) RAND_MAX - 0.5);
                uy[ID(x, y)] = 10. * (rand() / (double) RAND_MAX - 0.5);

                ux_p[ID(x, y)] = 0.0;
                uy_p[ID(x, y)] = 0.0;
                T[ID(x, y)] = 0.0;
                T_p[ID(x, y)] = 000.0;
            }
        }
       T[ID(200, 200)] += 500;
       for (int x = 0; x < nx_cells; x += 50) {
           for (int y = 0; y < ny_cells; y += 50) {
               T[ID(x, y)] += 500;
           }
       }

       // initialize particles in system
       for (int x = 0; x < nx_cells; x += 100) {
           for (int y = 0; y < ny_cells; y += 100) {
               particles_list.push_back(Fuel(Vector2D(x, y), 1.0, 1.0, 0.0, 0.0, this));
           }
       }

       // for (int y = 50; y < 150; y++) {
       //     for (int x = 100; x < 200; x++) {
       //         ux[ID(x, y)] = 34;
       //         uy[ID(x, y)] = -20;
       //     }
       // }
       // for (int y = 90; y < 150; y++) {
       //     for (int x = 250; x < 350; x++) {
       //         ux[ID(x, y)] = 0;
       //         uy[ID(x, y)] = 70;
       //     }
       // }
    }


    void Environment::simulate(float delta_t, Vector3D gravity, vector<InputItem> inputs) {
        // memset(ux_p, 0, nx_cells*ny_cells*sizeof(float));
        // memset(uy_p, 0, nx_cells*ny_cells*sizeof(float));
        // memset(T_p, 0, nx_cells*ny_cells*sizeof(float));
        // updates ux_p, uy_p, T_p from inputs and gravity

        get_from_UI(delta_t, gravity, inputs);
        thermal_buoyancy(uy_p, delta_t);
        simulate_vel(delta_t);
        simulate_temp(delta_t);
    }
    void Environment::get_from_UI(float delta_t, Vector3D gravity, vector<InputItem> inputs) {

        for (InputItem i : inputs) {
            int x = (int) i.pos.x / cell_width;
            int y = (int) i.pos.y / cell_height;
            if (x < nx_cells-1 && y < ny_cells-1 && x > 0 && y > 0) {
                switch (i.input_mode) {
                case temperature:
                    T_p[ID(x, y)] += 500;
                    break;
                }
            }
        }
    }


    void Environment::simulate_vel(float delta_t) {
        add_source(ux, ux_p, delta_t);
        add_source(uy, uy_p, delta_t);

        // adds vorticity
        calc_vorticity();

        add_source(ux, vort_f_x, delta_t);
        add_source(uy, vort_f_y, delta_t);

        project();
        SWAP(ux, ux_p);
        SWAP(uy, uy_p);
        advect(1, ux, ux_p, ux_p, uy_p, delta_t, false);
        advect(1, uy, uy_p, ux_p, uy_p, delta_t, false);
        project();

    }
    void Environment::thermal_buoyancy(float* f, float delta_t) {
        #pragma omp parallel for
        for (int j = 2; j < ny_cells - 1; j++) {
            for (int i = 1; i < nx_cells - 1; i++) {
                float T_here = T[ID(i, j)];
                float T_above = T[ID(i, j - 1)];
                if (T_above < T_here) {
                    f[ID(i, j)] -= 2000.*(T_here - T_above);
                }
            }
        }
    }


    void Environment::calc_vorticity() {

        // create vorticity field
        for (int i = 1; i <= (nx_cells - 2); i++) {
            for (int j = 1; j <= (ny_cells - 2); j++) {
               // calculate curl at point to obtain vorticity at each point
               vort[ID(i, j)] = uy[ID(i + 1, j)] - uy[ID(i - 1, j)] -  (ux[ID(i, j - 1)] - ux[ID(i, j + 1)]);
            }
        }

        //calculate vorticity force
        float epsilon = 2.0;
        for (int i = 1; i <= (nx_cells - 2); i++) {
            for (int j = 1; j <= (ny_cells - 2); j++) {

                // calculate gradient of vorticity field
                Vector2D eta = Vector2D(abs(vort[ID(i, j - 1)]) - abs(vort[ID(i, j + 1)]), abs(vort[ID(i + 1, j)]) - abs(vort[ID(i - 1, j)]));
                Vector2D N = eta / eta.norm();

                // calculate force, using Fedkiw 2001, not sure how to find h...
                Vector3D vort_f_vector = (epsilon * cross(Vector3D(N.x, N.y, 0), Vector3D(0, 0, vort[ID(i, j)])));
                vort_f_x[ID(i, j)] = vort_f_vector.x;
                vort_f_y[ID(i, j)] = vort_f_vector.y;
            }
        }

    }

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    // 4. project
    void Environment::simulate_temp(float delta_t) {
        add_source(T, T_p, delta_t);
        SWAP(T_p, T);
        diffuse(0, T, T_p, T_diff, delta_t, true);
        SWAP(T_p, T);
        advect(0, T, T_p, ux, uy, delta_t, true);
        temp_decay(T, delta_t);
    }

    void Environment::temp_decay(float * T, float delta_t) {
        // need to play with constant
        float k = 500;
        #pragma omp parallel for
        for (int j = 1; j <= (ny_cells - 2);j++) {
            for (int i = 1; i <= (nx_cells - 2); i++) {

                // calculates surrounding heat by averaging neighboring temperatures
                // can also try randomly sampling a direction
                float surround_temp =
                        (T[ID(i, j + 1)] + T[ID(i, j - 1)] + T[ID(i + 1, j)] +
                         T[ID(i - 1, j)]) / 4.0f;

                // Newton's law of cooling
                T[ID(i, j)] = min(T[ID(i, j)],
                                  (1 - delta_t) * surround_temp + (T[ID(i, j)] - surround_temp) * exp(-k * delta_t));
            }
        }
    }

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    void Environment::simulate_particle() {

    }

    void Environment::simulate_smoke() {

    }

    void Environment::add_source(float * curr, float * prev, float delta_t) {
        #pragma omp parallel for
        for (int j = 2; j < ny_cells-2; j++ ) {
            for (int i = 2; i < nx_cells-2; i++) {
                curr[ID(i, j)] += delta_t*prev[ID(i, j)];
            }
        }
            
    }

    void Environment::diffuse(int b, float * x, float * x0, float diff, float dt, bool isTemp) {
        int i, j, k;
        float a = dt * diff * (nx_cells - 2) * (ny_cells - 2);
        for (k = 0; k < 10; k++ ) {
            #pragma omp parallel for
            for (j = 1; j <= (ny_cells - 2); j++ ) {
                for (i = 1; i <= (nx_cells - 2); i++ ) {

                    x[ID(i, j)] = (x0[ID(i, j)] + a * (x[ID(i - 1, j)] + x[ID(i + 1, j)] +
                                                   x[ID(i, j - 1)] + x[ID(i, j + 1)])) / (1 + 4 * a);
                }
            }

            // does not apply boundary conditions for temperature
            if (!isTemp) {
                set_bnd(b, x);
            }
        }
    }

    void Environment::project() {
        int i, j, k;
        float h;
        float* p = ux_p;
        float* div = uy_p;
        h = cell_width;
        #pragma omp parallel for
        for (j = 1; j <= (ny_cells - 2); j++) {
            for (i = 1; i <= (nx_cells - 2); i++) {

                div[ID(i, j)] = -0.5f * h * (ux[ID(i + 1, j)] - ux[ID(i - 1, j)] +
                                uy[ID(i, j + 1)] - uy[ID(i, j - 1)]);
                p[ID(i,j)] = 0;
            }
        }
        set_bnd(0, div ); set_bnd(0, p );
        for (k = 0; k < 20; k++) {
#pragma omp parallel for schedule(dynamic, 4)
            for (j = 1 ; j <= (ny_cells - 2); j++ ) {
                for (i = 1 ; i <= (nx_cells - 2); i++ ) {

                    p[ID(i, j)] = (div[ID(i, j)] + p[ID(i - 1, j)] + p[ID(i + 1, j)] +
                                  p[ID(i, j - 1)] + p[ID(i, j + 1)]) / 4;
                }
            }
            set_bnd(0, p);
        }
        #pragma omp parallel for
        for (j = 1; j<=(ny_cells - 2); j++ ) {
            for (i = 1; i <= (nx_cells - 2); i++ ) {

                ux[ID(i, j)] -= 0.5 * (p[ID(i + 1, j)] - p[ID(i - 1, j)]) / h;
                uy[ID(i, j)] -= 0.5 * (p[ID(i, j + 1)] - p[ID(i, j - 1)]) / h;
            }
        }
        set_bnd(1, ux); set_bnd(2, uy);
    }

    void Environment::advect(int b, float * d, float * d0, float* u, float* v, float dt, bool isTemp) {
        int i, j, i0, j0, i1, j1;
        float x, y, s0, t0, s1, t1, dt0_x, dt0_y;
        dt0_x = dt * (nx_cells - 2);
        dt0_y = dt * (ny_cells - 2);
#pragma omp parallel for schedule(dynamic, 3)
        for (j = 1; j <= (ny_cells - 2); j++ ) {
            for (i = 1; i <= (nx_cells - 2); i++ ) {

                x = i - dt0_x * u[ID(i, j)];
                y = j - dt0_y * v[ID(i, j)];
                if (x < 0.5) {
                    x = 0.5;
                }
                if (x > (nx_cells - 2) + 0.5) {
                    x = (nx_cells - 2) + 0.5f;
                }
                i0 = (int) x;
                i1 = i0 + 1;
                if (y < 0.5) {
                    y = 0.5;
                }
                if (y > (ny_cells - 2) + 0.5) {
                    y = (ny_cells - 2) + 0.5f;
                }
                j0 = (int) y;
                j1 = j0 + 1;
                s1 = x - i0;
                s0 = 1 - s1;
                t1 = y - j0;
                t0 = 1 - t1;
                d[ID(i, j)] = s0 * (t0 * d0[ID(i0, j0)] + t1 * d0[ID(i0, j1)]) +
                    s1 * (t0 * d0[ID(i1, j0)] + t1 * d0[ID(i1, j1)]);
            }
        }

        // does not apply boundary conditions for temperature
        if (!isTemp) {
            set_bnd(b, d);
        }
    }

    void Environment::set_bnd(int b, float * x) {
        int i;
        int j;

        // handles y-coord bound setting
        for (i = 1; i <= (ny_cells - 2); i++) {
            if (b == 1) {
                x[ID(0, i)] = x[ID(1, i)] * -1;
                x[ID((nx_cells - 2) + 1, i)] = x[ID((nx_cells - 2), i)] * -1;
            } else {
                x[ID(0, i)] = x[ID(1, i)];
                x[ID((nx_cells - 2) + 1, i)] = x[ID((nx_cells - 2), i)];
            }

            // converted from following
            //x[ID(0, i)] = b == 1 ? –x[ID(1, i)] : x[ID(1, i)];
            //x[ID(nx_cells + 1, i)] = b == 1 ? –x[ID(nx_cells, i)] : x[ID(nx_cells, i)];
        }

        // handles x-coord bound setting
        for (j = 1; j <= (nx_cells - 2); j++) {
           if (b == 2) {
                x[ID(j, 0)] = x[ID(j, 1)] * -1;
                x[ID(j, (ny_cells - 2) + 1)] = x[ID(j, (ny_cells - 2))] * -1;
            } else {
                x[ID(j, 0)] = x[ID(j, 1)];
                x[ID(j, (ny_cells - 2) + 1)] = x[ID(j, (ny_cells - 2))];
            }

            // converted from following
            //x[ID(j, 0)] = b == 2 ? –x[ID(j, 1)] : x[ID(j, 1)];
            //x[ID(j, ny_cells + 1)] = b == 2 ? –x[ID(j, ny_cells)] : x[ID(j, ny_cells)];
        }


        x[ID(0, 0)] = 0.5f * (x[ID(1, 0)] + x[ID(0, 1)]);
        x[ID(0, (ny_cells - 2) + 1)] = 0.5f * (x[ID(1, (ny_cells - 2) + 1)] + x[ID(0, (ny_cells - 2))]);

        x[ID((nx_cells - 2) + 1, 0)] = 0.5f * (x[ID((nx_cells - 2), 0)] + x[ID((nx_cells - 2) + 1, 1)]);
        x[ID((nx_cells - 2) + 1, (ny_cells - 2) + 1)] = 0.5f * (x[ID((nx_cells - 2), (ny_cells - 2) + 1)] + x[ID((nx_cells - 2) + 1, (ny_cells - 2))]);

    }
}
