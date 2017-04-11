#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"


#include "environment.h"
#include <algorithm>

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

        // Initialize fluid field (u_field)
        for (int y = 0; y < ny_cells; y++) {
            for (int x = 0; x < nx_cells; x++) {
                ux[ID(x, y)] = 0.0;
                uy[ID(x, y)] = 0.0;
                ux_p[ID(x, y)] = 0.0;
                uy_p[ID(x, y)] = 0.0;
                T[ID(x, y)] = 0.0;
                T_p[ID(x, y)] = 0.0;
            }
        }
    }


    void Environment::simulate(float delta_t, Vector3D gravity) {
        // Update environment here

        // Navier Stokes for u_field
        // addForce step
        for (int y = 0; y < ny_cells; y++) {
            for (int x = 0; x < nx_cells; x++) {
                ux[ID(x, y)] += 0.0001;
                uy[ID(x, y)] += 0.0001;
                ux_p[ID(x, y)] += 0.0001;
                uy_p[ID(x, y)] += 0.0001;
                T[ID(x, y)] += 0.0001;
                T_p[ID(x, y)] += 0.0001;
            }
        }
        // Transport step
    }

    void Environment::sim_vel() {

    }

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    // 4. project
    void Environment::sim_temp() {

    }

    // Steps:
    // 1. add_source(force)
    // 2. diffuse
    // 3. advect
    void Environment::sim_particle() {

    }

    void Environment::sim_smoke() {

    }

    void Environment::add_source(float * x, float * s, float dt) {
        int i, size = ((nx_cells - 2) + 2) * ((ny_cells - 2) + 2);
        for (i = 0; i < size; i++ ) x[i] += dt * s[i];
    }

    void Environment::diffuse(int b, float * x, float * x0, float diff, float dt) {
        int i, j, k;
        float a = dt * diff * (nx_cells - 2) * (ny_cells - 2);
        for (k = 0; k < 20; k++ ) {
            for (i = 1; i <= (nx_cells - 2); i++ ) {
                for (j = 1; j <= (ny_cells - 2); j++ ) {
                    x[ID(i, j)] = (x0[ID(i, j)] + a * (x[ID(i - 1, j)] + x[ID(i + 1, j)] +
                                                   x[ID(i, j - 1)] + x[ID(i, j + 1)])) / (1 + 4 * a);
                }
            }
            set_bnd(b, x );
        }
    }

    void Environment::project( float * p, float * div ) {
        int i, j, k;
        float h;
        h = 1.0f / (nx_cells - 2);
        for (i = 1; i <= (nx_cells - 2); i++) {
            for (j = 1; j <= (ny_cells - 2); j++) {
                div[ID(i, j)] = -0.5f * h * (ux_p[ID(i + 1, j)] - ux_p[ID(i - 1, j)] +
                                uy_p[ID(i, j + 1)] - uy_p[ID(i, j - 1)]);
                p[ID(i,j)] = 0;
            }
        }
        set_bnd(0, div ); set_bnd(0, p );
        for (k = 0; k < 20; k++) {
            for (i = 1 ; i <= (nx_cells - 2); i++ ) {
                for (j = 1 ; j <= (ny_cells - 2); j++ ) {
                    p[ID(i, j)] = (div[ID(i, j)] + p[ID(i - 1, j)] + p[ID(i + 1, j)] +
                                  p[ID(i, j - 1)] + p[ID(i, j + 1)]) / 4;
                }
            }
            set_bnd(0, p);
        }
        for (i = 1; i <= (nx_cells - 2); i++ ) {
            for (j = 1; j<=(ny_cells - 2); j++ ) {
                ux_p[ID(i, j)] -= 0.5 * (p[ID(i + 1, j)] - p[ID(i - 1, j)]) / h;
                uy_p[ID(i, j)] -= 0.5 * (p[ID(i, j + 1)] - p[ID(i, j - 1)]) / h;
            }
        }
        set_bnd(1, ux_p); set_bnd(2, uy_p);
    }

    void Environment::advect(int b, float * d, float * d0, float dt ) {
        int i, j, i0, j0, i1, j1;
        float x, y, s0, t0, s1, t1, dt0_x, dt0_y;
        dt0_x = dt * (nx_cells - 2);
        dt0_y = dt * (ny_cells - 2);
        for (i = 1; i <= (nx_cells - 2); i++ ) {
            for (j = 1; j <= (ny_cells - 2); j++ ) {
                x = i - dt0_x * ux_p[ID(i, j)]; y = j - dt0_y * uy_p[ID(i, j)];
                if (x < 0.5) x = 0.5; if (x > (nx_cells - 2) + 0.5) x = (nx_cells - 2) + 0.5f; i0 = (int) x; i1 = i0 + 1;
                if (y < 0.5) y= 0.5; if (y > (ny_cells - 2) + 0.5) y = (ny_cells - 2) + 0.5f; j0 = (int) y; j1 = j0 + 1;
                s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
                d[ID(i, j)] = s0 * (t0 * d0[ID(i0, j0)] + t1 * d0[ID(i0, j1)])+
                                                                            s1 * (t0 * d0[ID(i1, j0)] + t1 * d0[ID(i1, j1)]);
            }
        }
        set_bnd (b, d);
    }

    void Environment::set_bnd(int b, float * x) {
        int i;
        int j;

        // handles y-coord bound setting
        for (i = 1; i <= (nx_cells - 2); i++) {
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
        for (j = 1; j <= (ny_cells - 2); j++) {
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
