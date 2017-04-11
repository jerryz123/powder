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

    void Environment::add_source(int N, int M, float * x, float * s, float dt) {
        int i, size = (N + 2) * (M + 2);
        for ( i = 0 ; i < size ; i++ ) x[i] += dt * s[i];
    }

    void Environment::diffuse(int N, int M, int b, float * x, float * x0, float diff, float dt) {
        int i, j, k;
        float a = dt * diff * N * M;
        for ( k = 0 ; k < 20 ; k++ ) {
            for ( i = 1 ; i <= N ; i++ ) {
                for ( j=1 ; j <= M ; j++ ) {
                    x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                                   x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
                }
            }
            set_bnd ( N, M, b, x );
        }
    }

    void Environment::project() {

    }

    void Environment::advect(int N, int M, int b, float * d, float * d0, float * u, float * v, float dt ) {

    }

    void Environment::set_bnd(int N, int M, int b, float * x) {
        int i;
        for ( i=1 ; i<=N ; i++ ) {
            if (b == 1) {
                x[IX(0 ,i)] = –1 * x[IX(1,i)];
                x[IX(N+1,i)] = b==1 ? –x[IX(N,i)] : x[IX(N,i)];
            } else {

            }
            if (b == 2) {

            } else {
                x[IX(0 ,i)] = x[IX(1,i)];
            }


            x[IX(0 ,i)] = b==1 ? –x[IX(1,i)] : x[IX(1,i)];
            x[IX(N+1,i)] = b==1 ? –x[IX(N,i)] : x[IX(N,i)];
            x[IX(i,0 )] = b==2 ? –x[IX(i,1)] : x[IX(i,1)];
            x[IX(i,N+1)] = b==2 ? –x[IX(i,N)] : x[IX(i,N)];
        }

        x[IX(0 ,0 )] = 0.5 * (x[IX(1,0 )]+x[IX(0 ,1)]);
        x[IX(0 ,N+1)] = 0.5 * (x[IX(1,N+1)]+x[IX(0 ,N )]);
        x[IX(N+1,0 )] = 0.5 * (x[IX(N,0 )]+x[IX(N+1,1)]);
        x[IX(N+1,N+1)] = 0.5 * (x[IX(N,N+1)]+x[IX(N+1,N )]);
    }
}
