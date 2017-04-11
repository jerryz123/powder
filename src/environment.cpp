#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"


#include "environment.h"
#include "cell.h"
#include <algorithm>

#define IX(i,j) ((i)+(N+2)*(j))

namespace CGL {


    Environment::Environment(size_t nx_cells, size_t ny_cells,
                             float cell_width, float cell_height) {

    }


    void Environment::simulate(float delta_t, Vector3D gravity) {
        // Update environment here

        // Navier Stokes for u_field
        // addForce step
//        for (Cell* c : u0_field) {
//            c->u += delta_t*gravity;
//        }
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
