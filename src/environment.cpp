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

    void Environment::simulate_vel(float delta_t) {
        // add_source(delta_t);
        // SWAP(ux, ux_p);
        // diffuse(1, ux, ux_p, delta_t);
        // diffuse(2, uy, uy_p, delta_t);
        // SWAP(ux, ux_p);
        // SWAP(uy, uy_p);
        // project();
        // SWAP(ux, ux_p);
        // SWAP(uy, uy_p);
        // advect(1, ux, ux_p, delta_t);
        // advect(1, uy, uy_p, delta_t);
        // project();
    }
}
