#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"


#include "environment.h"
#include "cell.h"
#include <algorithm>

namespace CGL {


    Environment::Environment(size_t nx_cells, size_t ny_cells,
                             float cell_width, float cell_height) {
        this->nx_cells = nx_cells;
        this->ny_cells = ny_cells;
        this->cell_width = cell_width;
        this->cell_height = cell_height;
        
        // Initialize fluid field (u_field)
        for (int y = 0; y < ny_cells; y++) {
            for (int x = 0; x < nx_cells; x++) {
                Cell* v = new Cell(Vector3D(0.0,
                                            0.0,
                                            0.0),
                                   0.0,
                                   0.0,
                                   Vector3D(x*cell_width,
                                            y*cell_height,
                                            0));
                u0_field.push_back(v);
                Cell* vp = new Cell(Vector3D(0.0,
                                              0.0,
                                              0.0),
                                    0.0,
                                    0.0,
                                    Vector3D(x*cell_width,
                                             y*cell_height,
                                             0));
                u1_field.push_back(vp);
                
            }
        }
    }


    void Environment::simulate(float delta_t, Vector3D gravity) {
        // Update environment here

        // Navier Stokes for u_field
        // addForce step
        for (Cell* c : u0_field) {
            c->u += delta_t*gravity;
        }
        // Transport step
    }
}
