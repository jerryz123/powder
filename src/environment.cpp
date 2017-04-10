#include <iostream>
#include <vector>

#include "CGL/vector2D.h"


#include "environment.h"
#include "voxel.h"
#include <algorithm>

namespace CGL {


    void Environment::init() {
        cout << width << height << endl;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Voxel* v = new Voxel(Vector3D(), 0.0);
                u_field.push_back(v);
            }
        }
        cout << u_field.size() << endl;
    }


    void Environment::simulate(float delta_t, Vector2D gravity) {
        // Update environment here
        int a = 1;
    }
}
