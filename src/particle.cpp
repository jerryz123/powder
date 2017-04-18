#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"


#include "particle.h"
#include "environment.h"
#include <algorithm>

#define ID(x, y) (x + (y)*env->nx_cells)
namespace CGL {
    Particle::Particle(Vector2D position, float radius, float density,
                       float ux, float uy, Environment *env) {
        this->position = position;
        this->radius = radius;
        this->density = density;
        this->ux = ux;
        this->uy = uy;
        this->env = env;
    }


    void Soot::simulate(float delta_t) {

    }

    void Fuel::simulate(float delta_t) {
        int x = (int) position.x;
        int y = (int) position.y;

        if (env->T[ID(x, y)] > ignition_T ||
            is_burning) {
            is_burning = true;
            radius -= burn_rate * delta_t;
            env->T_p[ID(x, y)] += 200;
            env->phi[ID(x, y)] += 5;
            env->smoke_p[ID(x, y)] += 10;
        }
    }
}
