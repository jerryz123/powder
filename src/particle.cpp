#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"


#include "particle.h"
#include "environment.h"
#include <algorithm>
//#include <omp.h>

namespace CGL {
    Particle::Particle(Vector2D position, float radius, float density, float ux, float uy, Environment *env) {
        this->position = position;
        this->radius = radius;
        this->density = density;
        this->ux = ux;
        this->uy = uy;
        this->env = env;
    }

    void Soot::simulate(float delta_t) {
        //mass = area * density
        float mass = PI * radius * radius * density;

        // calculate new radius
        radius = sqrt((mass - burn_rate) / (PI * density));
    }

    void Fuel::simulate(float delta_t) {
        //mass = area * density
        float mass = PI * radius * radius * density;

        // calculate new radius
        radius = sqrt((mass - burn_rate) / (PI * density));
    }
}
