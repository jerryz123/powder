#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"


#include "particle.h"
#include "environment.h"
#include <algorithm>

#define ID(x, y) (x + (y)*env->nx_cells)
namespace CGL {
    Particle::Particle(Vector2D position, float radius,
                       float ux, float uy, Environment *env) {
        this->position = position;
        this->radius = radius;

        this->ux = ux;
        this->uy = uy;
        this->env = env;
    }


    void Soot::simulate(float delta_t) {
        int x = (int) position.x;
        int y = (int) position.y;
        
    }

    void Fuel::simulate(float delta_t) {
        int x = (int) position.x;
        int y = (int) position.y;


        if (env->T[ID(x, y)] > ignition_T ||
            is_burning &&
            x > 1 &&
            x < env->nx_cells-1 &&
            y > 1 &&
            y < env->ny_cells-1) {
            if (!is_burning) {
                env->phi[ID(x, y)] += 100;
            }
            is_burning = true;
            radius -= burn_rate * delta_t;
            env->T_p[ID(x, y)] += 100;
            env->smoke_p[ID(x, y)] += 50;


            // particle splitting
            double powder_rand = ((double) rand() / RAND_MAX);
            if (powder_rand < 0.1) {
                double rand_direction = ((double) rand() / RAND_MAX);
                Vector2D dir(cos(2 * PI * rand_direction), 0.5 * sin(2 * PI * rand_direction));
                if (radius > 0) {
                    double rand_force = ((double) rand() / RAND_MAX) * 20;
                    Fuel *new_particle = new Fuel(Vector2D(position.x + dir.x, position.y + dir.y), radius * 0.5, 1, dir.x * rand_force, dir.y * rand_force, env);
                    radius -= radius * 0.5;
                    new_particle->is_burning = true;
                    Particle *temp = new_particle;
                    env->new_particles->push_back(temp);
                }

            }

        }
        if (radius < 0) {
            env->phi[ID(x, y)] -= 200;
        }
    }


}
