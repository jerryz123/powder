#ifndef PARTICLE_H
#define PARTICLE_H

#include "CGL/CGL.h"
#include "CGL/vector2D.h"
#include "environment.h"

using namespace std;

namespace CGL {
class Particle {
    public:
        Particle(Vector2D position, float radius, float density, float u_x, float u_y, Environment *env);

        Vector2D position;
        const double PI = 3.1415927;
        float burn_rate;
        float density;
        float radius;
        float ux;
        float uy;
        Environment *env;

};

class Soot : public Particle {
    public:
        Soot(Vector2D position, float radius, float density, float ux, float uy, Environment *env) : Particle (position, radius, density, ux, uy, env) {
            burn_rate = 0.1;
        }
        void simulate(float delta_t);

};

class Fuel : public Particle {
    public:
        Fuel(Vector2D position, float radius, float density, float ux, float uy, Environment *env) : Particle (position, radius, density, ux, uy, env) {
            burn_rate = 0.5;
        }
        void simulate(float delta_t);

};
}
#endif //PARTICLE_H
