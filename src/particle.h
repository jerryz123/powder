#ifndef PARTICLE_H
#define PARTICLE_H

#include "CGL/CGL.h"
#include "CGL/vector2D.h"
#include "environment.h"

using namespace std;

namespace CGL {
class Particle {
    public:
        Particle(Vector2D position, float radius,
                 float u_x, float u_y, Environment *env);

        Vector2D position;
        const double PI = 3.1415927;
        float density;
        float radius;
        float ux;
        float uy;
        Environment *env;
        virtual void simulate(float delta_t) {
            cout << "Should not be called" << endl;
        };

};

class Soot : public Particle {
    public:
        Soot(Vector2D position, float radius,
             float ux, float uy, Environment *env)
            : Particle (position, radius, ux, uy, env) {
        }
        virtual void simulate(float delta_t);

};

class Fuel : public Particle {
    public:
        Fuel(Vector2D position, float radius, 
             float ux, float uy, Environment *env)
            : Particle (position, radius, ux, uy, env) {
            burn_rate = 10;
            density = 0.5;
            is_burning = false;
            ignition_T = 0.5;
        }
    float burn_rate;
    bool is_burning;
    float ignition_T;
    
    virtual void simulate(float delta_t);


};
}
#endif //PARTICLE_H
