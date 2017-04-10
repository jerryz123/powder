#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "CGL/CGL.h"
#include "voxel.h"
using namespace std;

namespace CGL {

class Environment {
public:
  /* Rope(vector<Mass *> &masses, vector<Spring *> &springs) */
  /*     : masses(masses), springs(springs) {} */
  /* Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, */
  /*      vector<int> pinned_nodes); */

    Environment(size_t width, size_t height)
        : width(width), height(height) {}
    void init();
  void simulate(float delta_t, Vector2D gravity);

  /* Keep track of some variable shere */
  size_t width;
  size_t height;
  vector<Voxel*> u_field;
  /* vector<Mass *> masses; */
  /* vector<Spring *> springs; */
}; 
}
#endif 
