#include <iostream>

#include "application.h"
#include "environment.h"

namespace CGL {

Application::Application(AppConfig config) { this->config = config; }

Application::~Application() {}

void Application::init() {
  // Enable anti-aliasing and circular points.
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  glPointSize(1);
  glLineWidth(4);

  glColor3f(1.0, 1.0, 1.0);
  is_simulating = true;
  env = new Environment(config.nx_cells, config.ny_cells,
                        config.cell_width, config.cell_height);

}

void Application::render() {
    if (is_simulating) {
        for (int i = 0; i < config.steps_per_frame; i++) {
            // Simulate one step here
            env->simulate(1.0, Vector3D(0, -1, 0));
            // ropeEuler->simulateEuler(1 / config.steps_per_frame, config.gravity);
            // ropeVerlet->simulateVerlet(1 / config.steps_per_frame, config.gravity);
        }
    }
    // Draw to screen here
    
    
    glBegin(GL_POINTS);
    for (int y = 0; y < config.nx_cells; y++) {
        for (int x = 0; x < config.ny_cells; x++) {
            glColor3f(env->ux[x+y*config.nx_cells], env->uy[x+y*config.nx_cells], 0.0);
            glVertex2d(x, y);
        }
    }
    glEnd();
    glFlush();
}

void Application::resize(size_t w, size_t h) {
  screen_width = w;
  screen_height = h;

  float half_width = (float)screen_width / 2;
  float half_height = (float)screen_height / 2;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, screen_width, 0, screen_height, 1, 0);
}

void Application::keyboard_event(int key, int event, unsigned char mods) {
    switch (key) {
    case '-':
        if (config.steps_per_frame > 1) {
            config.steps_per_frame /= 2;
        }
        break;
    case '=':
        config.steps_per_frame *= 2;
        break;
    case 'p':
        is_simulating = !is_simulating;
        break;
    }
}

string Application::name() { return "Powder"; }

string Application::info() {
  ostringstream steps;
  steps << "Steps per frame: " << config.steps_per_frame;

  return steps.str();
}
}
