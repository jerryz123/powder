#include <iostream>

#include "application.h"
#include "environment.h"
#include "colormap.h"

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
        glClearColor(1, 1, 1, 1);
        env = new Environment(config.nx_cells, config.ny_cells,
                              config.cell_width, config.cell_height);

    }

    void Application::render() {
        if (config.is_simulating) {
            for (int i = 0; i < config.steps_per_frame; i++) {
                // Simulate one step here
                env->simulate(config.delta_t, config.gravity, inputs);
                // ropeEuler->simulateEuler(1 / config.steps_per_frame, config.gravity);
                // ropeVerlet->simulateVerlet(1 / config.steps_per_frame, config.gravity);
                inputs.clear();
            }
        }
        // Draw to screen here


        glBegin(GL_POINTS);
        float t;
        float t1;
        Vector3D color;
        for (int y = 0; y < config.ny_cells; y++) {
            for (int x = 0; x < config.nx_cells; x++) {
                switch (config.mode) {
                case temperature:
                    t = env->T[x+y*config.nx_cells];
                    color = jet(t);
                    glColor4f(color.x, color.y, color.z, 1);
                    break;
                case smoke:
                    t = env->smoke[x+y*config.nx_cells];
                    glColor4f(0, 0, 0, t);
                    break;
                case debug:
                    t = env->ux[x+y*config.nx_cells] + 0.5;
                    t1 = env->uy[x+y*config.nx_cells] + 0.5;
                    glColor4f(t, t1, 0, 1);
                    break;
                }
                
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
        glOrtho(0, screen_width, screen_height, 0, 1, 0);
    }

    void Application::key_event(char key) {
        switch (key) {
        case '-':
            if (config.steps_per_frame > 1) {
                config.steps_per_frame /= 2;
            }
            break;
        case '=':
            config.steps_per_frame *= 2;
            break;
        case 'P':
            config.is_simulating = !config.is_simulating;
            break;
        case 'S':
            config.mode = smoke;
            break;
        case 'T':
            config.mode = temperature;
            break;
        case 'D':
            config.mode = debug;
            break;
        }
    }

    void Application::cursor_event(float x, float y, unsigned char keys) {
        int left = keys >> 2;
        if (left > 0) {
            InputItem input;
            input.pos = Vector2D(x, y);
            input.input_mode = config.mode;
            inputs.push_back(input);
        }
    }



    string Application::name() { return "Powder"; }

    string Application::info() {
        ostringstream steps;
        steps << "Steps per frame: " << config.steps_per_frame;

        return steps.str();
    }
}
