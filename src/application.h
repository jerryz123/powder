#ifndef CGL_APPLICATION_H
#define CGL_APPLICATION_H

// STL
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// libCGL
#include "CGL/CGL.h"
#include "CGL/osdtext.h"
#include "CGL/renderer.h"

#include "environment.h"


using namespace std;

namespace CGL {

    struct AppConfig {
        AppConfig() {
            // Assign globla vars here
            /* // Rope config variables */
            /* mass = 1; */
            /* ks = 100; */

            // Environment variables
            steps_per_frame = 1;
            cell_width = 1.0;
            cell_height = 1.0;
            nx_cells = 600;
            ny_cells = 400;
            is_simulating = true;
            mode = temperature;
            delta_t = 0.001;

        }

        /* Define global vars here */
        bool is_simulating;
        size_t nx_cells;
        size_t ny_cells;
        float steps_per_frame;
        float delta_t;
        float cell_width;
        float cell_height;
        InputMode mode;

    };

    class Application : public Renderer {
    public:
        Application(AppConfig config);
        ~Application();

        void init();
        void render();
        void resize(size_t w, size_t h);

        std::string name();
        std::string info();

        void key_event(char key);
        // void cursor_event(float x, float y);
        // void scroll_event(float offset_x, float offset_y);
        void cursor_event(float x, float y, unsigned char keys);
        float interp(int i0, int i1, int j0, int j1, float* arr);

    private:
        AppConfig config;

        /* Track environment here */
        /* Rope *ropeEuler; */
        /* Rope *ropeVerlet; */
        bool is_simulating;
        size_t screen_width;
        size_t screen_height;
        Environment* env;
        vector<InputItem> inputs;



    }; // class Application

} // namespace CGL

#endif // CGL_APPLICATION_H
