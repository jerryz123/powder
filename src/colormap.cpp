#include <algorithm>
#include "colormap.h"

#define CLIP(x, y, z) (max(y, min(x, z)))

namespace CGL {
    Vector3D jet(float x) {
        float r, g, b;
        if (x < 0.7) {
            r =  4.0 * x - 1.5;
        } else {
            r = -4.0 * x + 4.5;
        }
        if (x < 0.5) {
            g = 4.0 * x - 0.5;
        } else {
            g = -4.0 * x + 3.5;
        }
        if (x < 0.3) {
            b = 4.0 * x + 0.5;
        } else {
            b = -4.0 * x + 2.5;
        }

        
        r = min(max(r, (float)0.), (float)1.);
        g = min(max(g, (float)0.), (float)1.);
        b = min(max(b, (float)0.), (float)1.);
        return Vector3D(r, g, b);
    }
    Vector3D hot(float x) {
        float r = CLIP(8.0 / 3.0 * x, 0.0, 1.0);
        float g = CLIP(8.0 / 3.0 * x - 1.0, 0.0, 1.0);
        float b = CLIP(4.0 * x - 3.0, 0.0, 1.0);
    
        return Vector3D(r, g, b);
    }
    
};
