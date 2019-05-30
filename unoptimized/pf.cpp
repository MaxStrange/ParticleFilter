/**
 * Much of this has been stolen most graciously from lazyfoo's SDL tutorial:
 * https://lazyfoo.net/tutorials/SDL/
 */
#include <iostream>
#include "graphics.h"

/** Exit macro: If err, print errormsg and exit cleanly. */
#define CHECK_EXIT(err, errormsg, g) do \
    {\
        if (err)\
        {\
            std::cout << errormsg << "\n";\
            g.exit();\
            exit(-1);\
        }\
    } while(0)

int main(void)
{
    int err;

    Graphics gfx;

    err = gfx.init();
    CHECK_EXIT(err, "SDL could not initialize!", gfx);

    bool done = false;
    while (!done)
    {
        // Update the graphics
        gfx.update();

        // Get user input
        done = gfx.isdone();

        // Derive new state (using particle filter)
    }

    gfx.exit();

    return 0;
}
