/**
 * Much of this has been stolen most graciously from lazyfoo's SDL tutorial:
 * https://lazyfoo.net/tutorials/SDL/
 */
#include <iostream>
#include <string>
#include <unistd.h>
#include "graphics.h"
#include "robot.h"
#include "unoptimizedparticlefilter.h"

/** The number of particles for this test */
#define NPARTICLES 100

/** Exit macro: If err, print errormsg and exit cleanly. */
#define CHECK_EXIT(err, errormsg, g) do \
    {\
        if (err)\
        {\
            std::cout << errormsg << " Error: " << err << "\n";\
            g.exit();\
            exit(-1);\
        }\
    } while(0)

int main(void)
{
    int err;
    Graphics gfx;
    UnoptimizedParticleFilter pf(NPARTICLES, gfx.get_screenheight(), gfx.get_screenwidth());
    Robot robot;

    err = gfx.init();
    CHECK_EXIT(err, "SDL could not initialize!", gfx);

    bool done = false;
    while (!done)
    {
        // Update the graphics
        gfx.update(robot);

        // Get user input
        done = gfx.isdone();

        // Derive new state of robot
        robot.update(gfx.get_screenheight(), gfx.get_screenwidth());

        // Derive new state of particles (using particle filter)
        pf.update(robot);
    }

    gfx.exit();

    return 0;
}
