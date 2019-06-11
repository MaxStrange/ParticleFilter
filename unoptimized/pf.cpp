/**
 * Much of this has been stolen most graciously from lazyfoo's SDL tutorial:
 * https://lazyfoo.net/tutorials/SDL/
 */
#include <iostream>
#include <string>
#include "graphics.h"
#include "robot.h"
#include "unoptimizedparticlefilter.h"

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

static int parse_args(int argc, char *argv[], unsigned int *nparticles)
{
    if (argc == 1)
    {
        return -1;
    }

    int n = atoi(argv[1]);
    if (n > 0)
    {
        *nparticles = n;
        return 0;
    }
    else
    {
        nparticles = nullptr;
        return -1;
    }
}

int main(int argc, char *argv[])
{
    int err;
    Graphics gfx;

    unsigned int nparticles;
    err = parse_args(argc, argv, &nparticles);
    CHECK_EXIT(err, "Need a number of particles, must be greater than zero.", gfx);

    UnoptimizedParticleFilter pf(nparticles, gfx.get_screenheight(), gfx.get_screenwidth());
    Robot robot;

    err = gfx.init();
    CHECK_EXIT(err, "SDL could not initialize!", gfx);

    bool done = false;
    while (!done)
    {
        // Derive new state of particles (using particle filter)
        pf.update_part1(robot);

        // Update the graphics
        gfx.update(robot, pf);

        // Do the rest of the particle filter algorithm
        pf.update_part2(robot);

        // Get user input
        done = gfx.isdone();

        // Derive new state of robot
        robot.update(gfx.get_screenheight(), gfx.get_screenwidth());
    }

    gfx.exit();

    return 0;
}
