#pragma once

#include "robot.h"

/**
 * This is a base class for all particle filters - unoptimized and CUDA-optimized.
 */
class ParticleFilter
{
public:
    ParticleFilter(unsigned int nparticles);
    virtual ~ParticleFilter(void);

    virtual unsigned int get_nparticles(void) const = 0;
    virtual int get_xpos(unsigned int index) const = 0;
    virtual int get_ypos(unsigned int index) const = 0;

    /**
     * Update the particle filter. Call this every time we get new measurements from the robot.
     * Then display, then call part2.
     */
    virtual void update_part1(Robot &robot) = 0;
    virtual void update_part2(Robot &robot) = 0;

protected:
    /** Particle's X locations. Malloc'd array. */
    int *particles_x;

    /** Particle's Y locations. Malloc'd array. */
    int *particles_y;

    /** Particle's weights. Malloc'd array. */
    double *particles_weights;

    /** Number of particles in the filter. */
    unsigned int nparticles;
};
