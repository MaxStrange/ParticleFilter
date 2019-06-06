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

    /** Update the particle filter. Call this every time we get new measurements from the robot. */
    virtual void update(const Robot &robot) = 0;

protected:
    /** Particle's X locations. Malloc'd array. */
    int *particles_x;

    /** Particle's Y locations. Malloc'd array. */
    int *particles_y;

    /** Particle's weights. Malloc'd array. */
    double *particles_weights;
};
