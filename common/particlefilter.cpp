#include <stdlib.h>
#include "particlefilter.h"
#include "robot.h"

/*
    Particle filter algorithm is like this:

    1. Sample uniformly from all posible locations of the robot
    2. Get (noisy) measurement of robot's position from the robot's sensors
    3. Assign weight to each particle based on how well that particle matches
       the measurement.
    4. Resample the particles according to their weights.
    5. Reset all weights.
    6. Move all particles according to a model of how the robot moves. Do this
       from a Gaussian centered on our estimate of motion, and update weights
       so that the ones that are unlikely based on our understanding of how the
       robot likely moved have small weights.
    7. Repeat from step 2.
 */

ParticleFilter::ParticleFilter(unsigned int nparticles)
{
    this->particles_x = (int *)malloc(sizeof(int) * nparticles);
    this->particles_y = (int *)malloc(sizeof(int) * nparticles);
    this->particles_weights = (double *)malloc(sizeof(double) * nparticles);
}

ParticleFilter::~ParticleFilter(void)
{
    if (this->particles_x != nullptr)
    {
        free(this->particles_x);
        this->particles_x = nullptr;
    }

    if (this->particles_y != nullptr)
    {
        free(this->particles_y);
        this->particles_y = nullptr;
    }

    if (this->particles_weights != nullptr)
    {
        free(this->particles_weights);
        this->particles_weights = nullptr;
    }
}
