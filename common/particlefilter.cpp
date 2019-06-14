#include <stdlib.h>
#include "particlefilter.h"
#include "robot.h"

ParticleFilter::ParticleFilter(unsigned int nparticles)
{
    this->particles_x = (int *)malloc(sizeof(int) * nparticles);
    this->particles_y = (int *)malloc(sizeof(int) * nparticles);
    this->particles_weights = (float *)malloc(sizeof(float) * nparticles);
    this->nparticles = nparticles;
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
