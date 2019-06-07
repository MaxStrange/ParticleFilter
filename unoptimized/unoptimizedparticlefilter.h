#pragma once

#include <random>
#include "particlefilter.h"

class UnoptimizedParticleFilter : public ParticleFilter
{
public:
    UnoptimizedParticleFilter(unsigned int nparticles, unsigned int screen_height, unsigned int screen_width);
    ~UnoptimizedParticleFilter(void);

    /** Override from parent class */
    void update(Robot &robot);

private:
    std::mt19937 rng;
    std::uniform_int_distribution<std::mt19937::result_type> height_distribution;
    std::uniform_int_distribution<std::mt19937::result_type> width_distribution;
};
