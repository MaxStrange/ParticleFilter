#pragma once

#include "particlefilter.h"

class UnoptimizedParticleFilter : public ParticleFilter
{
public:
    UnoptimizedParticleFilter(unsigned int nparticles);
    ~UnoptimizedParticleFilter(void);

    /** Override from parent class */
    void update(const Robot &robot);

private:
};
