#include "robot.h"
#include "unoptimizedparticlefilter.h"

UnoptimizedParticleFilter::UnoptimizedParticleFilter(unsigned int nparticles)
    : ParticleFilter(nparticles)
{
}

UnoptimizedParticleFilter::~UnoptimizedParticleFilter(void)
{
    // Parent's destructor called automatically
}

void UnoptimizedParticleFilter::update(const Robot &robot)
{
    // TODO
}
