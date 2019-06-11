#pragma once

#include <random>
#include "particlefilter.h"

class CudaParticleFilter: public ParticleFilter
{
public:
    CudaParticleFilter(unsigned int nparticles, unsigned int screen_height, unsigned int screen_width);
    ~CudaParticleFilter(void);

    /** Override from parent class */
    unsigned int get_nparticles(void) const;
    /** Override from parent class */
    int get_xpos(unsigned int index) const;
    /** Override from parent class */
    int get_ypos(unsigned int index) const;

    /** Override from parent class */
    void update_part1(Robot &robot);
    void update_part2(Robot &robot);

private:
    std::mt19937 rng;
    std::uniform_int_distribution<std::mt19937::result_type> height_distribution;
    std::uniform_int_distribution<std::mt19937::result_type> width_distribution;
    unsigned int *indices;
};
