#pragma once

#include <random>
#include "particlefilter.h"

class UnoptimizedParticleFilter : public ParticleFilter
{
public:
    UnoptimizedParticleFilter(unsigned int nparticles, unsigned int screen_height, unsigned int screen_width);
    ~UnoptimizedParticleFilter(void);

    /** Override from parent class */
    unsigned int get_nparticles(void) const;
    /** Override from parent class */
    int get_xpos(unsigned int index) const;
    /** Override from parent class */
    int get_ypos(unsigned int index) const;

    /** Override from parent class */
    void update(Robot &robot);

private:
    std::mt19937 rng;
    std::uniform_int_distribution<std::mt19937::result_type> height_distribution;
    std::uniform_int_distribution<std::mt19937::result_type> width_distribution;
    unsigned int *indices;

    /** Calculates P(position | measurement) for the particle at index i. */
    double calculate_likelihood(unsigned int i, int measured_x, int measured_y) const;

    /** Return a value drawn from the Gaussian distribution parameterized by mean and sigma. */
    double gaussian_noise(double mean, double sigma);

    /** Resample the particles from the distribution created from the weights. */
    void resample_particles(void);

    /** Sort the particles by weight. Sorts in place. */
    void sort_particles_by_weight_in_place(void);
};

/** Tiny utility class for sorting the particles' indices */
class SortIndices
{
   private:
     int *arr;

   public:
     SortIndices(int *arr) : arr(arr) {}
     bool operator()(int i, int j) const { return arr[i] < arr[j]; }
};
