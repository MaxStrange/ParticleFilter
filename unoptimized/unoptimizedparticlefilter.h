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
    void update_part1(Robot &robot);
    void update_part2(Robot &robot);

private:
    std::mt19937 rng;
    std::uniform_int_distribution<std::mt19937::result_type> height_distribution;
    std::uniform_int_distribution<std::mt19937::result_type> width_distribution;
    unsigned int *indices;

    /** Calculates P(position | measurement) for the particle at index i. */
    float calculate_likelihood(unsigned int i, int measured_x, int measured_y) const;

    /** Return a value drawn from the Gaussian distribution parameterized by mean and sigma. */
    float gaussian_noise(float mean, float sigma);

    /** Normalize the weights to between 0 and 1. */
    void normalize_weights(void);

    /**
     * Return the probability of x and y being drawn from the bivariate distribution described by mean and sigma.
     * Implementation currently assumes that the bivariate is composed of two (possibly different) Gaussians,
     * which are *independent*. If they are correlated, this will not work.
     */
    float probability_of_value_from_bivariate_gaussian(float x, float y, float mean_x, float mean_y, float sigma_x, float sigma_y) const;

    /** Resample the particles from the distribution created from the weights. */
    void resample_particles(void);

    /** Sort the particles by weight. Sorts in place. */
    void sort_particles_by_weight_in_place(void);
};

/** Tiny utility class for sorting the particles' indices */
class SortIndices
{
   private:
     float *weights;

   public:
     SortIndices(float *weights) : weights(weights) {}
     bool operator()(unsigned int i, unsigned int j) const { return weights[i] > weights[j]; }
};
