#pragma once

#include <random>

/**  Updates `weights` according to the likelihood of each particle, given their position and the given estimate. */
int device_calculate_likelihood(const int *particles_x, const int *particles_y, int estimate_x, int estimate_y, float *weights, unsigned int nparticles, int nthreads);

int device_resample_and_move(int estimated_vx, int estimated_vy, unsigned int nparticles, int *particles_x, int *particles_y, float *particles_weights, std::mt19937 &rng, unsigned int *indices, int nthreads);

/** Tiny utility class for sorting the particles' indices */
class SortIndices
{
   private:
     float *weights;

   public:
     SortIndices(float *weights) : weights(weights) {}
     bool operator()(unsigned int i, unsigned int j) const { return weights[i] > weights[j]; }
};
