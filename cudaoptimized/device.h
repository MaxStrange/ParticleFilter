#pragma once

#include <random>

/**  Updates `weights` according to the likelihood of each particle, given their position and the given estimate. */
void device_calculate_likelihood(const int *particles_x, const int *particles_y, int estimate_x, int estimate_y, double *weights, unsigned int nparticles);

void device_resample_and_move(int estimated_vx, int estimated_vy, unsigned int nparticles, int *particles_x, int *particles_y, double *particles_weights, std::mt19937 &rng, unsigned int *indices);

/** Tiny utility class for sorting the particles' indices */
class SortIndices
{
   private:
     double *weights;

   public:
     SortIndices(double *weights) : weights(weights) {}
     bool operator()(unsigned int i, unsigned int j) const { return weights[i] > weights[j]; }
};
