#include <algorithm>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <vector>
#include "device.h"

static double gaussian_noise(double mean, double sigma, std::mt19937 &rng)
{
    std::normal_distribution<double> gaussian(mean, sigma);
    return gaussian(rng);
}

static double probability_of_value_from_bivariate_gaussian(double x, double y, double mean_x, double mean_y, double sigma_x, double sigma_y)
{
    const double rho = 0.0; // cov / (sig1 * sig2); Covariance of two independent random variables is zero.
    double denom = 2.0 * M_PI * sigma_x * sigma_y * sqrt(1.0 - (rho * rho));
    double A = ((x - mean_x) * (x - mean_x)) / (sigma_x * sigma_x);
    double B = ((2.0 * rho * (x - mean_x) * (y - mean_y)) / (sigma_x * sigma_y));
    double C = ((y - mean_y) * (y - mean_y)) / (sigma_y * sigma_y);
    A /= 1000.0;  // For numerical stability
    C /= 1000.0;  // Ditto
    double z = A - B + C;
    double a = (-1.0 * z) / (2.0 * (1.0 - rho * rho));

    return exp(a) / denom;
}

void device_calculate_likelihood(const int *particles_x, const int *particles_y, int estimate_x, int estimate_y, double *weights, unsigned int nparticles)
{
    /*
        P(A | B) = P(B | A) * P(A)   /  P(B)
        P(location | measurement) = P(measurement | location) * P(location) / P(measurement)

        We can actually just use: P(measurement | lcoation) * P(location), discarding the probability of
        the measurement, since all particles are using the same measurement, and all I care about is
        the relative probability of each particle, not the true probability.

        So, using:
        P(measurement | location) = Gaussian(mean=location, std=who knows)
        P(location) = Uniform probability over the whole range (we don't have a reason to believe one location
                      is more probable in general than any other).

        Since P(location) is uniform, and therefore doesn't matter per particle, we can also do away with it,
        leaving just a Gaussian with mean centered on the location.

        So now we need the probability of measured_x and measured_y, given a Gaussian around location.
    */
    for (unsigned int i = 0; i < nparticles; i++)
    {
        double x = (double)particles_x[i];
        double y = (double)particles_y[i];
        weights[i] = probability_of_value_from_bivariate_gaussian(x, y, estimate_x, estimate_y, 2.5, 2.5);
    }
}

static void move_particles(int estimated_vx, int estimated_vy, unsigned int nparticles, int *particles_x, int *particles_y, double *particles_weights, std::mt19937 &rng)
{
    for (unsigned int i = 0; i < nparticles; i++)
    {
        static const double sigma = 2.5;
        double vx = gaussian_noise(estimated_vx, sigma, rng);
        double vy = gaussian_noise(estimated_vy, sigma, rng);
        particles_x[i] += vx;
        particles_y[i] += vy;
        particles_weights[i] = probability_of_value_from_bivariate_gaussian(vx, vy, estimated_vx, estimated_vy, sigma, sigma);
    }
}

static void sort_particles_by_weight_in_place(unsigned int *indices, unsigned int nparticles, double *particles_weights, int *particles_x, int *particles_y)
{
    // Sort the indices
    std::sort(indices, indices + nparticles, SortIndices(particles_weights));

    // Make copies of the three arrays (gross)
    int *xcpy = (int *)malloc(sizeof(int) * nparticles);
    int *ycpy = (int *)malloc(sizeof(int) * nparticles);
    double *wcpy = (double *)malloc(sizeof(double) * nparticles);
    memcpy(xcpy, particles_x, sizeof(int) * nparticles);
    memcpy(ycpy, particles_y, sizeof(int) * nparticles);
    memcpy(wcpy, particles_weights, sizeof(double) * nparticles);

    // Sort each array according to the sorted indices
    for (unsigned int i = 0; i < nparticles; i++)
    {
        particles_weights[i] = wcpy[indices[i]];
        particles_x[i] = xcpy[indices[i]];
        particles_y[i] = ycpy[indices[i]];
    }

   free(xcpy);
   free(ycpy);
   free(wcpy);
   xcpy = nullptr;
   ycpy = nullptr;
   wcpy = nullptr;
}

static void normalize_weights(unsigned int nparticles, double *particles_weights)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < nparticles; i++)
    {
        sum += particles_weights[i];
    }

    if (sum > 0.0)
    {
        for (unsigned int i = 0; i < nparticles; i++)
        {
            particles_weights[i] /= sum;
            assert((particles_weights[i] >= 0.0) && (particles_weights[i] <= 1.0));
        }
    }
}

static void resample_particles(unsigned int nparticles, double *particles_weights, std::mt19937 &rng, unsigned int *indices, int *particles_x, int *particles_y)
{
    // Create a distribution I will need
    auto dist = std::uniform_real_distribution<double>(0.0, 1.0);
    std::uniform_int_distribution<std::mt19937::result_type> height_distribution;
    std::uniform_int_distribution<std::mt19937::result_type> width_distribution;

    // Create the new particles in vectors
    std::vector<int> pxs;
    std::vector<int> pys;

    // Normalize the weights so that each one is between 0 and 1
    normalize_weights(nparticles, particles_weights);

    // Sort the particles by weight (in reverse - heaviest at the front of the array)
    sort_particles_by_weight_in_place(indices, nparticles, particles_weights, particles_x, particles_y);

    // Align a CMF (cumulative mass function) array, where each bin is the sum of all previous weights
    std::vector<double> cmf;
    double acc_prob_mass = 0.0;
    for (unsigned int i = 0; i < nparticles; i++)
    {
        acc_prob_mass += particles_weights[i];
        cmf.push_back(acc_prob_mass);
    }

    // Do a search into the CMF to find the place where our randomly generated probability (0 to 1) fits
    for (unsigned int i = 0; i < nparticles; i++)
    {
        double p = dist(rng);
        assert((p <= 1.0) && (p >= 0.0));

        int cmf_index = -1;
        for (unsigned int j = 0; j < nparticles; j++)
        {
            // Search for where the generated probability belongs
            if (p <= cmf[j])
            {
                cmf_index = j;
                break;
            }
        }

        if (cmf_index >= 0)
        {
            pxs.push_back(particles_x[cmf_index]);
            pys.push_back(particles_y[cmf_index]);
        }
        else
        {
            // Probabilities are all zero. Resample from uniform.
            pxs.push_back(width_distribution(rng));
            pys.push_back(height_distribution(rng));
        }
    }

    // Now overwrite the current batch of particles with the new ones
    for (unsigned int i = 0; i < nparticles; i++)
    {
        particles_x[i] = pxs[i];
        particles_y[i] = pys[i];
    }
}

void device_resample_and_move(int estimated_vx, int estimated_vy, unsigned int nparticles, int *particles_x, int *particles_y, double *particles_weights, std::mt19937 &rng, unsigned int *indices)
{
    // Resample from weights
    resample_particles(nparticles, particles_weights, rng, indices, particles_x, particles_y);

    // Reset all weights
    for (unsigned int i = 0; i < nparticles; i++)
    {
        particles_weights[i] = 0.0;
    }

    // Move all particles according to our movement model (plus Gaussian noise)
    // Also update weights based on how likely the movements were
    move_particles(estimated_vx, estimated_vy, nparticles, particles_x, particles_y, particles_weights, rng);
}
