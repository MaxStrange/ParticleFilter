#include <algorithm>
#include <assert.h>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "robot.h"
#include "unoptimizedparticlefilter.h"

#define DEBUG 1
#define DEBUG_PRINTF(fmt, ...) do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)

UnoptimizedParticleFilter::UnoptimizedParticleFilter(unsigned int nparticles, unsigned int screen_height, unsigned int screen_width)
    : ParticleFilter(nparticles)
{
    std::random_device dev;
    this->rng = std::mt19937(dev());
    this->height_distribution = std::uniform_int_distribution<std::mt19937::result_type>(0, screen_height);
    this->width_distribution = std::uniform_int_distribution<std::mt19937::result_type>(0, screen_width);

    // Allocate the indices
    this->indices = (unsigned int *)malloc(sizeof(unsigned int) * nparticles);

    // Initialize the particles
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        int x = this->width_distribution(this->rng);
        int y = this->height_distribution(this->rng);

        this->particles_x[i] = x;
        this->particles_y[i] = y;
        this->particles_weights[i] = 0.0;
        this->indices[i] = i;

        assert(x <= (int)screen_width);
        assert(y <= (int)screen_height);
    }
}

UnoptimizedParticleFilter::~UnoptimizedParticleFilter(void)
{
    // Parent's destructor called automatically

    // Free the indices we allocated in our own constructor
    if (this->indices != nullptr)
    {
        free(this->indices);
        this->indices = nullptr;
    }
}

unsigned int UnoptimizedParticleFilter::get_nparticles(void) const
{
    return this->nparticles;
}

int UnoptimizedParticleFilter::get_xpos(unsigned int index) const
{
    return this->particles_x[index];
}

int UnoptimizedParticleFilter::get_ypos(unsigned int index) const
{
    return this->particles_y[index];
}

void UnoptimizedParticleFilter::update_part1(Robot &robot)
{
    /*
        Particle filter algorithm is like this:

        1. Sample uniformly from all posible locations of the robot
            --- Already done in the initialization of this object

        2. Get (noisy) measurement of robot's position from the robot's sensors

        3. Assign weight to each particle based on how well that particle matches
        the measurement.

        4. Resample the particles according to their weights.

        5. Reset all weights.

        6. Move all particles according to a model of how the robot moves. Do this
        from a Gaussian centered on our estimate of motion, and update weights
        so that the ones that are unlikely based on our understanding of how the
        robot likely moved have small weights.

        7. Repeat from step 2.
    */

    // Get measurements of robot's position from its sensors
    int estimate_x;
    int estimate_y;
    robot.get_xy_estimate(&estimate_x, &estimate_y);

    // Assign the weights
    double accumulated = 0.0;
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        this->particles_weights[i] = this->calculate_likelihood(i, estimate_x, estimate_y);
        accumulated += this->particles_weights[i];
    }
    assert(accumulated > 0.0);

    // Take a break now to update the graphics in the main loop
}

void UnoptimizedParticleFilter::update_part2(Robot &robot)
{
    // Resample from weights
    this->resample_particles();

    // Reset all weights
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        this->particles_weights[i] = 0.0;
    }

    // Move all particles according to our movement model (plus Gaussian noise)
    // Also update weights based on how likely the movements were
    int estimate_vx;
    int estimate_vy;
    robot.get_v_estimate(&estimate_vx, &estimate_vy);
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        static const double sigma = 2.5;
        double vx = this->gaussian_noise(estimate_vx, sigma);
        double vy = this->gaussian_noise(estimate_vy, sigma);
        this->particles_x[i] += vx;
        this->particles_y[i] += vy;
        this->particles_weights[i] = this->probability_of_value_from_bivariate_gaussian(vx, vy, estimate_vx, estimate_vy, sigma, sigma);
    }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// PRIVATE FUNCTIONS ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////

double UnoptimizedParticleFilter::calculate_likelihood(unsigned int i, int measured_x, int measured_y) const
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
    double x = (double)this->particles_x[i];
    double y = (double)this->particles_y[i];

    return this->probability_of_value_from_bivariate_gaussian(x, y, measured_x, measured_y, 2.5, 2.5);
}

double UnoptimizedParticleFilter::gaussian_noise(double mean, double sigma)
{
    std::normal_distribution<double> gaussian(mean, sigma);
    return gaussian(this->rng);
}

void UnoptimizedParticleFilter::normalize_weights(void)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        sum += this->particles_weights[i];
    }

    if (sum > 0.0)
    {
        for (unsigned int i = 0; i < this->nparticles; i++)
        {
            this->particles_weights[i] /= sum;
            assert((this->particles_weights[i] >= 0.0) && (this->particles_weights[i] <= 1.0));
        }
    }
}

double UnoptimizedParticleFilter::probability_of_value_from_bivariate_gaussian(double x, double y, double mean_x, double mean_y, double sigma_x, double sigma_y) const
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

void UnoptimizedParticleFilter::resample_particles(void)
{
    // Create a distribution I will need
    auto dist = std::uniform_real_distribution<double>(0.0, 1.0);

    // Create the new particles in vectors
    std::vector<int> pxs;
    std::vector<int> pys;

    // Normalize the weights so that each one is between 0 and 1
    this->normalize_weights();

    printf("Before sorting: ");
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        printf("(%d, %d, %f) ", this->particles_x[i], this->particles_y[i], this->particles_weights[i]);
    }
    printf("\n");

    // Sort the particles by weight (in reverse - heaviest at the front of the array)
    this->sort_particles_by_weight_in_place();

    printf("After sorting: ");
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        printf("(%d, %d, %f) ", this->particles_x[i], this->particles_y[i], this->particles_weights[i]);
    }
    printf("\n");

    // Align a CMF (cumulative mass function) array, where each bin is the sum of all previous weights
    std::vector<double> cmf;
    double acc_prob_mass = 0.0;
    for (unsigned int i = 0; i < this->nparticles; i--)
    {
        acc_prob_mass += this->particles_weights[i];
        cmf.push_back(acc_prob_mass);
    }

    exit(0);

    // Do a search into the CMF to find the place where our randomly generated probability (0 to 1) fits
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        double p = dist(this->rng);
        int cmf_index = -1;
        for (unsigned int j = 0; j < this->nparticles; j++)
        {
            // Search for where the generated probability belongs
            if (p >= cmf[j])
            {
                cmf_index = j;
            }
        }

        if (cmf_index >= 0)
        {
            pxs.push_back(this->particles_x[cmf_index]);
            pys.push_back(this->particles_y[cmf_index]);
        }
        else
        {
            // Probabilities are all zero. Resample from uniform.
            pxs.push_back(this->width_distribution(this->rng));
            pys.push_back(this->height_distribution(this->rng));
        }
    }

    // Now overwrite the current batch of particles with the new ones
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        this->particles_x[i] = pxs[i];
        this->particles_y[i] = pys[i];
    }
}

void UnoptimizedParticleFilter::sort_particles_by_weight_in_place(void)
{
    // Sort the indices
    std::sort(this->indices, this->indices + this->nparticles, SortIndices(this->particles_weights));

    // Make copies of the three arrays (gross)
    int *xcpy = (int *)malloc(sizeof(int) * this->nparticles);
    int *ycpy = (int *)malloc(sizeof(int) * this->nparticles);
    double *wcpy = (double *)malloc(sizeof(double) * this->nparticles);
    memcpy(xcpy, this->particles_x, sizeof(int) * this->nparticles);
    memcpy(ycpy, this->particles_y, sizeof(int) * this->nparticles);
    memcpy(wcpy, this->particles_weights, sizeof(double) * this->nparticles);

    // Sort each array according to the sorted indices
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        this->particles_weights[i] = wcpy[this->indices[i]];
        this->particles_x[i] = xcpy[this->indices[i]];
        this->particles_y[i] = ycpy[this->indices[i]];
    }

   free(xcpy);
   free(ycpy);
   free(wcpy);
   xcpy = nullptr;
   ycpy = nullptr;
   wcpy = nullptr;
}
