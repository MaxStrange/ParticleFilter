#include <algorithm>
#include <random>
#include <stdlib.h>
#include "robot.h"
#include "unoptimizedparticlefilter.h"

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

void UnoptimizedParticleFilter::update(Robot &robot)
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
   for (unsigned int i = 0; i < this->nparticles; i++)
   {
       this->particles_weights[i] = this->calculate_likelihood(i, estimate_x, estimate_y);
   }

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
       this->particles_x[i] += this->gaussian_noise(estimate_vx, sigma);
       this->particles_y[i] += this->gaussian_noise(estimate_vy, sigma);
       // TODO: Update weight
   }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// PRIVATE FUNCTIONS ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////

double UnoptimizedParticleFilter::calculate_likelihood(unsigned int i, int measured_x, int measured_y) const
{
    // TODO
    /*
        P(A | B) = P(B | A) * P(A)   /  P(B)
        P(location | measurement) = P(measurement | location) * P(location) / P(measurement)
    */
   return (double)(i * measured_x * measured_y);
}

double UnoptimizedParticleFilter::gaussian_noise(double mean, double sigma)
{
   std::normal_distribution<double> gaussian(mean, sigma);
   return gaussian(this->rng);
}

void UnoptimizedParticleFilter::resample_particles(void)
{
    // Create a distribution I will need
    auto dist = std::uniform_real_distribution<std::mt19937::result_type>(0.0, 1.0);

    // Create the new particles in vectors
    std::vector<int> pxs;
    std::vector<int> pys;

    // Sort the particles by weight
    this->sort_particles_by_weight_in_place();

    // Align a CMF (cumulative mass function) array, where each bin is the sum of all previous weights
    std::vector<double> cmf;
    double acc_prob_mass = 0.0;
    for (int i = this->nparticles - 1; i >= 0; i--)
    {
        acc_prob_mass += this->particles_weights[i];
        cmf.push_back(acc_prob_mass);
    }

    // Do a search into the CMF to find the place where our randomly generated probability (0 to 1) fits
    for (unsigned int i = 0; i < this->nparticles; i++)
    {
        double p = dist(this->rng);
        unsigned int cmf_index = 0;
        for (unsigned int j = 0; j < this->nparticles; j++)
        {
            // Search for where the generated probability belongs
            if (p >= cmf[j])
            {
                cmf_index = j;
            }
        }
        pxs.push_back(this->particles_x[cmf_index]);
        pys.push_back(this->particles_y[cmf_index]);
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
    std::sort(this->indices, this->indices + this->nparticles, SortIndices(this->indices));

    // Sort each array according to the sorted indices
}
