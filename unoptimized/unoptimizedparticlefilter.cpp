#include <random>
#include "robot.h"
#include "unoptimizedparticlefilter.h"

UnoptimizedParticleFilter::UnoptimizedParticleFilter(unsigned int nparticles, unsigned int screen_height, unsigned int screen_width)
    : ParticleFilter(nparticles)
{
    std::random_device dev;
    this->rng = std::mt19937(dev());
    this->height_distribution = std::uniform_int_distribution<std::mt19937::result_type>(0, screen_height);
    this->width_distribution = std::uniform_int_distribution<std::mt19937::result_type>(0, screen_width);

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
   // TODO

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

double UnoptimizedParticleFilter::gaussian_noise(double mean, double sigma) const
{
    // TODO
    /*
        (1 / sqrt(2*pi*sigma*sigma)) * exp(-1 * (x - mu)^2 / 2*sigma*sigma)
    */
   return (double)(mean * sigma);
}
