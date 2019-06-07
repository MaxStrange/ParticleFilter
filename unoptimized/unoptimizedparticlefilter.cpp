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
   // TODO

   // Resample from weights
   // TODO

   // Reset all weights
   for (unsigned int i = 0; i < this->nparticles; i++)
   {
       this->particles_weights[i] = 0.0;
   }

   // Move all particles according to our movement model (plus Gaussian noise)
   // TODO

   // Update weights based on how likely the movements were
   // TODO
}
