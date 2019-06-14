#include <algorithm>
#include <assert.h>
#include <cuda_runtime.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <vector>
#include "device.h"

static float gaussian_noise(float mean, float sigma, std::mt19937 &rng)
{
    std::normal_distribution<float> gaussian(mean, sigma);
    return gaussian(rng);
}

static float probability_of_value_from_bivariate_gaussian(float x, float y, float mean_x, float mean_y, float sigma_x, float sigma_y)
{
    const float rho = 0.0; // cov / (sig1 * sig2); Covariance of two independent random variables is zero.
    float denom = 2.0 * M_PI * sigma_x * sigma_y * sqrt(1.0 - (rho * rho));
    float A = ((x - mean_x) * (x - mean_x)) / (sigma_x * sigma_x);
    float B = ((2.0 * rho * (x - mean_x) * (y - mean_y)) / (sigma_x * sigma_y));
    float C = ((y - mean_y) * (y - mean_y)) / (sigma_y * sigma_y);
    A /= 1000.0;  // For numerical stability
    C /= 1000.0;  // Ditto
    float z = A - B + C;
    float a = (-1.0 * z) / (2.0 * (1.0 - rho * rho));

    return exp(a) / denom;
}

__global__ void kernel_calculate_likelihood(int *particles_x, int *particles_y, float *weights, unsigned int nparticles, float estimate_x, float estimate_y)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < nparticles)
    {
        float x = particles_x[index];
        float y = particles_y[index];

        const float sigma_x = 2.5;
        const float sigma_y = 2.5;
        float mean_x = estimate_x;
        float mean_y = estimate_y;

        // Compute the probability of getting this x,y combo from a distribution centered at estimate_x, estimte_y.
        const float rho = 0.0; // cov / (sig1 * sig2); Covariance of two independent random variables is zero.
        float denom = 2.0 * M_PI * sigma_x * sigma_y * sqrt(1.0 - (rho * rho));
        float A = ((x - mean_x) * (x - mean_x)) / (sigma_x * sigma_x);
        float B = ((2.0 * rho * (x - mean_x) * (y - mean_y)) / (sigma_x * sigma_y));
        float C = ((y - mean_y) * (y - mean_y)) / (sigma_y * sigma_y);
        A /= 1000.0;  // For numerical stability
        C /= 1000.0;  // Ditto
        float z = A - B + C;
        float a = (-1.0 * z) / (2.0 * (1.0 - rho * rho));

        float prob = exp(a) / denom;

        weights[index] = prob;
    }
}

int device_calculate_likelihood(const int *particles_x, const int *particles_y, int estimate_x, int estimate_y, float *weights, unsigned int nparticles)
{
    cudaError_t err;
    int *dev_particles_x = nullptr;
    int *dev_particles_y = nullptr;
    float *dev_weights = nullptr;

    unsigned int i;

#define CHECK_CUDA_ERR(err) do { if (err != cudaSuccess) { err = (cudaError_t)__LINE__; goto fail; }} while (0)

    /* Malloc all the device memory we need */
    err = cudaMalloc(&dev_particles_x, nparticles * sizeof(int));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_particles_y, nparticles * sizeof(int));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_weights, nparticles * sizeof(float));
    CHECK_CUDA_ERR(err);

    /* Copy arrays onto device */
    err = cudaMemcpy(dev_particles_x, particles_x, nparticles * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(dev_particles_y, particles_y, nparticles * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(dev_weights, weights, nparticles * sizeof(float), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERR(err);

    /* Call the kernel */
    //kernel_calculate_likelihood<<<ceil(nparticles / 524.0), 524>>>(dev_particles_x, dev_particles_y, dev_weights, estimate_x, estimate_y, nparticles);
    for (i = 0; i < nparticles; i++)
    {
        float x = (float)particles_x[i];
        float y = (float)particles_y[i];

        weights[i] = probability_of_value_from_bivariate_gaussian(x, y, estimate_x, estimate_y, 2.5, 2.5);
    }

    /* Copy array back onto host */
    err = cudaMemcpy(weights, dev_weights, nparticles * sizeof(float), cudaMemcpyDeviceToHost);
    CHECK_CUDA_ERR(err);

    /* Deallocate the device arrays */
    err = cudaFree(dev_particles_x);
    CHECK_CUDA_ERR(err);

    err = cudaFree(dev_particles_y);
    CHECK_CUDA_ERR(err);

    err = cudaFree(dev_weights);
    CHECK_CUDA_ERR(err);

#undef CHECK_CUDA_ERR

fail:
    assert(err == cudaSuccess);
    return (int)err;
}

static void move_particles(int estimated_vx, int estimated_vy, unsigned int nparticles, int *particles_x, int *particles_y, float *particles_weights, std::mt19937 &rng)
{
    for (unsigned int i = 0; i < nparticles; i++)
    {
        static const float sigma = 2.5;
        float vx = gaussian_noise(estimated_vx, sigma, rng);
        float vy = gaussian_noise(estimated_vy, sigma, rng);
        particles_x[i] += vx;
        particles_y[i] += vy;
        particles_weights[i] = probability_of_value_from_bivariate_gaussian(vx, vy, estimated_vx, estimated_vy, sigma, sigma);
    }
}

static void sort_particles_by_weight_in_place(unsigned int *indices, unsigned int nparticles, float *particles_weights, int *particles_x, int *particles_y)
{
    // Sort the indices
    std::sort(indices, indices + nparticles, SortIndices(particles_weights));

    // Make copies of the three arrays (gross)
    int *xcpy = (int *)malloc(sizeof(int) * nparticles);
    int *ycpy = (int *)malloc(sizeof(int) * nparticles);
    float *wcpy = (float *)malloc(sizeof(float) * nparticles);
    memcpy(xcpy, particles_x, sizeof(int) * nparticles);
    memcpy(ycpy, particles_y, sizeof(int) * nparticles);
    memcpy(wcpy, particles_weights, sizeof(float) * nparticles);

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

static void normalize_weights(unsigned int nparticles, float *particles_weights)
{
    float sum = 0.0;
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

static void resample_particles(unsigned int nparticles, float *particles_weights, std::mt19937 &rng, unsigned int *indices, int *particles_x, int *particles_y)
{
    // Create a distribution I will need
    auto dist = std::uniform_real_distribution<float>(0.0, 1.0);
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
    std::vector<float> cmf;
    float acc_prob_mass = 0.0;
    for (unsigned int i = 0; i < nparticles; i++)
    {
        acc_prob_mass += particles_weights[i];
        cmf.push_back(acc_prob_mass);
    }

    // Do a search into the CMF to find the place where our randomly generated probability (0 to 1) fits
    for (unsigned int i = 0; i < nparticles; i++)
    {
        float p = dist(rng);
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

int device_resample_and_move(int estimated_vx, int estimated_vy, unsigned int nparticles, int *particles_x, int *particles_y, float *particles_weights, std::mt19937 &rng, unsigned int *indices)
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

    return 0;
}
