#include <algorithm>
#include <assert.h>
#include <cuda_runtime.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <stdio.h>
#include <vector>
#include "device.h"

/** Good ol' MIN macro */
#define MIN(a, b)   ((a) < (b) ? (a) : (b))

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
        float x = (float)particles_x[index];
        float y = (float)particles_y[index];

        const float sigma_x = 2.5;
        const float sigma_y = 2.5;
        float mean_x = estimate_x;
        float mean_y = estimate_y;

        // Compute the probability of getting this x,y combo from a distribution centered at estimate_x, estimte_y.
        const float rho = 0.0; // cov / (sig1 * sig2); Covariance of two independent random variables is zero.
        float denom = 2.0f * M_PI * sigma_x * sigma_y * sqrt(1.0f - (rho * rho));
        float A = ((x - mean_x) * (x - mean_x)) / (sigma_x * sigma_x);
        float B = ((2.0f * rho * (x - mean_x) * (y - mean_y)) / (sigma_x * sigma_y));
        float C = ((y - mean_y) * (y - mean_y)) / (sigma_y * sigma_y);
        A /= 1000.0f;  // For numerical stability
        C /= 1000.0f;  // Ditto
        float z = A - B + C;
        float a = (-1.0f * z) / (2.0f * (1.0f - rho * rho));
        float prob = exp(a) / denom;
        weights[index] = prob;
    }
}

int device_calculate_likelihood(const int *particles_x, const int *particles_y, int estimate_x, int estimate_y, float *weights, unsigned int nparticles, int nthreads_per_block)
{
    cudaError_t err;
    int *dev_particles_x = nullptr;
    int *dev_particles_y = nullptr;
    float *dev_weights = nullptr;

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
    kernel_calculate_likelihood<<<ceil(nparticles / (float)nthreads_per_block), nthreads_per_block>>>(dev_particles_x, dev_particles_y, dev_weights, nparticles, estimate_x, estimate_y);
    err = cudaDeviceSynchronize();
    CHECK_CUDA_ERR(err);

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

static void complete_resample_and_move_step(unsigned int nparticles, float *particles_weights, std::mt19937 &rng, unsigned int *indices, int *particles_x, int *particles_y, int estimated_vx, int estimated_vy)
{
    // Create a distribution I will need
    auto dist = std::uniform_real_distribution<float>(0.0, 1.0);
    std::uniform_int_distribution<std::mt19937::result_type> height_distribution;
    std::uniform_int_distribution<std::mt19937::result_type> width_distribution;

    // Create the new particles in vectors
    std::vector<int> pxs;
    std::vector<int> pys;

    // Sort the particles by weight (in reverse - heaviest at the front of the array)
    //sort_particles_by_weight_in_place(indices, nparticles, particles_weights, particles_x, particles_y);

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

    // Reset all weights
    for (unsigned int i = 0; i < nparticles; i++)
    {
        particles_weights[i] = 0.0;
    }

    // Move particles
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

__global__ void kernel_normalize_weights_reduction(unsigned int nparticles, float *dev_weights, float *intermediate)
{
    // Dynamically-sized shared memory buffer for the reduction (this should be no smaller than blockDim.x)
    extern __shared__ float tmp[];

    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

    // load all weights in this block into temp array
    if (index < nparticles)
    {
        tmp[threadIdx.x] = dev_weights[index];
    }
    __syncthreads();

    // Now do a binary sum tree to reduce to a single accumulated total weight
    for (unsigned int stride = 1; stride < nparticles; stride *= 2)
    {
        if ((index < nparticles) && (threadIdx.x >= stride))
        {
            tmp[threadIdx.x] += tmp[threadIdx.x - stride];
        }
        __syncthreads();
    }

    // Each block now needs to add its total to its index in intermediate
    // So determine which thread should do this, since we only need one
    // item from each block
    bool lastusefulthread;
    if (blockIdx.x == (gridDim.x - 1))
    {
        // If my block index is that of the final block, then I am
        // the thread responsible for the last useful item if
        // my index is that of the final particle
        lastusefulthread = (index == (nparticles - 1));
    }
    else
    {
        // If my block is not the final one, then I am
        // the thread responsible for the last useful item if
        // my index is that of the final item in this block
        lastusefulthread = (threadIdx.x == (blockDim.x - 1));
    }

    if (lastusefulthread)
    {
        intermediate[blockIdx.x] = tmp[threadIdx.x];
    }
}

__global__ void kernel_normalize_weights_complete(unsigned int nparticles, float *dev_weights, float summed_weights)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    // Divide all weights by sum in parallel
    if ((index < nparticles) && summed_weights > 0.0f)
    {
        dev_weights[index] /= summed_weights;
    }
}

__device__ void kernel_sequential_merge(float *tmpbuf_weights, float *weights_a, float *weights_b,
                                        int *tmpbuf_x,         int *x_a,         int *x_b,
                                        int *tmpbuf_y,         int *y_a,         int *y_b,
                                        unsigned int len_arr_a, unsigned int len_arr_b)
{
    // Sorts backwards (largest first)

    unsigned int i = 0;
    unsigned int j = 0;
    while ((i < len_arr_a) && (j < len_arr_b))
    {
        float wa = weights_a[i];
        float wb = weights_b[j];
        if (wa > wb)
        {
            tmpbuf_weights[i + j] = weights_a[i];
            tmpbuf_x[i + j] = x_a[i];
            tmpbuf_y[i + j] = y_a[i];
            i++;
        }
        else
        {
            tmpbuf_weights[i + j] = weights_b[j];
            tmpbuf_x[i + j] = x_b[j];
            tmpbuf_y[i + j] = y_b[j];
            j++;
        }
    }

    // Now add the rest from whichever array is not done
    if (j < len_arr_b)
    {
        memcpy(&tmpbuf_weights[i + j], &weights_b[j], sizeof(unsigned int) * (len_arr_b - j));
        memcpy(&tmpbuf_x[i + j], &x_b[j], sizeof(int) * (len_arr_b - j));
        memcpy(&tmpbuf_y[i + j], &y_b[j], sizeof(int) * (len_arr_b - j));
    }
    else if (i < len_arr_a)
    {
        memcpy(&tmpbuf_weights[i + j], &weights_a[i], sizeof(unsigned int) * (len_arr_a - i));
        memcpy(&tmpbuf_x[i + j], &x_a[i], sizeof(int) * (len_arr_a - i));
        memcpy(&tmpbuf_y[i + j], &y_b[j], sizeof(int) * (len_arr_b - j));
    }
}

// The most naive parallel merge sort possible - quite possibly worse than sequential // TODO do better
__global__ void kernel_sort_particles(unsigned int nparticles, int *particles_x, int *particles_y, float *particles_weights,
                                                               int *tmpbuf_x,    int *tmpbuf_y,    float *tmpbuf_weights)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < nparticles)
    {
        // All threads grab their corresponding input element(s)
        float weight = particles_weights[index];

        // Every other thread merges their input with their neighbor
        // Binary reduction merge
        for (unsigned int stride = 1; stride < nparticles; stride *= 2)
        {
            if (index >= stride)
            {
                // The first (stride / 2) elements are sorted and the second (stride / 2) elements are sorted
                // The second half though may be less than stride / 2, if we are at the end of the reduction.
                unsigned int len_arr_a = ceil(stride / 2.0f);
                unsigned int len_arr_b = MIN(ceil(stride / 2.0f), nparticles - len_arr_a);
                unsigned int start_a = index - stride;
                unsigned int start_b = start_a + len_arr_a;

                // Merge
                float *weights_a = &particles_weights[start_a];
                float *weights_b = &particles_weights[start_b];
                int *x_a = &particles_x[start_a];
                int *x_b = &particles_x[start_b];
                int *y_a = &particles_y[start_a];
                int *y_b = &particles_y[start_b];
                // Since each thread is writing to the same global array, we need to make sure they are only
                // writing to their appropriate subsection.
                // The start of each thread's output array should be given by the following formula.
                unsigned int tmpbuf_start = (index - 1) * (2 * stride);
                kernel_sequential_merge(&tmpbuf_weights[tmpbuf_start], weights_a, weights_b, &tmpbuf_x[tmpbuf_start], x_a, x_b, &tmpbuf_y[tmpbuf_start], y_a, y_b, len_arr_a, len_arr_b);
            }
            // Since we are doing a reduction, we need to make sure each thread is done before moving on.
            __syncthreads();
        }
    }
}

int device_resample_and_move(int estimated_vx, int estimated_vy, unsigned int nparticles, int *particles_x, int *particles_y, float *particles_weights, std::mt19937 &rng, unsigned int *indices, int nthreads_per_block)
{
    #if 1
    cudaError_t err;
    int *dev_particles_x = nullptr;
    int *dev_particles_y = nullptr;
    float *dev_weights = nullptr;
    unsigned int *dev_indices = nullptr;
    float *dev_sum_tmp = nullptr;   // The temporary results from each block during sum
    float *dev_sort_weights_tmp = nullptr;
    int *dev_sort_x_tmp = nullptr;
    int *dev_sort_y_tmp = nullptr;
    float *sum_tmp = nullptr;
    float summed_weights = 0.0;
    int nblocks = ceil(nparticles / (float)nthreads_per_block);

    #define CHECK_CUDA_ERR(err) do { if (err != cudaSuccess) { err = (cudaError_t)__LINE__; goto fail; }} while (0)

    /* Allocate everything we need */
    err = cudaMalloc(&dev_particles_x, nparticles * sizeof(int));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_particles_y, nparticles * sizeof(int));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_weights, nparticles * sizeof(float));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_indices, nparticles * sizeof(unsigned int));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_sum_tmp, nblocks * sizeof(float));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_sort_weights_tmp, nparticles * sizeof(float));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_sort_x_tmp, nparticles * sizeof(int));
    CHECK_CUDA_ERR(err);

    err = cudaMalloc(&dev_sort_y_tmp, nparticles * sizeof(int));
    CHECK_CUDA_ERR(err);

    /* Copy everything to the device */
    err = cudaMemcpy(dev_particles_x, particles_x, nparticles * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(dev_particles_y, particles_y, nparticles * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(dev_weights, particles_weights, nparticles * sizeof(float), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(dev_indices, indices, nparticles * sizeof(unsigned int), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERR(err);

    /* Launch kernels */
    kernel_normalize_weights_reduction<<<nblocks, nthreads_per_block, (sizeof(float) * nthreads_per_block)>>>(nparticles, dev_weights, dev_sum_tmp);
    err = cudaDeviceSynchronize();
    CHECK_CUDA_ERR(err);

    // Sequential sum of the intermediate results in dev_sum_tmp
    sum_tmp = (float *)malloc(nblocks * sizeof(float));
    err = cudaMemcpy(sum_tmp, dev_sum_tmp, nblocks * sizeof(float), cudaMemcpyDeviceToHost);
    CHECK_CUDA_ERR(err);
    for (unsigned int i = 0; i < nblocks; i++)
    {
        summed_weights += sum_tmp[i];
    }
    free(sum_tmp);
    sum_tmp = nullptr;

    kernel_normalize_weights_complete<<<nblocks, nthreads_per_block>>>(nparticles, dev_weights, summed_weights);
    err = cudaDeviceSynchronize();
    CHECK_CUDA_ERR(err);

    kernel_sort_particles<<<nblocks, nthreads_per_block>>>(nparticles, dev_particles_x, dev_particles_y, dev_weights, dev_sort_x_tmp, dev_sort_y_tmp, dev_sort_weights_tmp);
    err = cudaDeviceSynchronize();
    CHECK_CUDA_ERR(err);

    free(dev_sort_y_tmp);
    free(dev_sort_x_tmp);
    free(dev_sort_weights_tmp);
    dev_sort_y_tmp = nullptr;
    dev_sort_x_tmp = nullptr;
    dev_sort_weights_tmp = nullptr;

    //kernel_resample_particles
    //kernel_reset_all_weights
    //kernel_move_particles

    /* Transfer results back to host */
    err = cudaMemcpy(particles_x, dev_particles_x, nparticles * sizeof(int), cudaMemcpyDeviceToHost);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(particles_y, dev_particles_y, nparticles * sizeof(int), cudaMemcpyDeviceToHost);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(particles_weights, dev_weights, nparticles * sizeof(float), cudaMemcpyDeviceToHost);
    CHECK_CUDA_ERR(err);

    err = cudaMemcpy(indices, dev_indices, nparticles * sizeof(unsigned int), cudaMemcpyDeviceToHost);
    CHECK_CUDA_ERR(err);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
    //Remove the logic here as you convert it to CUDA
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
    complete_resample_and_move_step(nparticles, particles_weights, rng, indices, particles_x, particles_y, estimated_vx, estimated_vy);
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
    // End
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//

    /* Free up memory */
    err = cudaFree(dev_particles_x);
    dev_particles_x = nullptr;
    CHECK_CUDA_ERR(err);

    err = cudaFree(dev_particles_y);
    dev_particles_y = nullptr;
    CHECK_CUDA_ERR(err);

    err = cudaFree(dev_weights);
    dev_weights = nullptr;
    CHECK_CUDA_ERR(err);

    err = cudaFree(dev_indices);
    dev_indices = nullptr;
    CHECK_CUDA_ERR(err);

    err = cudaFree(dev_sum_tmp);
    dev_sum_tmp = nullptr;
    CHECK_CUDA_ERR(err);

    #undef CHECK_CUDA_ERR

fail:
    if (err != cudaSuccess)
    {
        std::cout << "Error at line " << err << std::endl;
        assert(false);
    }
    return err;
#else
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
#endif
}
