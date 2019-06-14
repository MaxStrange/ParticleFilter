#include <algorithm>
#include <assert.h>
#include <cuda_runtime.h>
#include <random>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "cudaparticlefilter.h"
#include "device.h"
#include "robot.h"

#define DEBUG 1
#define DEBUG_PRINTF(fmt, ...) do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)

CudaParticleFilter::CudaParticleFilter(unsigned int nparticles, unsigned int screen_height, unsigned int screen_width)
    : ParticleFilter(nparticles)
{
    std::random_device rnddev;
    this->rng = std::mt19937(rnddev());
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

CudaParticleFilter::~CudaParticleFilter(void)
{
    // Parent's destructor called automatically

    // Free the indices we allocated in our own constructor
    if (this->indices != nullptr)
    {
        free(this->indices);
        this->indices = nullptr;
    }
}

unsigned int CudaParticleFilter::get_nparticles(void) const
{
    return this->nparticles;
}

int CudaParticleFilter::get_xpos(unsigned int index) const
{
    return this->particles_x[index];
}

int CudaParticleFilter::get_ypos(unsigned int index) const
{
    return this->particles_y[index];
}

int CudaParticleFilter::init(void)
{
    cudaError_t err;
    cudaDeviceProp dev_prop;
    int dev_id;

    // Check if there is a CUDA-enabled device and panic if not
    err = cudaGetDevice(&dev_id);
    if (err != cudaSuccess)
    {
        goto fail;
    }

    err = cudaGetDeviceProperties(&dev_prop, dev_id);
    if (err != cudaSuccess)
    {
        goto fail;
    }

    this->print_dev_props(dev_id, dev_prop);

fail:
    return err;
}

void CudaParticleFilter::update_part1(Robot &robot)
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
    device_calculate_likelihood(this->particles_x, this->particles_y, estimate_x, estimate_y, this->particles_weights, this->nparticles);

    // Take a break now to update the graphics in the main loop
}

void CudaParticleFilter::update_part2(Robot &robot)
{
    int estimate_vx;
    int estimate_vy;
    robot.get_v_estimate(&estimate_vx, &estimate_vy);
    device_resample_and_move(estimate_vx, estimate_vy, this->nparticles, this->particles_x, this->particles_y, this->particles_weights, this->rng, this->indices);
}

void CudaParticleFilter::print_dev_props(int dev_id, const cudaDeviceProp &dev_prop) const
{
    printf("GPU device %d (%s):\n", dev_id, dev_prop.name);
    printf("    Total global memory: %zu GB\n", dev_prop.totalGlobalMem / 1000000000);
    printf("    Shared memory per block: %zu KB\n", dev_prop.sharedMemPerBlock / 1000);
    printf("    Registers per block: %d registers\n", dev_prop.regsPerBlock);
    printf("    Warp size: %d threads\n", dev_prop.warpSize);
    printf("    Max threads per block: %d threads\n", dev_prop.maxThreadsPerBlock);
    printf("    Max threads per dimension: %d (x) %d (y) and %d (z)\n", dev_prop.maxThreadsDim[0], dev_prop.maxThreadsDim[1], dev_prop.maxThreadsDim[2]);
    printf("    Max blocks per grid dimension: %d %d and %d\n", dev_prop.maxGridSize[0], dev_prop.maxGridSize[1], dev_prop.maxGridSize[2]);
    printf("    Clock rate: %d MHz\n", dev_prop.clockRate / 1000);
    printf("    Total constant memory: %zu bytes\n", dev_prop.totalConstMem);
    printf("    Major and Minor version: %d.%d\n", dev_prop.major, dev_prop.minor);
    printf("    Multiprocessor count: %d\n", dev_prop.multiProcessorCount);
    printf("    Concurrent Kernels: %d\n", dev_prop.concurrentKernels);
    printf("    Unified addressing: %s\n", dev_prop.unifiedAddressing ? "Yes" : "No");
    printf("    Memory clock rate: %d MHz\n", dev_prop.memoryClockRate / 1000);
    printf("    Memory bus width: %d bits\n", dev_prop.memoryBusWidth);
    printf("    L2 Cache Size: %d KB\n", dev_prop.l2CacheSize / 1000);
    printf("    Threads per multiprocessor: %d threads\n", dev_prop.maxThreadsPerMultiProcessor);
    printf("    Global L1 Cache supported: %s\n", dev_prop.globalL1CacheSupported ? "Yes" : "No");
    printf("    Local L1 Cache supported: %s\n", dev_prop.localL1CacheSupported ? "Yes" : "No");
    printf("    Registers per multiprocessor: %d registers\n", dev_prop.regsPerMultiprocessor);
    printf("    Single to float precision performance ratio: %d\n", dev_prop.singleToDoublePrecisionPerfRatio);
    printf("    Pageable memory access: %s\n", dev_prop.pageableMemoryAccess ? "Yes" : "No");
    printf("    Concurrent managed access: %s\n", dev_prop.concurrentManagedAccess ? "Yes" : "No");
}
