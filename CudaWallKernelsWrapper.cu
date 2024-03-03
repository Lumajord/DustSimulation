#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <helper_cuda.h>

#include "cuda.h"
#include <helper_functions.h>
#include "CudaDefines.cuh"
//#include "CudaWallKernelsWrapper.cuh"
#include "CudaWallKernels.cuh"


//Round a / b to nearest higher integer value
static unsigned int iDivUp(const unsigned int a, const unsigned int b)
{
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

static void computeGridSize(
	const unsigned int n,
	const unsigned int blockSize,
	unsigned int &numBlocks,
	unsigned int &numThreads
)
{
	numThreads = min(blockSize, n);
	numBlocks = iDivUp(n, numThreads);
}



extern "C"
{

    void allocateArray(void **devPtr, size_t size)
    {
        checkCudaErrors(cudaMalloc(devPtr, size));
    }

    void callocateArray(void **devPtr, size_t size)
    {
        checkCudaErrors(cudaMalloc(devPtr, size));
        checkCudaErrors(cudaMemset(*devPtr, 0, size));
    }

    void freeArray(void *devPtr)
    {
        checkCudaErrors(cudaFree(devPtr));
    }


    void memsetArray(void *devPtr, int value, size_t size)
    {
        checkCudaErrors(cudaMemset(devPtr, value, size));
    }

    void copyToGPU(void *devPtr, void *hostPtr, size_t size)
    {
        checkCudaErrors(cudaMemcpy(devPtr, hostPtr, size, cudaMemcpyHostToDevice));
    }

    void copyFromGPU(void *hostPtr, void *devPtr, size_t size)
    {
        checkCudaErrors(cudaMemcpy(hostPtr, devPtr, size, cudaMemcpyDeviceToHost));
    }

    void init_GPUWallConstansWrapper(const GPUWallConstants &cpu_constants)
    {
        init_GPUWallConstansKernel(cpu_constants);
    }


    void updateWallContactNormalsWrapper(
        int* RESTRICT contact_wall_particle_ids,
        double * RESTRICT dist_to_wall,
        double3* RESTRICT contact_normals_new,
        const double* RESTRICT positions,
        const double3* RESTRICT wall_initial_contact,
        const double3* RESTRICT wall_positions,
        const int number_of_particles
            )
    {

        unsigned int numThreads, numBlocks;
        computeGridSize(MAX_WALL_CONTACTS*number_of_particles, BLOCK_SIZE, numBlocks, numThreads);



        updateWallContactNormalsKernel <<< numBlocks, numThreads >>>(
            contact_wall_particle_ids,
            dist_to_wall,
            contact_normals_new,
            positions,
            wall_initial_contact,
            wall_positions,
            number_of_particles
                );

        cudaDeviceSynchronize();
        getLastCudaError("updateWallContactNormalsKernel failed");

    }




    void updateWallStickingWrapper(
        int* RESTRICT contact_wall_particle_ids,
        const double* RESTRICT positions,
        double * RESTRICT dist_to_wall,
        double3* RESTRICT wall_contact_vector_new,
        double3* RESTRICT wall_contact_vector_old,
        double3* RESTRICT wall_contact_n_intitial,
        double4* RESTRICT wall_contact_rot,
        double3* RESTRICT wall_contact_n,
        double3* RESTRICT wall_initial_contact,
        double * RESTRICT wall_twisting_displacement,
        const double3* RESTRICT wall_pos,
        const int number_of_particles,

		const double* RESTRICT  angular_velocities,
		const double* RESTRICT  torques,
		const double			timestep
            )
    {

        unsigned int numThreads, numBlocks;
        computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);


        updateRotation << < numBlocks, numThreads >> > (
			wall_contact_rot,
			angular_velocities,
			torques,
			contact_wall_particle_ids,
			timestep,
			number_of_particles
			);


        updateWallStickingKernel<6><<< numBlocks, numThreads>>>(
            contact_wall_particle_ids,
            positions,
            dist_to_wall,
            wall_contact_vector_new,
            wall_contact_vector_old,
            wall_contact_n_intitial,
            wall_contact_rot,
            wall_contact_n,
            wall_initial_contact,
            wall_twisting_displacement,
            wall_pos,
            number_of_particles
                );

        cudaDeviceSynchronize();
        getLastCudaError("updateWallStickingKernel failed");

    }



    void updateWallStickingWrapperNoSw(
        int* RESTRICT contact_wall_particle_ids,
        const double* RESTRICT positions,
        double * RESTRICT dist_to_wall,
        double3* RESTRICT wall_contact_vector_new,
        double3* RESTRICT wall_contact_vector_old,
        double3* RESTRICT wall_contact_n_intitial,
        double4* RESTRICT wall_contact_rot,
        double3* RESTRICT wall_contact_n,
        double3* RESTRICT wall_initial_contact,
        double * RESTRICT wall_twisting_displacement,
        const double3* RESTRICT wall_pos,
        const int number_of_particles,

        const double* RESTRICT  angular_velocities,
        const double* RESTRICT  torques,
        const double			timestep
            )
    {

        unsigned int numThreads, numBlocks;
        computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);


        updateRotation << < numBlocks, numThreads >> > (
            wall_contact_rot,
            angular_velocities,
            torques,
            contact_wall_particle_ids,
            timestep,
            number_of_particles
            );


        updateWallStickingKernel<2><<< numBlocks, numThreads>>>(
            contact_wall_particle_ids,
            positions,
            dist_to_wall,
            wall_contact_vector_new,
            wall_contact_vector_old,
            wall_contact_n_intitial,
            wall_contact_rot,
            wall_contact_n,
            wall_initial_contact,
            wall_twisting_displacement,
            wall_pos,
            number_of_particles
                );

        cudaDeviceSynchronize();
        getLastCudaError("updateWallStickingKernelNoSw failed");

    }

    void updateWallContactsWrapper(
        double4* RESTRICT       rot,
        double3* RESTRICT       n1_initial,
        double3* RESTRICT       n_1,
        double * RESTRICT       forces_new,
        double * RESTRICT       torques_new,
        double * RESTRICT       wall_top_force, // force the particle exert on the top wall
        double * RESTRICT       wall_bot_force, // force the particle exert on the bottom wall
        const double * RESTRICT dist_to_wall,
        const double3* RESTRICT contact_normals_new,
        const double3* RESTRICT contact_normals_old,
              double3* RESTRICT initial_contact, // = n2_initial in alex code
        const double* RESTRICT  velocities,
        const double* RESTRICT  angular_velocities,
        const double* RESTRICT  torques,
        const int* RESTRICT     contact_wall_particle_ids,
        double* RESTRICT        twisting_displacement_arr,
        const double            timestep,
        const uint              number_of_particles
        )
    {

        unsigned int numThreads, numBlocks;
        computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);


        updateWallContactsKernel<<< numBlocks, numThreads >>>(
            rot,
            n1_initial,
            n_1,
            forces_new,
            torques_new,
            wall_top_force,
            wall_bot_force,
            dist_to_wall,
            contact_normals_new,
            contact_normals_old,
            initial_contact,
            velocities,
            angular_velocities,
            torques,
            contact_wall_particle_ids,
            twisting_displacement_arr,
            timestep,
            number_of_particles
            );

        cudaDeviceSynchronize();
        getLastCudaError("updateWallContactsKernel failed");

    }



} // extern "C"
