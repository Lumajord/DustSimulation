#ifndef __CUDA_WALL_KERNEL_WRAPPER_H__
#define __CUDA_WALL_KERNEL_WRAPPER_H__


#include "CudaDefines.cuh"
#include "my_helper_math.h"

extern "C"
{

    void allocateArray(void **devPtr, int size_t);
    void callocateArray(void **devPtr, int size_t);
    void freeArray(void *devPtr);
    void memsetArray(void *devPtr, int value, int size_t);
    void copyToGPU(void *devPtr, void *hostPtr, size_t size);
    void copyFromGPU(void *hostPtr, void *devPtr, size_t size);

    void init_GPUWallConstansWrapper(const GPUWallConstants &cpu_constants);

    void updateWallContactNormalsWrapper(
        int* RESTRICT contact_wall_particle_ids,
        double * RESTRICT dist_to_wall,
        double3* RESTRICT contact_normals_new,
        const double* RESTRICT positions,
        const double3* RESTRICT wall_initial_contact,
        const double3* RESTRICT wall_positions,
        const int number_of_particles
            );


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
            );

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
            );

    void updateWallContactsWrapper(
        double4* RESTRICT       rot,
        double3* RESTRICT       n1_initial,
        double3* RESTRICT       n_1,
        double * RESTRICT       forces_new,
        double * RESTRICT       torques_new,
        double * RESTRICT       wall_top_force, // force the particle exert on the top wall
        double * RESTRICT       wall_bot_force,
        const double * RESTRICT dist_to_wall,
        const double3* RESTRICT contact_normals_new,
        const double3* RESTRICT contact_normals_old,
              double3* RESTRICT initial_contact, // = n2_initial in alex code
        const double* RESTRICT  velocities,
        const double* RESTRICT  angular_velocities,
        const double* RESTRICT  torques,
        const int* RESTRICT     contact_wall_particle_ids,
        const double* RESTRICT  twisting_displacement_arr,
        const double            timestep,
        const uint              number_of_particles
        );


} // extern "C"




#endif // __CUDA_WALL_KERNEL_WRAPPER_H__
