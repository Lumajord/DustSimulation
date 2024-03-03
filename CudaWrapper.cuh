#ifndef CUDA_WRAPPER_H
#define CUDA_WRAPPER_H

#include "CudaDefines.cuh"

class NormalInteraction;

extern "C"
{


void createUnsignedIntegerTextureObject(
        unsigned int* devPtr,
        cudaTextureObject_t &texObj,
        unsigned int sizeInBytes
        );


    void initGPUConstants(
            const NormalInteraction normalInteraction,
            const unsigned int gpu_num_grid_cells
            );

    void getGPUConstants(GPUConstants &d_gpu_constants);
    void setGPUConstants(const GPUConstants &d_gpu_constants);


    void cudaPosPredictor(
            double* new_pos,
            const double *old_pos,
            const double *velocity,
            const double *force,
            const double timestep,
            unsigned int *grid_particle_cell,
            unsigned int *grid_particle_index,
            const unsigned int number_of_particles
            );

    void cudaCorrector(
            double* velocity,
            double* angular_velocity,
            double* old_force,
            double* new_force,
            double* old_torque,
            double* new_torque,
            const double timestep,
            const unsigned int number_of_particles
            );

    void cudaDampingCorrector(
            double* velocity,
            double* angular_velocity,
            double* old_force,
            double *new_force,
            double* old_torque,
            double* new_torque,
            const double damping_factor,
            const double timestep,
            const unsigned int number_of_particles
            );


    void cudaFindCellStart(
            unsigned int *cell_start,
            unsigned int *cell_end,
            unsigned int *grid_particle_cell,
            unsigned int *grid_particle_index,
            double* gpu_positions_new,
            double* gpu_positions_sorted,
            const unsigned int number_of_particles,
			const unsigned int num_grid_cells
            );

/*
    void updateContactNormalsKernel(
        #ifdef GPU_TRACK_NUMBER_OF_CONTACTS
            int *contacts_broken,
        #endif
            int *next_free_contact_id,
            int *free_contacts_list,
            int2 *contact_particle_ids,
            double4 *contact_normals_new,
            double *positions,
            unsigned int number_of_particles
            );
*/


    void cudaUpdateSticking(
            unsigned int &total_contacts_created,
            unsigned int &total_contacts_broken,
        #ifdef GPU_TRACK_DISSIPATED_ENERGY
            double  *dissipated_contact_energy,
            double4 *n_1,
            double4 *n_2,
        #endif // GPU_TRACK_DISSIPATED_ENERGY
            int* update_local_contact_list,
            int *next_free_contact_id,
            int *free_contacts_list,
            int2 *new_contact_particle_ids,
            int *number_of_new_contacts,
            int* particle_number_of_contacts,
            int* particle_particle_ids,
            int* particle_contact_ids,
            int2* contact_particle_ids,
            const unsigned int *cell_start,
            const unsigned int *cell_end,
            const unsigned int *grid_particle_index,
            double4 *contact_normals_new,
            double *positions,
            double *positions_sorted,
            const int number_of_particles
            );


    void cudaAddNewContacts(int *particle_number_of_contacts,
            int *particle_particle_ids,
            int *particle_contact_ids,
            int *next_free_contact_id,
            int *free_contacts_list,
            int2 *new_contact_particle_ids,
            const int number_of_new_contacts,
            int2 *contact_particle_ids,
            double4* contact_normals_new,
            double4* contact_normals_old,
            double4 *rot1,
            double4 *rot2,
            double3 *n1_initial,
            double3 *n2_initial,
            double4 *n1,
            double *positions,
            const int number_of_particles);

    void cudaUpdateContacts(
            int2 *contact_particle_ids,
            double4 *rot1,
            double4 *rot2,
            double3 *n1_initial,
            double3 *n2_initial,
            double4 *n_1,
            double4 *n_2,
            double4 *contact_normals_new,
            double4 *contact_normals_old,
            double *positions,
            double *angular_velocities,
            double *torques,
            const double timestep,
            const unsigned int number_of_particles
        #ifdef GPU_TRACK_DISSIPATED_ENERGY
            ,
            double* dissipated_rolling_energy,
            double* dissipated_sliding_energy,
            double* dissipated_twisting_energy
        #endif // GPU_TRACK_DISSIPATED_ENERGY
            );

    void cudaUpdateInteraction(
            double *forces_new,
            double *torques_new,
            int *particle_number_of_contacts,
            int *contact_ids,
            int2 *contact_particle_ids,
            double4 *n_1,
            double4 *n_2,
            double4 *contact_normals_new,
            const unsigned int number_of_particles
        #ifdef GPU_TRACK_DISSIPATED_ENERGY
            ,
            double* dissipated_damping_energy
        #endif // GPU_TRACK_DISSIPATED_ENERGY
            );

    void cudaUpdateBoxInteraction(
            double *forces_new,
            double *positions,
            const double3 lower_pos,
            const double3 upper_pos,
			double* top_wall_force,
			double* bot_wall_force,
            const unsigned int number_of_particles
            );

    void cudaUpdateOpenBoxInteraction(
            double *forces_new,
            double *positions,
            const double3 lower_pos,
            const double3 upper_pos,
            const unsigned int number_of_particles
            );

    void cudaUpdateEnclosingSphereWallInteractionKernel(
            double *forces_new,
            double *positions,
            const double sphere_radius,
            const unsigned int number_of_particles
            );

    void cudaUpdateBottomWallInteraction(
            double *forces_new,
            double *positions,
            const double3 lower_pos,
            const unsigned int number_of_particles
            );

#ifdef GPU_TRACK_PARTICLE_ORIENTATION
    void cudaUpdateParticleOrientation(
            double4 *particle_orientations,
            double *angular_velocities,
            double *torques,
            const double timestep,
            const unsigned int number_of_particles
            );
#endif // GPU_TRACK_PARTICLE_ORIENTATION

#ifdef GPU_TRACK_DISSIPATED_ENERGY
    class CubPlan;
    void cudaGetEnergy(CubPlan &sortPlan,
            double *gpu_velocities,
            double *gpu_dissipated_damping_energy,
            double *gpu_dissipated_rolling_energy,
            double *gpu_dissipated_contact_energy,
            double *gpu_dissipated_sliding_energy,
            double *gpu_dissipated_twisting_energy,
            double *kinetic_energies,
            double *rotation_energies,
            double *normal_energies,
            double *rolling_energies,
            double *sliding_energies,
            double *twisting_energies,
            double *velocities,
            double *angular_velocities,
            int2 *contact_particle_ids,
            double4 *contact_normals,
            double4 *n_1,
            double4 *n_2,
            double &E_tot,
            double &E_kin,
            double &E_rot,
            double &V_tot,
            double &V_normal,
            double &V_rolling,
            double &V_sliding,
            double &V_twisting,
            double &diss_normal,
            double &diss_rolling,
            double &diss_sliding,
            double &diss_twisting,
            double &diss_contact,
            double3 &v1,
            double3 &v2,
            const double timestep,
            const unsigned int number_of_particles
            );
#endif // GPU_TRACK_DISSIPATED_ENERGY
} // extern "C"

#endif // CUDA_WRAPPER_H
