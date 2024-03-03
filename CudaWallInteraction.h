#ifndef __CUDA_WALL_INTERACTION_H__
#define __CUDA_WALL_INTERACTION_H__

#include "my_helper_math.h"
#include "CudaDefines.cuh"
#include <vector>
#include <cstring>
#include "Wall.h"
#include "CompressionInterpolation.h"

#include "myCubWrapper.h"
#include "Simulation.h"

class particleWallInteraction
{
public:

    particleWallInteraction();
    ~particleWallInteraction(void);

    void initParticleWallInteraction(const int number_of_particles);
    void deleteParticleWallInteraction();

    double getGPUTopWallForce(CubPlan &cubPlan);
    double getGPUBotWallForce(CubPlan &cubPlan);

    void get_box_pos_from_gpu(Simulation &sim);
    void load_box_from_gpu(Simulation &sim);
    void load_box_to_gpu(const Simulation &sim);

    void initMaterialConstants(const NormalInteraction &normal_interaction);

    void storeToFile(FILE* file);
    void loadFromFile(FILE* file);

    void updateWallContacts(
		double * RESTRICT gpu_forces_new,
		double * RESTRICT gpu_torques_new,
        const double* RESTRICT  velocities,
		const double* RESTRICT gpu_angular_velocities,
		const double* RESTRICT gpu_torques,
		const double timestep
	);


    void updateWallContactNormals(
            const double* RESTRICT positions
            );

    void updateWallSticking(
            const double* RESTRICT positions,
			const double* RESTRICT gpu_angular_velocities,
			const double* RESTRICT gpu_torques,
            const double timestep
            );

    void updateWallStickingNoSw(
            const double* RESTRICT positions,
            const double* RESTRICT gpu_angular_velocities,
            const double* RESTRICT gpu_torques,
            const double timestep
            );


    bool is_initialized(void){ return this->m_wall_interaction_initialized;}

    double get_height(void){ return this->m_cpu_wall_positions[GPU_WALL_TOP].y - this->m_cpu_wall_positions[GPU_WALL_BOTTOM].y  - 2.0 * particle_radius;}



private:

    int m_number_of_particles;

    bool m_wall_interaction_allocated;
    bool m_wall_gpu_constants_initialized;
    bool m_wall_interaction_initialized;
    bool m_wall_box_loaded;

    // for every particle
    int* m_gpu_wall_particle_ids;		// stores the ids of the walls a certain particle is in contact with

    double* m_gpu_top_wall_force;				// storage for the force the monomers exert on the top wall on the compaction box
    double* m_gpu_bot_wall_force;				// sotrage for the force the monimers exert on the bottom wall of the compation box



    // for every contact
    double*  m_gpu_dist_to_wall;
    double3* m_gpu_wall_contact_vector_new;     // stores contact vector
    double3* m_gpu_wall_contact_vector_old;
    double3* m_gpu_wall_contact_n_initial;		// initial contact pointer
    double4* m_gpu_wall_contact_rot;			// rotation state
    double3* m_gpu_wall_contact_n;				// current normal
	double3* m_gpu_wall_initial_contact;		// position of the contact relative to the edge of the wall

    double*  m_gpu_wall_twisting_displacement;


    double3* m_cpu_wall_positions;              // cpu positions of the walls
    double3* m_gpu_wall_positions;              // gpu positions of the walls
public:
    GPUWallConstants m_wallConstants;

};




#endif // __CUDA_WALL_INTERACTION_H__
