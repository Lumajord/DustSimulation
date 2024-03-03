#include <math_constants.h>
#include "helper_cuda.h"

#include "Constants.h"
#include "CudaDefines.cuh"
#include "CudaWrapper.cuh"
#include "CudaKernels.cuh"

#include <iostream>
using namespace std;

extern double ENERGY_UNIT;


extern "C"
void createUnsignedIntegerTextureObject(
        unsigned int* devPtr,
        cudaTextureObject_t &texObj,
        unsigned int sizeInBytes
        )
{
    // create texture objects for unsigned integer type data
    cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = devPtr;
    resDesc.res.linear.desc.f = cudaChannelFormatKindUnsigned;
    resDesc.res.linear.desc.x = 32; // bits per channel
    resDesc.res.linear.sizeInBytes = sizeInBytes;

    cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.readMode = cudaReadModeElementType;

    // create texture object: we only have to do this once!
    cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL);

    cudaDeviceSynchronize();
    getLastCudaError("Creating Texture Failed");

}




static double getContactRadiusBuild(const double compression_length, double particle_radius, double c1_contact_radius, double c2_contact_radius)
{
	// contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
	// use equilibrium contact radius as starting value
	const double k = particle_radius * compression_length / 6.0;
	double x_pow3;
	double x_new;
	double x_old = c1_contact_radius;

	// use Newton-Raphson method to find root
    for (int i = 0; i < 200; ++i)
	{
		x_pow3 = x_old * x_old * x_old;
		x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius);

        if (fabs(x_new - x_old) / particle_radius < 1.e-14)
			return x_new * x_new;

		x_old = x_new;
	}

	return x_new * x_new;
}


#include "CompressionInterpolation.h"


extern double particle_radius;		// in cm
extern double particle_radius_inv;	// in cm^-1
extern double reduced_radius;		// in cm

extern double young_mod_reduced;


// material parameters
extern double mass;				// in g
extern double mass_inv;			// in g^-1
extern double moment_of_inertia;	// in g cm^2
extern double moment_of_inertia_inv;
extern double osc_damping_factor;	// additional damping for normal oscillations
extern double viscoelastic_damping_constant;  // viscoelastic damping constant as proposed by Sebastian Krijt

extern double equilibrium_contact_radius;	// radius of contact area when no normal is acting on the particles
extern double delta_0;				// in cm
extern double delta_c;				// in cm
extern double F_c;					// in 10^-5 N = 10^-2 mJ/m

extern double k_s;	// constant for sliding potential
extern double k_r;	// constant for rolling potential
extern double k_t;	// constant for twisting potential

extern double contact_breaking_dist;			// in cm
extern double contact_breaking_dist_squared;	// in cm^2
extern double contact_making_dist;				// in cm
extern double contact_making_dist_squared;		// in cm^2

extern double crit_rolling_displacement;			// in cm
extern double crit_rolling_displacement_squared;	// in cm^2
extern double crit_sliding_displacement;
extern double crit_sliding_displacement_squared;
extern double crit_twisting_displacement;
extern double crit_twisting_displacement_squared;

extern "C"
void initGPUConstants(
        const NormalInteraction normalInteraction,
        const unsigned int num_grid_cells
        )
{


	// grid specifications
	
    double gpu_grid_cell_width = (2.001 * particle_radius); // CHANGE: 2.00001 -> 2.0
	double gpu_grid_cell_width_inv = 1.0 / gpu_grid_cell_width;
	double gpu_grid_shift = (0.5 * gpu_grid_cell_width * (double)num_grid_cells);


    GPUConstants temp;
    temp.eq_contact_radius = equilibrium_contact_radius;


    temp.c1_contact_radius = normalInteraction.c1_contact_radius;
    temp.c2_contact_radius = normalInteraction.c2_contact_radius;

    temp.c1_normal_force = normalInteraction.k1;
    temp.c2_normal_force = -normalInteraction.k2;
    temp.delta_c = delta_c;
    temp.delta_0 = delta_0;

    temp.contact_breaking_dist = contact_breaking_dist;
    temp.contact_breaking_dist_squared = contact_breaking_dist_squared;

    temp.contact_making_dist_squared = contact_making_dist_squared;

    temp.particle_radius = particle_radius;
#ifdef CPU_EQUIVALENT
    temp.reduced_radius = reduced_radius;
#endif
    temp.particle_radius_inv = particle_radius_inv;
    temp.mass = mass;
    temp.mass_inv = mass_inv;
    temp.moment_of_inertia = moment_of_inertia;
    temp.moment_of_inertia_inv = moment_of_inertia_inv;

    temp.F_c = F_c;

    temp.k_r = k_r * reduced_radius * reduced_radius;
    temp.k_s = k_s * particle_radius * particle_radius;
    temp.k_t = k_t;

    temp.crit_rolling_displacement = crit_rolling_displacement;


    temp.crit_rolling_displacement_squared = crit_rolling_displacement_squared;
    temp.crit_sliding_displacement = crit_sliding_displacement;
    temp.crit_sliding_displacement_squared = crit_sliding_displacement_squared;
    temp.crit_twisting_displacement = crit_twisting_displacement;

#ifdef GPU_USE_CONSTANT_DAMPING
    temp.osc_damping_factor = osc_damping_factor;
#endif

#ifdef GPU_USE_VISCOELASTIC_DAMPING
    temp.osc_damping_factor = viscoelastic_damping_constant;
#endif

    temp.gpu_grid_cell_width_inv = gpu_grid_cell_width_inv;
    temp.gpu_num_grid_cells = num_grid_cells;
    temp.gpu_grid_shift = gpu_grid_shift;

	

	double get_contact_radius_c1 = getContactRadiusBuild(-temp.delta_c, temp.particle_radius, temp.c1_contact_radius, temp.c2_contact_radius);
	double get_contact_radius_c2 = (temp.eq_contact_radius - get_contact_radius_c1) / sqrt(temp.delta_c + temp.delta_0);

	temp.get_contact_radius_c1 = get_contact_radius_c1;
	temp.get_contact_radius_c2 = get_contact_radius_c2;

    temp.k_hertz = 4.0 / 3.0 * young_mod_reduced * sqrt(reduced_radius);


    /*
	cout << endl;
	cout << "eq_contact_radius = " << temp.eq_contact_radius << endl;
	cout << "c1_contact_radius = " << temp.c1_contact_radius << endl;
	cout << "c2_contact_radius = " << temp.c2_contact_radius << endl;

	cout << "c1_normal_force = " << temp.c1_normal_force << endl;
	cout << "c2_normal_force = " << temp.c2_normal_force << endl;
	cout << "delta_c = " << temp.delta_c << endl;
	cout << "delta_0 = " << temp.delta_0 << endl;

	cout << "contact_breaking_dist = " << temp.contact_breaking_dist << endl;
	cout << "contact_breaking_dist_squared = " << temp.contact_breaking_dist_squared << endl;

	cout << "particle_radius_inv = " << temp.particle_radius_inv << endl;
	cout << "mass = " << temp.mass << endl;
	cout << "mass_inv = " << temp.mass_inv << endl;
	cout << "moment_of_inertia = " << temp.moment_of_inertia << endl;

	cout << "k_r = " << temp.k_r << endl;
	cout << "k_s = " << temp.k_s << endl;
	cout << "k_t = " << temp.k_t << endl;

	cout << "crit_rolling_displacement = " << temp.crit_rolling_displacement << endl;
	cout << "crit_rolling_displacement_squared = " << temp.crit_rolling_displacement_squared << endl;
	cout << "crit_sliding_displacement_squared = " << temp.crit_sliding_displacement_squared << endl;
	cout << "crit_twisting_displacement = " << temp.crit_twisting_displacement << endl;

	cout << "osc_damping_factor = " << temp.osc_damping_factor << endl;
	cout << endl << endl;
    */
	checkCudaErrors(cudaMemcpyToSymbol(gpu_constants, &temp, sizeof(GPUConstants)));
}



extern "C"
void getGPUConstants(GPUConstants &d_gpu_constants)
{
    checkCudaErrors(cudaMemcpyFromSymbol(&d_gpu_constants, gpu_constants, sizeof(GPUConstants)));
}


extern "C"
void setGPUConstants(const GPUConstants &d_gpu_constants)
{
    checkCudaErrors(cudaMemcpyToSymbol(gpu_constants, &d_gpu_constants, sizeof(GPUConstants)));
}



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
void cudaPosPredictor(
        double* new_pos,
        const double* old_pos,
        const double* velocity,
        const double* force,
        const double timestep,
        unsigned int *grid_particle_cell,
        unsigned int *grid_particle_index,
        const unsigned int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    posPredictionKernel<<< numBlocks, numThreads >>>(
                                                       new_pos,
                                                       old_pos,
                                                       velocity,
                                                       force,
                                                       timestep,
                                                       grid_particle_cell,
                                                       grid_particle_index,
                                                       number_of_particles
                                                       );
    cudaDeviceSynchronize();
    getLastCudaError("posPredictionKernel failed");


}

extern "C"
void cudaCorrector(
        double* velocity,
        double* angular_velocity,
        double* old_force,
        double *new_force,
        double* old_torque,
        double* new_torque,
        const double timestep,
        const unsigned int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    correctionKernel<<< numBlocks, numThreads >>>(
                                                    velocity,
                                                    angular_velocity,
                                                    old_force,
                                                    new_force,
                                                    old_torque,
                                                    new_torque,
                                                    timestep,
                                                    number_of_particles
                                                    );
    cudaDeviceSynchronize();
	getLastCudaError("correctionKernel failed");


}

extern "C"
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
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    dampingCorrectionKernel<<< numBlocks, numThreads >>>(
                                                           velocity,
                                                           angular_velocity,
                                                           old_force,
                                                           new_force,
                                                           old_torque,
                                                           new_torque,
                                                           damping_factor,
                                                           timestep,
                                                           number_of_particles
                                                           );
    cudaDeviceSynchronize();
	getLastCudaError("dampingKernel failed");
}


extern "C"
void cudaFindCellStart(
        unsigned int *cell_start,
        unsigned int *cell_end,
        unsigned int *grid_particle_cell,
        unsigned int *grid_particle_index,
        double *gpu_positions_new,
        double *gpu_positions_sorted,
        const unsigned int number_of_particles,
		const unsigned int num_grid_cells
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);
	
	// set all cells to be empty
	checkCudaErrors(cudaMemset(cell_start, NO_PARTICLE, num_grid_cells * num_grid_cells * num_grid_cells * sizeof(unsigned int)));

	unsigned int shared_mem_size = sizeof(unsigned int)*(numThreads+1);

    findCellStartKernel<<< numBlocks, numThreads, shared_mem_size >>>(
                                                                        cell_start,
                                                                        cell_end,
                                                                        grid_particle_cell,
                                                                        grid_particle_index,
                                                                        gpu_positions_new,
                                                                        gpu_positions_sorted,
                                                                        number_of_particles
                                                                        );
    cudaDeviceSynchronize();
	getLastCudaError("findCellStartKernel failed");

	
}

extern "C"
void cudaUpdateSticking(
        unsigned int &total_contacts_created,
        unsigned int &total_contacts_broken,
    #ifdef GPU_TRACK_DISSIPATED_ENERGY
        double *dissipated_contact_energy,
        double4* n_1,
        double4* n_2,
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
        const unsigned int* cell_start,
        const unsigned int* cell_end,
        const unsigned int* grid_particle_index,
        double4 *contact_normals_new,
        double *positions,
        double *positions_sorted,
        const int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(6 * number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate new contact pointers/particle distances and detect broken contacts
	////////////////////////////////////////////////////////////////////////////////////////////////////////

    updateContactNormalsKernel<<< numBlocks, numThreads >>>(
                                                              contact_particle_ids,
                                                              next_free_contact_id,
                                                              free_contacts_list,
                                                              number_of_new_contacts,
                                                              update_local_contact_list,
                                                              contact_normals_new,
                                                              positions,
                                                              number_of_particles);
    cudaDeviceSynchronize();
    getLastCudaError("updateContactNormalsKernel failed");

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// delete broken contacts from local contact lists
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

	int broken_contacts;
	checkCudaErrors(cudaMemcpy(&broken_contacts, number_of_new_contacts, sizeof(int), cudaMemcpyDeviceToHost));

	if(broken_contacts > 0)
	{
        total_contacts_broken += broken_contacts;

        deleteBrokenContactsKernel<<< numBlocks, numThreads >>>(
                                                              #ifdef GPU_TRACK_DISSIPATED_ENERGY
                                                                  dissipated_contact_energy,
                                                                  contact_normals_new,
                                                                  n_1,
                                                                  n_2,
                                                              #endif // GPU_TRACK_DISSIPATED_ENERGY
                                                                  particle_number_of_contacts,
                                                                  particle_particle_ids,
                                                                  particle_contact_ids,
                                                                  contact_particle_ids,
                                                                  update_local_contact_list,
                                                                  number_of_particles
                                                                  );
        cudaDeviceSynchronize();
		getLastCudaError("deleteBrokenContactsKernel failed");


		checkCudaErrors(cudaMemset(number_of_new_contacts, 0, sizeof(int)));
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// check for new contacts
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    updateStickingKernel<<< numBlocks, numThreads >>>(
                                                        new_contact_particle_ids,
                                                        number_of_new_contacts,
                                                        particle_particle_ids,
                                                        cell_start,
                                                        cell_end,
                                                        grid_particle_index,
                                                        positions_sorted,
                                                        number_of_particles
		);

    int new_contacts;
    checkCudaErrors(cudaMemcpy(&new_contacts, number_of_new_contacts, sizeof(int), cudaMemcpyDeviceToHost));
    total_contacts_created += new_contacts;

    cudaDeviceSynchronize();
	getLastCudaError("updateStickingKernel failed");

}

extern "C"
void cudaAddNewContacts(
        int* particle_number_of_contacts,
        int* particle_particle_ids,
        int* particle_contact_ids,
        int* next_free_contact_id,
        int* free_contacts_list,
        int2* new_contact_particle_ids,
        const int number_of_new_contacts,
        int2* contact_particle_ids,
        double4* contact_normals_new,
        double4* contact_normals_old,
        double4* rot1,
        double4* rot2,
        double3* n1_initial,
        double3* n2_initial,
        double4* n1,
        double* positions,
        const int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_new_contacts, BLOCK_SIZE, numBlocks, numThreads);

    addNewContactsKernel<<< numBlocks, numThreads >>>(
                                                        particle_number_of_contacts,
                                                        particle_particle_ids,
                                                        particle_contact_ids,
                                                        next_free_contact_id,
                                                        free_contacts_list,
                                                        new_contact_particle_ids,
                                                        number_of_new_contacts,
                                                        contact_particle_ids,
                                                        contact_normals_new,
                                                        contact_normals_old,
                                                        rot1,
                                                        rot2,
                                                        n1_initial,
                                                        n2_initial,
                                                        n1,
                                                        positions,
                                                        number_of_particles
                                                        );
    cudaDeviceSynchronize();
	getLastCudaError("addNewContactsKernel failed");

}

extern "C"
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
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(6 * number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    updateContactsKernel<<< numBlocks, numThreads >>>(
                                                        contact_particle_ids,
                                                        rot1,
                                                        rot2,
                                                        n1_initial,
                                                        n2_initial,
                                                        n_1,
                                                        n_2,
                                                        contact_normals_new,
                                                        contact_normals_old,
                                                        angular_velocities,
                                                        torques,
														timestep,
                                                        number_of_particles
                                                    #ifdef GPU_TRACK_DISSIPATED_ENERGY
                                                        ,
                                                        dissipated_rolling_energy,
                                                        dissipated_sliding_energy,
                                                        dissipated_twisting_energy
                                                    #endif // GPU_TRACK_DISSIPATED_ENERGY
                                                        );
    cudaDeviceSynchronize();
	getLastCudaError("updateContactsKernel failed");

}

extern "C"
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
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    updateInteractionKernel<<< numBlocks, numThreads >>>(
                                                           forces_new,
                                                           torques_new,
                                                           particle_number_of_contacts,
                                                           contact_ids,
                                                           contact_particle_ids,
                                                           n_1,
                                                           n_2,
                                                           contact_normals_new,
                                                           number_of_particles
                                                       #ifdef GPU_TRACK_DISSIPATED_ENERGY
                                                           ,
                                                           dissipated_damping_energy
                                                       #endif // GPU_TRACK_DISSIPATED_ENERGY
                                                           );
    cudaDeviceSynchronize();
	getLastCudaError("updateInteractionKernel failed");

}

extern "C"
void cudaUpdateBoxInteraction(
        double *forces_new,
        double *positions,
        const double3 lower_pos,
        const double3 upper_pos,
		double* top_wall_force,
		double* bot_wall_force,
        const unsigned int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    updateBoxInteractionKernel<<< numBlocks, numThreads >>>(
                                                              forces_new,
                                                              positions,
                                                              lower_pos,
                                                              upper_pos,
															  top_wall_force,
															  bot_wall_force,
                                                              number_of_particles
                                                              );
	//cudaDeviceSynchronize();
	getLastCudaError("updateBoxInteractionKernel failed");
}

extern "C"
void cudaUpdateOpenBoxInteraction(
        double *forces_new,
        double *positions,
        const double3 lower_pos,
        const double3 upper_pos,
        const unsigned int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    updateOpenBoxInteractionKernel<<< numBlocks, numThreads >>>(
                                                                  forces_new,
                                                                  positions,
                                                                  lower_pos,
                                                                  upper_pos,
                                                                  number_of_particles
                                                                  );
	//cudaDeviceSynchronize();
	getLastCudaError("updateOpenBoxInteractionKernel failed");
}


extern "C"
void cudaUpdateEnclosingSphereWallInteractionKernel(
        double *forces_new,
        double *positions,
        const double sphere_radius,
        const unsigned int number_of_particles
        )
{
    unsigned int numThreads, numBlocks;
    computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    updateEnclosingSphereWallInteractionKernel<<< numBlocks, numThreads >>>(
                                                                  forces_new,
                                                                  positions,
                                                                  sphere_radius,
                                                                  number_of_particles
                                                                  );
    //cudaDeviceSynchronize();
    getLastCudaError("updateEnclosingSphereWallInteractionKernel failed");
}


extern "C"
void cudaUpdateBottomWallInteraction(
        double *forces_new,
        double *positions,
        const double3 lower_pos,
        const unsigned int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    updateBottomWallInteractionKernel<<< numBlocks, numThreads >>>(
                                                                     forces_new,
                                                                     positions,
                                                                     lower_pos,
                                                                     number_of_particles
                                                                     );
	//cudaDeviceSynchronize();
	getLastCudaError("updateOpenBoxInteractionKernel failed");
}

#ifdef GPU_TRACK_PARTICLE_ORIENTATION
extern "C"
void cudaUpdateParticleOrientation(
        double4 *particle_orientations,
        double *angular_velocities,
        double *torques,
        const double timestep,
        const unsigned int number_of_particles
        )
{
	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    updateParticleOrientationKernel<<< numBlocks, numThreads >>>(
                                                                   particle_orientations,
                                                                   angular_velocities,
                                                                   torques,
                                                                   timestep,
                                                                   number_of_particles);
	//cudaDeviceSynchronize();
	getLastCudaError("updateParticleOrientationKernel failed");
}
#endif // GPU_TRACK_PARTICLE_ORIENTATION

#ifdef GPU_TRACK_DISSIPATED_ENERGY



#include "myCubWrapper.h"
extern "C"
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
        )
{


    diss_normal = 0.5 * timestep * ENERGY_UNIT * sortPlan.reduce(gpu_dissipated_damping_energy, number_of_particles); // factor 0.5 because GPU version does particle interaction for both contact partners induvidually and thus adds the damping energy twice
    diss_rolling = k_r * crit_rolling_displacement * ENERGY_UNIT * sortPlan.reduce(gpu_dissipated_rolling_energy, 6*number_of_particles);
    diss_sliding = k_s * crit_sliding_displacement * ENERGY_UNIT * sortPlan.reduce(gpu_dissipated_sliding_energy, 6*number_of_particles);
    diss_twisting = k_t * crit_twisting_displacement * ENERGY_UNIT * sortPlan.reduce(gpu_dissipated_twisting_energy, 6*number_of_particles);
    diss_contact = ENERGY_UNIT * sortPlan.reduce(gpu_dissipated_contact_energy, 6*number_of_particles);

    int half_number_of_particles = number_of_particles/2;
    v1.x = sortPlan.reduce(gpu_velocities, half_number_of_particles) / double(half_number_of_particles);
    v2.x = sortPlan.reduce(gpu_velocities+half_number_of_particles, half_number_of_particles) / double(half_number_of_particles);

    v1.y = sortPlan.reduce(gpu_velocities + number_of_particles, half_number_of_particles) / double(half_number_of_particles);
    v2.y = sortPlan.reduce(gpu_velocities + number_of_particles + half_number_of_particles, half_number_of_particles) / double(half_number_of_particles);

    v1.z = sortPlan.reduce(gpu_velocities + 2 * number_of_particles, half_number_of_particles) / double(half_number_of_particles);
    v2.z = sortPlan.reduce(gpu_velocities + 2 * number_of_particles + half_number_of_particles, half_number_of_particles) / double(half_number_of_particles);



	unsigned int numThreads, numBlocks;
	computeGridSize(number_of_particles, BLOCK_SIZE, numBlocks, numThreads);


    calculateKineticEnergyKernel<<< numBlocks, numThreads >>>(
                                                                 kinetic_energies,
                                                                 velocities,
                                                                 number_of_particles
                                                                 );
    cudaDeviceSynchronize();
	getLastCudaError("calculateParticleEnergyKernel failed");

    E_kin = ENERGY_UNIT * sortPlan.reduce(kinetic_energies, number_of_particles);



    calculateRollingEnergyKernel<<< numBlocks, numThreads >>>(
                                                                 rotation_energies,
                                                                 angular_velocities,
                                                                 number_of_particles
                                                                 );
    cudaDeviceSynchronize();
    getLastCudaError("calculateParticleEnergyKernel failed");

    E_rot = ENERGY_UNIT * sortPlan.reduce(rotation_energies, number_of_particles);


	computeGridSize(6*number_of_particles, BLOCK_SIZE, numBlocks, numThreads);

    calculateContactEnergyKernel<<< numBlocks, numThreads >>>(
                                                                normal_energies,
                                                                rolling_energies,
                                                                sliding_energies,
                                                                twisting_energies,
                                                                contact_particle_ids,
                                                                contact_normals,
                                                                n_1,
                                                                n_2,
                                                                number_of_particles
                                                                );
    cudaDeviceSynchronize();
	getLastCudaError("calculateContactEnergyKernel failed");

    V_normal = ENERGY_UNIT * sortPlan.reduce(normal_energies, 6*number_of_particles);
    V_rolling = ENERGY_UNIT * sortPlan.reduce(rolling_energies, 6*number_of_particles);
    V_sliding = ENERGY_UNIT * sortPlan.reduce(sliding_energies, 6*number_of_particles);
    V_twisting = ENERGY_UNIT * sortPlan.reduce(twisting_energies, 6*number_of_particles);


    V_tot = V_normal + V_rolling + V_sliding + V_twisting;

    E_tot = E_kin + E_rot + V_tot + diss_normal + diss_rolling + diss_sliding + diss_twisting + diss_contact;






    return;
}
#endif // GPU_TRACK_DISSIPATED_ENERGY
