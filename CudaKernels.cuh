
#ifndef CUDA_KERNELS_H
#define CUDA_KERNELS_H


#include "CudaDefines.cuh"
#include "my_device_helper_math.h"

__device__ __constant__ GPUConstants gpu_constants;



inline __device__
int3 getGridCell(double3 pos)
{
	int3 grid_pos;

	grid_pos.x = floor((pos.x - gpu_constants.gpu_grid_shift) * gpu_constants.gpu_grid_cell_width_inv);
	grid_pos.y = floor((pos.y - gpu_constants.gpu_grid_shift) * gpu_constants.gpu_grid_cell_width_inv);
	grid_pos.z = floor((pos.z - gpu_constants.gpu_grid_shift) * gpu_constants.gpu_grid_cell_width_inv);

	return grid_pos;
}

inline __device__
uint getGridHash(int3 grid_pos)
{
	// wrap grid
	grid_pos.x = grid_pos.x & (gpu_constants.gpu_num_grid_cells - 1);
	grid_pos.y = grid_pos.y & (gpu_constants.gpu_num_grid_cells - 1);
	grid_pos.z = grid_pos.z & (gpu_constants.gpu_num_grid_cells - 1);

	return  gpu_constants.gpu_num_grid_cells * ( grid_pos.z * gpu_constants.gpu_num_grid_cells + grid_pos.y ) + grid_pos.x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Integration
////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__
void posPredictionKernel(
	double* RESTRICT new_pos,
	const double* RESTRICT old_pos,
	const double* RESTRICT velocity,
	const double* RESTRICT force,
	const double timestep,
	unsigned int* RESTRICT grid_particle_cell,
	unsigned int* RESTRICT grid_particle_index,
	const unsigned int number_of_particles
        )
{
	const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size


	double3 pos = make_double3(old_pos, index, number_of_particles);
	double3 vel = make_double3(velocity, index, number_of_particles);
	double3 _force = make_double3(force, index, number_of_particles);

    // integrate position with: pos_new = pos_old + vel*dt + force*dt*dt/mass
	pos = fma(fma(_force, 0.5 * timestep * gpu_constants.mass_inv, vel), timestep, pos);

	store_double3(new_pos, index, number_of_particles, pos);

	// get address in grid
	int cell = getGridHash(getGridCell(pos));

	// store grid hash and particle index
	grid_particle_cell[index] = cell;
    grid_particle_index[index] = index;

}

__global__
void correctionKernel(
	double* RESTRICT velocity,
	double* RESTRICT angular_velocity,
	const double* RESTRICT old_force,
	const double* RESTRICT new_force,
	const double* RESTRICT old_torque,
	const double* RESTRICT new_torque,
	const double timestep,
	const unsigned int number_of_particles
        )
{
	const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size


	// load data
	double3 vel = make_double3(velocity, index, number_of_particles);
	double3 oldForce = make_double3(old_force, index, number_of_particles);
	double3 newForce = make_double3(new_force, index, number_of_particles);

	// update velocity
    vel += (0.5 * timestep * gpu_constants.mass_inv) * (newForce + oldForce);

	// store data
	store_double3(velocity, index, number_of_particles, vel);

	// load data
	double3 angularVel = make_double3(angular_velocity, index, number_of_particles);
	double3 oldTorque = make_double3(old_torque, index, number_of_particles);
	double3 newTorque = make_double3(new_torque, index, number_of_particles);


	// update angular velocity
    angularVel += (0.5 * timestep * gpu_constants.moment_of_inertia_inv) * (newTorque + oldTorque);

	// store data
	store_double3(angular_velocity, index, number_of_particles, angularVel);

	return;
}

__global__
void dampingCorrectionKernel(
	double* RESTRICT velocity,
	double* RESTRICT angular_velocity,
	const double* RESTRICT old_force,
	const double* RESTRICT new_force,
	const double* RESTRICT old_torque,
	const double* RESTRICT new_torque,
	const double damping_factor,
	const double timestep,
	const unsigned int number_of_particles
        )
{
	const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size
	
	// load data
	double3 vel = make_double3(velocity, index, number_of_particles);
	double3 oldForce = make_double3(old_force, index, number_of_particles);
	double3 newForce = make_double3(new_force, index, number_of_particles);

	// update velocity
	vel = damping_factor * (vel + (0.5 * timestep * gpu_constants.mass_inv) * (oldForce + newForce));

	// store data
	store_double3(velocity, index, number_of_particles, vel);



	// load data
	double3 angularVel = make_double3(angular_velocity, index, number_of_particles);
	double3 oldTorque = make_double3(old_torque, index, number_of_particles);
	double3 newTorque = make_double3(new_torque, index, number_of_particles);

	// update angular velocity
	angularVel = damping_factor * (angularVel + (0.5 * timestep * gpu_constants.moment_of_inertia_inv) * (oldTorque + newTorque));

	// store data
	store_double3(angular_velocity, index, number_of_particles, angularVel);

	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Grid
////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__
void findCellStartKernel(
	uint* RESTRICT cell_start,
	uint* RESTRICT cell_end,
	const uint* RESTRICT grid_particle_hash,
	const uint* RESTRICT grid_particle_index,
	double* RESTRICT positions_new,
	double* RESTRICT positions_sorted,
	const unsigned int number_of_particles
        )
{
	extern __shared__ uint shared_hash[];    // blockSize + 1 elements

    const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	uint hash;

    if(index < number_of_particles) 
	{

        hash = grid_particle_hash[index];

        // Load hash data into shared memory so that we can look 
        // at neighboring particle's hash value without loading
        // two hash values per thread
	    shared_hash[threadIdx.x+1] = hash;

		// first thread in block must load neighbor particle hash
	    if(index > 0 && threadIdx.x == 0)
		    shared_hash[0] = grid_particle_hash[index-1];
	}

	__syncthreads();
	
	if(index < number_of_particles) 
	{
		// If this particle has a different cell index to the previous
		// particle then it must be the first particle in the cell,
		// so store the index of this particle in the cell.
		// As it isn't the first particle, it must also be the cell end of
		// the previous particle's cell

	    if(index == 0 || hash != shared_hash[threadIdx.x])
	    {
		    cell_start[hash] = index;

            if(index > 0)
                cell_end[shared_hash[threadIdx.x]] = index;
	    }

        if (index == number_of_particles - 1)
            cell_end[hash] = index + 1;

        // Now use the sorted index to reorder the position data
		const uint sortedIndex = grid_particle_index[index];

		const double3 pos = make_double3(positions_new, sortedIndex, number_of_particles);

        store_double3(positions_sorted, index, number_of_particles, pos);

        return;

	}

}



////////////////////////////////////////////////////////////////////////////////////////////////////////
// Particle Interaction
////////////////////////////////////////////////////////////////////////////////////////////////////////

__device__
double getApproxContactRadius(const double compression)
{
    double a = gpu_constants.get_contact_radius_c1 + gpu_constants.get_contact_radius_c2*sqrt(compression);
	return sqrt(a);
}

__device__
double getContactRadius(const double compression_length)
{
    // contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
    // use equilibrium contact radius as starting value
#ifdef CPU_EQUIVALENT
    const double k = gpu_constants.reduced_radius * compression_length / 3.0;
#else
    const double k = gpu_constants.particle_radius * compression_length * 0.166666666666666666666666666666;
#endif
	double x_pow3;
	double x_new;
    double compression = gpu_constants.delta_c + compression_length;
    double x_old = 0.0;

    if(compression  >= 0.0)
    {
#ifdef CPU_EQUIVALENT
        x_old = gpu_constants.c1_contact_radius;
#else
        x_old = getApproxContactRadius(compression);
#endif
    }
    else
    {
        printf("Error getContactRadius: too long bond not correctly broken\n");
        return 0.0;
    }

    // use Newton-Raphson method to find root
    for(int i = 0; i < 20; ++i)
                                
	{
		x_pow3 = x_old * x_old * x_old;
        x_new = 0.75 * (__dadd_rn(__dmul_rn(x_pow3, x_old), k)) / (x_pow3 - gpu_constants.c2_contact_radius);

        if (fabs(x_new - x_old) / x_new < 1.e-14)
            break;

		x_old = x_new;
	}



    return x_new * x_new;
}

__device__
double getNormalForce(double contact_radius)
{
	double x = sqrt(contact_radius*contact_radius*contact_radius);

    return x * (gpu_constants.c1_normal_force * x + gpu_constants.c2_normal_force);

}

__device__
double getNormalHertzForce(double compression_length)
{
    // F_hertz = 4/3 * young_mod_reduced * sqrt(reduced_particle_radius) * compression_length**(3/2)
    double F = gpu_constants.k_hertz * compression_length * sqrt(compression_length);

    return F;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contact Management
////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__
void updateContactNormalsKernel(
	int2* RESTRICT contact_particle_ids,
	int* RESTRICT next_free_contact_id,
	int* RESTRICT free_contacts_list,
	int* RESTRICT number_of_broken_contacts,
	int* RESTRICT update_local_contact_list,
	double4* RESTRICT contact_normals_new,
    const double* RESTRICT positions,
	const int number_of_particles
        )
{
	const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;


	if(index >= 6 * number_of_particles)
		return;


      // determine particle ids
      int2 particle_ids = contact_particle_ids[index];

      if (particle_ids.x == NO_PARTICLE)
          return;


	//////////////////////////////////////////////////////////////////////////////
	// update contact normals
	//////////////////////////////////////////////////////////////////////////////
	double3 pos1 = make_double3(positions, particle_ids.x, number_of_particles);
	double3 pos2 = make_double3(positions, particle_ids.y, number_of_particles);


	double3 deltaPos = pos1 - pos2;

	double4 n_c = make_double4(deltaPos, dot(deltaPos, deltaPos));

	// normalize n_c and calculate particle distance d
	n_c.w = sqrt(n_c.w);

    double temp = 1.0 / n_c.w;
    n_c.x *= temp;
    n_c.y *= temp;
    n_c.z *= temp;

	contact_normals_new[index] = n_c;


	//////////////////////////////////////////////////////////////////////////////
	// check if contact has broken
	//////////////////////////////////////////////////////////////////////////////

	if (n_c.w > gpu_constants.contact_breaking_dist)
	{
		// tell particles to update their local contact lists
		update_local_contact_list[particle_ids.x] = 1;
		update_local_contact_list[particle_ids.y] = 1;

		// remove contact from global contact list
		contact_particle_ids[index] = make_int2(NO_PARTICLE, NO_PARTICLE);
		int id = atomicSub(next_free_contact_id, 1) - 1;
		free_contacts_list[id] = index;

        // all 4 lines do the job. Which one is the fastest though?
        atomicAdd(number_of_broken_contacts, 1);  // this line also helps tracking number of broken contacts
        //atomicCAS(number_of_broken_contacts, 0, 1);
        //if(*number_of_broken_contacts == 0) *number_of_broken_contacts = 1;
        //*number_of_broken_contacts = 1;

		return;
	}

}

#ifdef GPU_TRACK_DISSIPATED_ENERGY
__device__
void calculateDissipatedContactEnergyKernel(
        double * RESTRICT dissipated_contact_energy,
        double4* RESTRICT contact_normals,
        double4* RESTRICT n_1,
        double4* RESTRICT n_2,
        const int index
        )
{

    //////////////////////////////////////////////////////////////////////////////
    // determine current contact pointers
    //////////////////////////////////////////////////////////////////////////////

    double4 n_c = contact_normals[index];
    double4 n1 = n_1[index];
    double4 n2 = n_2[index];

    //////////////////////////////////////////////////////////////////////////////
    // calculate potential energy
    //////////////////////////////////////////////////////////////////////////////

    double3 displacement;

#ifdef GPU_ROLLING
    displacement = make_double3(n1 + n2);
    dissipated_contact_energy[index] += 0.5 * gpu_constants.k_r * dot(displacement, displacement);
#endif // GPU_ROLLING

#ifdef GPU_SLIDING
    displacement = (make_double3(n1 - n2 + 2.0 * n_c));

    double temp = dot(displacement, n_c);

    displacement -= temp * make_double3(n_c);

    dissipated_contact_energy[index] += 0.5 * gpu_constants.k_s * dot(displacement, displacement);

#endif // GPU_SLIDING

#ifdef GPU_TWISTING
    dissipated_contact_energy[index] += 0.5 * gpu_constants.k_t * n1.w * n1.w;
#endif // GPU_TWISTING

}
#endif // GPU_TRACK_DISSIPATED_ENERGY




__global__
void deleteBrokenContactsKernel(
#ifdef GPU_TRACK_DISSIPATED_ENERGY
    double * RESTRICT dissipated_contact_energy,
    double4* RESTRICT contact_normals,
    double4* RESTRICT n_1,
    double4* RESTRICT n_2,
#endif // GPU_TRACK_DISSIPATED_ENERGY
	int* RESTRICT particle_number_of_contacts,
	int* RESTRICT particle_particle_ids,
	int* RESTRICT particle_contact_ids,
	const int2* RESTRICT contact_particle_ids,
	int* RESTRICT update_local_contact_list,
	const int number_of_particles
        )
{
	const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// check if update of local contact list is necessary
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int update = update_local_contact_list[index];
	if(update == 0)
		return;

	update_local_contact_list[index] = 0;
	
	int my_contacts = particle_number_of_contacts[index];

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// check contact list
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int c = 0; c < MAX_CONTACTS; ++c)
	{
		int contact_id = particle_contact_ids[CONTACT_ID(c, index, number_of_particles)];

		if(contact_id >= 0)
		{
			if(contact_particle_ids[contact_id].x == NO_PARTICLE)
			{

#ifdef GPU_TRACK_DISSIPATED_ENERGY
                calculateDissipatedContactEnergyKernel(
                            dissipated_contact_energy,
                            contact_normals,
                            n_1,
                            n_2,
                            contact_id);
#endif // GPU_TRACK_DISSIPATED_ENERGY

				///////////////////////////////////////////////////////////////////////////////////////////////////////
				// remove contact from contact list of the particle
				///////////////////////////////////////////////////////////////////////////////////////////////////////

				int last_id = my_contacts - 1;

				// overwrite with data of last contact to avoid gaps in the contact list of the particle
				if(last_id != c)
				{
					particle_particle_ids[ CONTACT_ID(c, index, number_of_particles) ] = particle_particle_ids[CONTACT_ID(last_id, index, number_of_particles)];
					particle_contact_ids[CONTACT_ID(c, index, number_of_particles)] = particle_contact_ids[CONTACT_ID(last_id, index, number_of_particles)];
					--c; // decrement counter to make sure the copied contact is checked as well
				}

				particle_particle_ids[CONTACT_ID(last_id, index, number_of_particles)] = NO_PARTICLE;
				particle_contact_ids[CONTACT_ID(last_id, index, number_of_particles)] = -1;

				--my_contacts;
			}
		}
	}

	particle_number_of_contacts[index] = my_contacts;
}




__global__
void updateStickingKernel(
	int2* RESTRICT new_contact_particle_ids,
	int* RESTRICT number_of_new_contacts,
	const int* RESTRICT particle_particle_ids,
	const uint* RESTRICT cell_start,
	const uint* RESTRICT cell_end,
    const uint* RESTRICT grid_particle_index,
	const double* RESTRICT positions_sorted,
	const int number_of_particles
        )
{
    const uint tid = blockIdx.x * blockDim.x + threadIdx.x;

	if(tid >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size

    uint index =  __ldg(&grid_particle_index[tid]);

	// load position
    double3 pos = make_double3(positions_sorted, tid, number_of_particles);

	// load contact ids
	__shared__ int particle_ids[BLOCK_SIZE * MAX_CONTACTS];

	for(int c = 0; c < MAX_CONTACTS; ++c)
		particle_ids[c * BLOCK_SIZE + threadIdx.x] = particle_particle_ids[CONTACT_ID(c, index, number_of_particles)];


	int3 grid_particle_pos = getGridCell(pos);

	// scan neighbouring cells
	for(int z = -1; z <= 1; ++z)
	{
		for(int y = -1; y <= 1; ++y)
		{
			for(int x = -1; x <= 1; ++x)
			{

				uint cell = getGridHash(grid_particle_pos + make_int3(x, y, z));


				// get index of the first particle in this cell
				uint start_index = cell_start[cell];


				// check if cell is not empty
				if(start_index != NO_PARTICLE) 
				{
					uint end_index = cell_end[cell];

                    for(unsigned int j = start_index; j < end_index && j != tid; ++j)
					{
                        int index2 = __ldg(&grid_particle_index[j]);

						if (index <= index2)
							continue;

						// look for the id of the other particle in the contact list
						bool contact_list_pos  = true;
						for(int c = 0; c < MAX_CONTACTS; ++c)
						{
							if( particle_ids[c * BLOCK_SIZE + threadIdx.x] == index2 )
							{
								contact_list_pos = false;
								break;
							}
						}

						if(contact_list_pos)
						{
							// determine distance to other particle
                            double3 pos2 = make_double3(positions_sorted, j, number_of_particles);
                            double3 delta_pos = pos - pos2;

							double dist_squared = dot(delta_pos, delta_pos);


							///////////////////////////////////////////////////////////////////////////////////////////////////////
							// add new contact if particles are close enough
							///////////////////////////////////////////////////////////////////////////////////////////////////////

                            if(dist_squared < gpu_constants.contact_making_dist_squared)    // contact if delta >= 0
                                                                                            // delta = r1+r2 - |x1-x2|                                                                                     // delta >= 0  <=> |x1-x2|^2 <= (r1+r2)^2
                            {
                                // put contact on the list of contacts that will be added by a separate kernel
								int id = atomicAdd(number_of_new_contacts, 1);
								new_contact_particle_ids[id] = make_int2(index, index2);

							}
						}
					}
				}
			}
		}
	}
}





__global__
void addNewContactsKernel(
	int* RESTRICT particle_number_of_contacts,
	int* RESTRICT particle_particle_ids,
	int* RESTRICT particle_contact_ids,
	int* RESTRICT  next_free_contact_id,
	const int* RESTRICT free_contacts_list,
	const int2* RESTRICT new_contact_particle_ids,
	const int number_of_new_contacts,
	int2* RESTRICT contact_particle_ids,
	double4* RESTRICT contact_normals_new,
	double4* RESTRICT contact_normals_old,
	double4* RESTRICT rot1,
	double4* RESTRICT rot2,
    double3* RESTRICT n1_initial,
    double3* RESTRICT n2_initial,
	double4* RESTRICT n1,
    const double* RESTRICT positions,
	const int number_of_particles
        )
{
	const unsigned int index = blockIdx.x *blockDim.x + threadIdx.x;

	if(index >= number_of_new_contacts)
		return;


	int2 particle_ids = new_contact_particle_ids[index];

	int contact_list_pos = atomicAdd(next_free_contact_id, 1);
	int global_contact_id = free_contacts_list[contact_list_pos];

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// add new contact to global contact list
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	// determine n_c/n1/n2
	double3 pos1 = make_double3(positions, particle_ids.x, number_of_particles);
	double3 pos2 = make_double3(positions, particle_ids.y, number_of_particles);

	double3 delta_pos = pos1 - pos2;

	double dist_squared = dot(delta_pos, delta_pos);
	double dist = rsqrt(dist_squared);
	delta_pos *= dist;
    dist *= dist_squared;

	contact_particle_ids[global_contact_id] = particle_ids;
	rot1[global_contact_id] = make_double4(0.0, 0.0, 0.0, 1.0);
	rot2[global_contact_id] = make_double4(0.0, 0.0, 0.0, 1.0);
    n1_initial[global_contact_id] = make_double3(-delta_pos.x, -delta_pos.y, -delta_pos.z);
    n2_initial[global_contact_id] = make_double3(delta_pos.x, delta_pos.y, delta_pos.z);
    contact_normals_old[global_contact_id] = make_double4(delta_pos.x, delta_pos.y, delta_pos.z, 2.0*gpu_constants.particle_radius);
	contact_normals_new[global_contact_id] = make_double4(delta_pos.x, delta_pos.y, delta_pos.z, dist);
	n1[global_contact_id] = make_double4(0.0, 0.0, 0.0, 0.0);

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// add new contact to local contact list of the particles
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	// first particle
	int local_contact_id = atomicAdd( &(particle_number_of_contacts[particle_ids.x]), 1);
	particle_particle_ids[CONTACT_ID(local_contact_id, particle_ids.x, number_of_particles)] = particle_ids.y;
	particle_contact_ids[CONTACT_ID(local_contact_id, particle_ids.x, number_of_particles)] = global_contact_id;

	// second particle
	local_contact_id = atomicAdd( &(particle_number_of_contacts[particle_ids.y]), 1);
	particle_particle_ids[CONTACT_ID(local_contact_id, particle_ids.y, number_of_particles)] = particle_ids.x;
	particle_contact_ids[CONTACT_ID(local_contact_id, particle_ids.y, number_of_particles)] = global_contact_id;

	return;
}






__global__
void updateContactsKernel(
    const int2* RESTRICT contact_particle_ids,
    double4* rot1,
    double4* rot2,
    double3* n1_initial,
    double3* n2_initial,
    double4* n_1,
    double4* n_2,
    const double4* RESTRICT contact_normals_new,
    const double4* RESTRICT contact_normals_old,
    const double* RESTRICT angular_velocities,
    const double* RESTRICT torques,
	const double timestep,
	const uint number_of_particles
#ifdef GPU_TRACK_DISSIPATED_ENERGY
    ,
    double* RESTRICT        dissipated_rolling_energy,
    double* RESTRICT        dissipated_sliding_energy,
    double* RESTRICT        dissipated_twisting_energy
#endif // GPU_TRACK_DISSIPATED_ENERGY
	)
{
	const unsigned int index = blockIdx.x *blockDim.x + threadIdx.x;

	if (index >= 6 * number_of_particles)
		return;

	// determine particle ids
	int2 particle_ids = contact_particle_ids[index];

	if (particle_ids.x == NO_PARTICLE)
		return;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load torque and angular velocity
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double3 omega1 = make_double3(angular_velocities, particle_ids.x, number_of_particles);
	double3 omega2 = make_double3(angular_velocities, particle_ids.y, number_of_particles);

    double3 omega1_dot = gpu_constants.moment_of_inertia_inv * make_double3(torques, particle_ids.x, number_of_particles);
    double3 omega2_dot = gpu_constants.moment_of_inertia_inv * make_double3(torques, particle_ids.y, number_of_particles);

	//////////////////////////////////////////////////////////////////////////////
	// update twisting displacement
	//////////////////////////////////////////////////////////////////////////////

	// use second order integration: Phi^n+1 = Phi^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
	double4 n_c = contact_normals_new[index];
	double4 n_c_old = contact_normals_old[index];


#ifdef GPU_TWISTING
    double3 delta_omega_old;
    delta_omega_old = omega1 - omega2;


    double3 delta_omega_new;
    delta_omega_new = delta_omega_old + timestep * (omega1_dot - omega2_dot);

#ifdef CPU_EQUIVALENT
    double twisting_displacement = n_1[index].w + 0.5 * timestep * (dot(delta_omega_old, make_double3(n_c_old)) + dot(delta_omega_new, make_double3(n_c)));
#else
    double twisting_displacement = n_1[index].w + 0.5 * timestep * (
                fma(delta_omega_old.x, n_c_old.x,
                    fma(delta_omega_old.y, n_c_old.y,
                        fma(delta_omega_old.z, n_c_old.z,
                            fma(delta_omega_new.x, n_c.x,
                                fma(delta_omega_new.y, n_c.y, delta_omega_new.z * n_c.z))))));
#endif // CPU_EQUIVALENT



#endif // GPU_TWISTING


#ifdef GPU_INELASTIC_TWISTING
    if (twisting_displacement > gpu_constants.crit_twisting_displacement)
    {

#ifdef GPU_TRACK_DISSIPATED_ENERGY
        dissipated_twisting_energy[index] += (twisting_displacement - gpu_constants.crit_twisting_displacement); // missing: k_t * crit_twisting_displacement * ENERGY_UNIT
#endif// GPU_TRACK_DISSIPATED_ENERGY

        twisting_displacement = gpu_constants.crit_twisting_displacement;

    }
#endif

	// cache to store later
	double v_rel = (n_c.w - n_c_old.w) / timestep;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// update rotation parameters
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double4 e_dot;
	double4 e_ddot;
	double4 rot_param1 = rot1[index];

    e_dot.w = -0.5 * dot(make_double3(rot_param1), omega1);
    e_dot.x = 0.5 * (diff_of_products_simple(rot_param1.w, omega1.x, rot_param1.y, omega1.z) + rot_param1.z * omega1.y);
    e_dot.y = 0.5 * (diff_of_products_simple(rot_param1.w, omega1.y, rot_param1.z, omega1.x) + rot_param1.x * omega1.z);
    e_dot.z = 0.5 * (diff_of_products_simple(rot_param1.w, omega1.z, rot_param1.x, omega1.y) + rot_param1.y * omega1.x);

	double temp = 0.5 * e_dot.w;

    e_ddot.w = -0.25 * (rot_param1.w * dot(omega1, omega1) + 2.0 * dot(rot_param1, omega1_dot));
    e_ddot.x = temp * omega1.x + 0.5 * (diff_of_products_simple(rot_param1.w, omega1_dot.x, rot_param1.y, omega1_dot.z) + rot_param1.z * omega1_dot.y);
    e_ddot.y = temp * omega1.y + 0.5 * (diff_of_products_simple(rot_param1.w, omega1_dot.y, rot_param1.z, omega1_dot.x) + rot_param1.x * omega1_dot.z);
    e_ddot.z = temp * omega1.z + 0.5 * (diff_of_products_simple(rot_param1.w, omega1_dot.z, rot_param1.x, omega1_dot.y) + rot_param1.y * omega1_dot.x);

	rot_param1 += timestep * e_dot + 0.5 * timestep * timestep * e_ddot;


	// normalize
	temp = rsqrt(rot_param1.w*rot_param1.w + rot_param1.x*rot_param1.x + rot_param1.y*rot_param1.y + rot_param1.z*rot_param1.z);
	rot_param1 *= temp;


	rot1[index] = rot_param1;

	//////////////////////////////////////////////////////////////////////////////
	// update rotation parameters
	//////////////////////////////////////////////////////////////////////////////

	double4 rot_param2 = rot2[index];

    e_dot.w = -0.5 * dot(make_double3(rot_param2), omega2);
    e_dot.x = 0.5 * (rot_param2.w * omega2.x - rot_param2.y * omega2.z + rot_param2.z * omega2.y);
    e_dot.y = 0.5 * (rot_param2.w * omega2.y - rot_param2.z * omega2.x + rot_param2.x * omega2.z);
    e_dot.z = 0.5 * (rot_param2.w * omega2.z - rot_param2.x * omega2.y + rot_param2.y * omega2.x);

	temp = 0.5 * e_dot.w;

	e_ddot.w = -0.25 * (rot_param2.w * dot(omega2, omega2) + 2.0 * dot(rot_param2,  omega2_dot));
    e_ddot.x = temp * omega2.x + 0.5 * (rot_param2.w * omega2_dot.x - rot_param2.y * omega2_dot.z + rot_param2.z * omega2_dot.y);
    e_ddot.y = temp * omega2.y + 0.5 * (rot_param2.w * omega2_dot.y - rot_param2.z * omega2_dot.x + rot_param2.x * omega2_dot.z);
    e_ddot.z = temp * omega2.z + 0.5 * (rot_param2.w * omega2_dot.z - rot_param2.x * omega2_dot.y + rot_param2.y * omega2_dot.x);

	rot_param2 += timestep * e_dot + 0.5 * timestep * timestep * e_ddot;


	// normalize
	temp = rsqrt(rot_param2.w*rot_param2.w + rot_param2.x*rot_param2.x + rot_param2.y*rot_param2.y + rot_param2.z*rot_param2.z);
	rot_param2 *= temp;


	rot2[index] = rot_param2;

	//////////////////////////////////////////////////////////////////////////////
	// determine current contact pointers
	//////////////////////////////////////////////////////////////////////////////

	double4 n1;
#ifdef GPU_TWISTING
	n1.w = twisting_displacement;
#endif // GPU_TWISTING
    double3 n_initial = n1_initial[index];
    n1.x = 2.0 * ((rot_param1.w * rot_param1.w + rot_param1.x * rot_param1.x - 0.5) * n_initial.x + diff_of_products_simple(rot_param1.x, rot_param1.y, rot_param1.z, rot_param1.w) * n_initial.y + (rot_param1.x * rot_param1.z + rot_param1.y * rot_param1.w) * n_initial.z);
    n1.y = 2.0 * ((rot_param1.x * rot_param1.y + rot_param1.z * rot_param1.w) * n_initial.x + (rot_param1.w * rot_param1.w + rot_param1.y * rot_param1.y - 0.5) * n_initial.y + diff_of_products_simple(rot_param1.y, rot_param1.z, rot_param1.x, rot_param1.w) * n_initial.z);
    n1.z = 2.0 * (diff_of_products_simple(rot_param1.x, rot_param1.z, rot_param1.y, rot_param1.w) * n_initial.x + (rot_param1.y * rot_param1.z + rot_param1.x * rot_param1.w) * n_initial.y + (rot_param1.w * rot_param1.w + rot_param1.z * rot_param1.z - 0.5) * n_initial.z);

	double4 n2;
	n2.w = v_rel;
	n_initial = n2_initial[index];
    n2.x = 2.0 * ((rot_param2.w * rot_param2.w + rot_param2.x * rot_param2.x - 0.5) * n_initial.x + diff_of_products_simple(rot_param2.x, rot_param2.y, rot_param2.z, rot_param2.w) * n_initial.y + (rot_param2.x * rot_param2.z + rot_param2.y * rot_param2.w) * n_initial.z);
    n2.y = 2.0 * ((rot_param2.x * rot_param2.y + rot_param2.z * rot_param2.w) * n_initial.x + (rot_param2.w * rot_param2.w + rot_param2.y * rot_param2.y - 0.5) * n_initial.y + diff_of_products_simple(rot_param2.y, rot_param2.z, rot_param2.x, rot_param2.w) * n_initial.z);
    n2.z = 2.0 * (diff_of_products_simple(rot_param2.x, rot_param2.z, rot_param2.y, rot_param2.w) * n_initial.x + (rot_param2.y * rot_param2.z + rot_param2.x * rot_param2.w) * n_initial.y + (rot_param2.w * rot_param2.w + rot_param2.z * rot_param2.z - 0.5) * n_initial.z);


	//////////////////////////////////////////////////////////////////////////////
	// apply inelastic corrections
	//////////////////////////////////////////////////////////////////////////////

	bool contact_pointers_modified = false;

#ifdef GPU_INELASTIC_SLIDING
	{
		double3 displacement;
		displacement = make_double3(n1 - n2);


		double temp = dot(displacement, n_c);

		displacement -= make_double3(temp * n_c);


		displacement *= gpu_constants.particle_radius;


        double displacement_norm = dot(displacement, displacement);

		// check if we are in the inelastic regime
        if (displacement_norm > gpu_constants.crit_sliding_displacement_squared)
		{
            displacement_norm = sqrt(displacement_norm);

            // determine correction of contact pointers
            // Problem: sometimes norm == gpu_constants.crit_sliding_displacement; reason for if statement below
            displacement = (1.0 - gpu_constants.crit_sliding_displacement / displacement_norm) * displacement;


            double norm = dot(displacement, displacement);

            if(norm > 0.0) // norm == 0 -> displacement == 0 -> division by 0 = nan -> crash
            {
                // correct contact pointers
                temp = dot(n1, displacement);
                double alpha = 1.0 / (1.0 - temp * temp / norm);

                n1.x -= 0.5 * gpu_constants.particle_radius_inv * alpha * displacement.x;
                n1.y -= 0.5 * gpu_constants.particle_radius_inv * alpha * displacement.y;
                n1.z -= 0.5 * gpu_constants.particle_radius_inv * alpha * displacement.z;

                temp = dot(n2, displacement);
                alpha = 1.0 / (1.0 - temp * temp / norm);

                n2.x += 0.5 * gpu_constants.particle_radius_inv * alpha * displacement.x;
                n2.y += 0.5 * gpu_constants.particle_radius_inv * alpha * displacement.y;
                n2.z += 0.5 * gpu_constants.particle_radius_inv * alpha * displacement.z;

                // normalize contact pointers
                norm = rsqrt(n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);
                n1.x *= norm;
                n1.y *= norm;
                n1.z *= norm;

                norm = rsqrt(n2.x*n2.x + n2.y*n2.y + n2.z*n2.z);
                n2.x *= norm;
                n2.y *= norm;
                n2.z *= norm;

                contact_pointers_modified = true;


#ifdef GPU_TRACK_DISSIPATED_ENERGY
                dissipated_sliding_energy[index] += (displacement_norm - gpu_constants.crit_sliding_displacement); // missing: k_s * crit_sliding_displacement * ENERGY_UNIT
#endif // GPU_TRACK_DISSIPATED_ENERGY

            }
		}
	}
#endif	// GPU_INELASTIC_SLIDING

#ifdef GPU_INELASTIC_ROLLING
	{
		double3 displacement;
		displacement = make_double3(0.5 * gpu_constants.particle_radius * (n1 + n2));


        double displacement_norm = dot(displacement, displacement);


        if (displacement_norm > gpu_constants.crit_rolling_displacement_squared)
		{
            displacement_norm = sqrt(displacement_norm);

			// determine correction of contact pointers
            // Problem: sometimes norm == gpu_constants.crit_rolling_displacement; reason for if statement below
            displacement *= (1.0 - gpu_constants.crit_rolling_displacement / displacement_norm);


            double norm = dot(displacement, displacement);

            if(norm > 0.0) // norm == 0 -> displacement == 0 -> division by 0 = nan -> crash
            {

                // correct contact pointers
                // calculate correction factor alpha (see Wada et al. 2007 appendix for details)

                double temp = dot(n1, displacement);
                double alpha = 1.0 / (1.0 - temp * temp / norm);


                n1.x -= gpu_constants.particle_radius_inv * alpha * displacement.x;
                n1.y -= gpu_constants.particle_radius_inv * alpha * displacement.y;
                n1.z -= gpu_constants.particle_radius_inv * alpha * displacement.z;


                temp = dot(n2, displacement);
                alpha = 1.0 / (1.0 - temp * temp / norm);

                n2.x -= gpu_constants.particle_radius_inv * alpha * displacement.x;
                n2.y -= gpu_constants.particle_radius_inv * alpha * displacement.y;
                n2.z -= gpu_constants.particle_radius_inv * alpha * displacement.z;


                // normalize contact pointers
                norm = rsqrt(n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);
                n1.x *= norm;
                n1.y *= norm;
                n1.z *= norm;

                norm = rsqrt(n2.x*n2.x + n2.y*n2.y + n2.z*n2.z);
                n2.x *= norm;
                n2.y *= norm;
                n2.z *= norm;

                contact_pointers_modified = true;


#ifdef GPU_TRACK_DISSIPATED_ENERGY
                dissipated_rolling_energy[index] += (displacement_norm - gpu_constants.crit_rolling_displacement); // missing: k_r * crit_rolling_displacement * ENERGY_UNIT
#endif // GPU_TRACK_DISSIPATED_ENERGY

            }
		}
	}
#endif	// GPU_INELASTIC_ROLLING

#if defined(GPU_INELASTIC_ROLLING) || defined(GPU_INELASTIC_SLIDING)
	if (contact_pointers_modified)
	{
		// update initial contact pointers
        double3 n1_initial_buffer;
		double4 rot_param1 = rot1[index];
        n1_initial_buffer.x = 2.0 * ((0.5 - rot_param1.y * rot_param1.y - rot_param1.z * rot_param1.z) * n1.x + (rot_param1.x * rot_param1.y + rot_param1.z * rot_param1.w) * n1.y + diff_of_products_simple(rot_param1.x, rot_param1.z, rot_param1.y, rot_param1.w) * n1.z);
        n1_initial_buffer.y = 2.0 * (diff_of_products_simple(rot_param1.x, rot_param1.y, rot_param1.z, rot_param1.w) * n1.x + (0.5 - rot_param1.x * rot_param1.x - rot_param1.z * rot_param1.z) * n1.y + (rot_param1.y * rot_param1.z + rot_param1.x * rot_param1.w) * n1.z);
        n1_initial_buffer.z = 2.0 * ((rot_param1.x * rot_param1.z + rot_param1.y * rot_param1.w) * n1.x + diff_of_products_simple(rot_param1.y, rot_param1.z, rot_param1.x, rot_param1.w) * n1.y + (0.5 - rot_param1.x * rot_param1.x - rot_param1.y * rot_param1.y) * n1.z);
		n1_initial[index] = n1_initial_buffer;


        double3 n2_initial_buffer;
		double4 rot_param2 = rot2[index];
        n2_initial_buffer.x = 2.0 * ((0.5 - rot_param2.y * rot_param2.y - rot_param2.z * rot_param2.z) * n2.x + (rot_param2.x * rot_param2.y + rot_param2.z * rot_param2.w) * n2.y + diff_of_products_simple(rot_param2.x, rot_param2.z, rot_param2.y, rot_param2.w) * n2.z);
        n2_initial_buffer.y = 2.0 * (diff_of_products_simple(rot_param2.x, rot_param2.y, rot_param2.z, rot_param2.w) * n2.x + (0.5 - rot_param2.x * rot_param2.x - rot_param2.z * rot_param2.z) * n2.y + (rot_param2.y * rot_param2.z + rot_param2.x * rot_param2.w) * n2.z);
        n2_initial_buffer.z = 2.0 * ((rot_param2.x * rot_param2.z + rot_param2.y * rot_param2.w) * n2.x + diff_of_products_simple(rot_param2.y, rot_param2.z, rot_param2.x, rot_param2.w) * n2.y + (0.5 - rot_param2.x * rot_param2.x - rot_param2.y * rot_param2.y) * n2.z);
		n2_initial[index] = n2_initial_buffer;

	}
#endif


	// write back n1/n2
	n_1[index] = n1;
    n_2[index] = n2;

	return;
}



__global__
void updateInteractionKernel(
        double* RESTRICT forces_new,
        double* RESTRICT torques_new,
        const int* RESTRICT particle_number_of_contacts,
        const int* RESTRICT contact_ids,
        const int2* RESTRICT contact_particle_ids,
        const double4* RESTRICT n_1,
        const double4* RESTRICT n_2,
        const double4* RESTRICT contact_normals_new,
        const unsigned int number_of_particles
        #ifdef GPU_TRACK_DISSIPATED_ENERGY
            ,
            double* RESTRICT        dissipated_damping_energy
        #endif // GPU_TRACK_DISSIPATED_ENERGY
        )
{
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size

    double3 force = make_double3(0.0);
    double3 torque = make_double3(0.0);

	int num_contacts = particle_number_of_contacts[index];

	//for(int c = 0; c < MAX_CONTACTS; ++c)
	for(int c = 0; c < num_contacts; ++c)
	{
		// load contact ids
		int contact_id = contact_ids[CONTACT_ID(c, index, number_of_particles)];

		//if(contact_id >= 0)
		{
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// load data
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// load id of other particle
			int particle_id = contact_particle_ids[contact_id].y;

			// load n1, n2, n_c, phil, v_rel, dist
			double4 n1 = n_1[contact_id];
			double4 n2 = n_2[contact_id];
            double4 n_c = contact_normals_new[contact_id];

			// apply correct sign/change contact pointers if necessary (F_ij = - F_ji
			if(particle_id == index)
			{
				double temp;
				temp = n1.x;
				n1.x = n2.x;
				n2.x = temp;

				temp = n1.y;
				n1.y = n2.y;
				n2.y = temp;

				temp = n1.z;
				n1.z = n2.z;
				n2.z = temp;

				n_c.x *= -1.0;
				n_c.y *= -1.0;
				n_c.z *= -1.0;
			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// normal force (+ damping)
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            double compression_length = 2.0 * gpu_constants.particle_radius - n_c.w;

            double contact_radius = getContactRadius(compression_length);

            double interaction_strength = getNormalForce(contact_radius);

#ifdef GPU_USE_CONSTANT_DAMPING
			// relative velocity is stored in n2.w
            interaction_strength -= gpu_constants.osc_damping_factor * n2.w;

#ifdef GPU_TRACK_DISSIPATED_ENERGY
            dissipated_damping_energy[index] += gpu_constants.osc_damping_factor * n2.w * n2.w; // missing  *timestep*ENERGY_UNIT
#endif // GPU_TRACK_DISSIPATED_ENERGY
#endif // GPU_USE_CONSTANT_DAMPING

#ifdef GPU_USE_VISCOELASTIC_DAMPING
			// relative velocity is stored in n2.w
            double damping_strength = gpu_constants.osc_damping_factor * contact_radius;
			interaction_strength -= damping_strength * n2.w;

#ifdef GPU_TRACK_DISSIPATED_ENERGY
            dissipated_damping_energy[index] += damping_strength * n2.w * n2.w; // missing: timestep*ENERGY_UNIT
#endif // GPU_TRACK_DISSIPATED_ENERGY
#endif // GPU_USE_VISCOELASTIC_DAMPING

            force += interaction_strength * make_double3(n_c);


			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// sliding
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef GPU_SLIDING
            double3 delta_n = make_double3(n1 - n2);

            // sliding force
            double temp = dot(delta_n, make_double3(n_c));


            // temp_vec = delta_n - temp * n_c
#ifdef CPU_EQUIVALENT
            double3 temp_vec; // numerics equal to cpu version
            temp_vec.x = __dsub_rn(delta_n.x, __dmul_rn(temp, n_c.x));
            temp_vec.y = __dsub_rn(delta_n.y, __dmul_rn(temp, n_c.y));
            temp_vec.z = __dsub_rn(delta_n.z, __dmul_rn(temp, n_c.z));
#else
            double3 temp_vec = fma(make_double3(n_c), -temp, delta_n);
#endif
            // Note k_s is actually k_s * particle_radius ^2
            interaction_strength = gpu_constants.k_s / n_c.w * temp;
            force += interaction_strength * temp_vec;

            // sliding torque
#ifdef CPU_EQUIVALENT
            torque.x -= gpu_constants.k_s * diff_of_products_simple(n1.y, temp_vec.z, n1.z, temp_vec.y);
            torque.y -= gpu_constants.k_s * diff_of_products_simple(n1.z, temp_vec.x, n1.x, temp_vec.z);
            torque.z -= gpu_constants.k_s * diff_of_products_simple(n1.x, temp_vec.y, n1.y, temp_vec.x);
#else
            torque.x -= gpu_constants.k_s * diff_of_products(n1.y, temp_vec.z, n1.z, temp_vec.y);
            torque.y -= gpu_constants.k_s * diff_of_products(n1.z, temp_vec.x, n1.x, temp_vec.z);
            torque.z -= gpu_constants.k_s * diff_of_products(n1.x, temp_vec.y, n1.y, temp_vec.x);
#endif


#endif // GPU_SLIDING




            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // rolling
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef GPU_ROLLING
            // Note: k_r is actually k_r * reduced_radius^2
#ifdef CPU_EQUIVALENT
            torque.x -= gpu_constants.k_r * diff_of_products_simple(n1.y, n2.z, n1.z, n2.y);
            torque.y -= gpu_constants.k_r * diff_of_products_simple(n1.z, n2.x, n1.x, n2.z);
            torque.z -= gpu_constants.k_r * diff_of_products_simple(n1.x, n2.y, n1.y, n2.x);
#else
            torque.x -= gpu_constants.k_r * diff_of_products(n1.y, n2.z, n1.z, n2.y);
            torque.y -= gpu_constants.k_r * diff_of_products(n1.z, n2.x, n1.x, n2.z);
            torque.z -= gpu_constants.k_r * diff_of_products(n1.x, n2.y, n1.y, n2.x);
#endif  // CPU_CPU_EQUIVALENT
#endif  // GPU_ROLLING

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// twisting
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef GPU_TWISTING
            // Note n1.w = twisting_displacement
			torque -= (gpu_constants.k_t * n1.w) * make_double3(n_c);

#endif

		}
	}
/*
    if(isnan(dot(force, torque))) // TODO: remove this   debug purpose only
        *isNaN = 1;
*/

	// write back total force/torque
	store_double3(forces_new, index, number_of_particles, force);
	store_double3(torques_new, index, number_of_particles, torque);

	return;
}

__global__
void updateBoxInteractionKernel(
	double* RESTRICT forces_new,
	const double* RESTRICT positions,
	const double3 lower_pos,
	const double3 upper_pos,
	double* top_wall_force,
	double* bot_wall_force,
	const unsigned int number_of_particles
        )
{
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size

	double3 pos = make_double3(positions, index, number_of_particles);

	// determine distance to wall
	double3 wall_force = make_double3(0.0);

	double dist = pos.x - lower_pos.x;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		wall_force.x += interaction_strength;

	}

	dist = pos.y - lower_pos.y;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		wall_force.y += interaction_strength;
		bot_wall_force[index] += interaction_strength;

	}

	dist = pos.z - lower_pos.z;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		wall_force.z += interaction_strength;

	}

	dist = upper_pos.x - pos.x;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		wall_force.x -= interaction_strength;

	}

	dist = upper_pos.y - pos.y;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		wall_force.y -= interaction_strength;
        top_wall_force[index] += interaction_strength;

	}

	dist = upper_pos.z - pos.z;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		wall_force.z -= interaction_strength;

	}

	if(wall_force.x != 0)
		forces_new[X_COORD_GPU(index, number_of_particles)] += wall_force.x;

	if(wall_force.y != 0)
		forces_new[Y_COORD_GPU(index, number_of_particles)] += wall_force.y;

	if(wall_force.z != 0)
		forces_new[Z_COORD_GPU(index, number_of_particles)] += wall_force.z;


}



__global__
void updateOpenBoxInteractionKernel(
	double* RESTRICT forces_new,
	const double* RESTRICT positions,
	const double3 lower_pos,
	const double3 upper_pos,
	const unsigned int number_of_particles
        )
{
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size

	double3 pos = make_double3(positions, index, number_of_particles);

	// bottom wall
	double dist = pos.y - lower_pos.y;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		forces_new[Y_COORD_GPU(index, number_of_particles)] += interaction_strength;
	}

	// side walls
	if(pos.y < upper_pos.y)
	{
		dist = pos.x - lower_pos.x;
		if(fabs(dist) < gpu_constants.particle_radius)
		{
			double delta = gpu_constants.particle_radius - dist;
			if(dist < 0)
				delta += 2.0 * dist;

			double contact_radius = getContactRadius(delta);
			double interaction_strength = getNormalForce(contact_radius);
			forces_new[X_COORD_GPU(index, number_of_particles)] += interaction_strength;
		}

		dist = pos.z - lower_pos.z;
		if(fabs(dist) < gpu_constants.particle_radius)
		{
			double delta = gpu_constants.particle_radius - dist;
			if(dist < 0)
				delta += 2.0 * dist;

			double contact_radius = getContactRadius(delta);
			double interaction_strength = getNormalForce(contact_radius);
			forces_new[Z_COORD_GPU(index, number_of_particles)] += interaction_strength;
		}

		dist = upper_pos.x - pos.x;
		if(fabs(dist) < gpu_constants.particle_radius)
		{
			double delta = gpu_constants.particle_radius - dist;
			if(dist < 0)
				delta += 2.0 * dist;

			double contact_radius = getContactRadius(delta);
			double interaction_strength = getNormalForce(contact_radius);
			forces_new[X_COORD_GPU(index, number_of_particles)] -= interaction_strength;
		}

		dist = upper_pos.z - pos.z;
		if(fabs(dist) < gpu_constants.particle_radius)
		{
			double delta = gpu_constants.particle_radius - dist;
			if(dist < 0)
				delta += 2.0 * dist;

			double contact_radius = getContactRadius(delta);
			double interaction_strength = getNormalForce(contact_radius);
			forces_new[Z_COORD_GPU(index, number_of_particles)] -= interaction_strength;
		}
	}
}

__global__ void updateBottomWallInteractionKernel(
        double* forces_new,
        const double* RESTRICT positions,
        double3 lower_pos,
        unsigned int number_of_particles
        )
{
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size

	double3 pos = make_double3(positions, index, number_of_particles);

	// bottom wall
	double dist = pos.y - lower_pos.y;
	if(dist < gpu_constants.particle_radius)
	{
		double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
		double interaction_strength = getNormalForce(contact_radius);
		forces_new[Y_COORD_GPU(index, number_of_particles)] += interaction_strength;
	}
}



__global__ void updateEnclosingSphereWallInteractionKernel(
        double* forces_new,
        const double* RESTRICT positions,
        const double sphere_radius,
        unsigned int number_of_particles
        )
{
    unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

    if(index >= number_of_particles)
        return;          // handle case when no. of particles not multiple of block size

    double3 pos = make_double3(positions, index, number_of_particles);

    // enclosing sphere
    double dist = sqrt(dot(pos, pos));
    pos /= dist;

    dist = sphere_radius - dist;

    if(dist < gpu_constants.particle_radius)
    {
        double contact_radius = getContactRadius(gpu_constants.particle_radius - dist);
        double interaction_strength = -getNormalForce(contact_radius);

        add_double3(forces_new, index, number_of_particles, interaction_strength*pos);
    }
}





#ifdef GPU_TRACK_PARTICLE_ORIENTATION
__global__
void updateParticleOrientationKernel(
	double4* RESTRICT orientations,
	const double* RESTRICT angular_velocities,
	const double* RESTRICT torques,
	const double timestep,
	const unsigned int number_of_particles
        )
{
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;          // handle case when no. of particles not multiple of block size

	double4 orientation = orientations[index];

	double3 omega = make_double3(angular_velocities, index, number_of_particles);

	double3 omega_dot = make_double3(torques, index, number_of_particles);


	omega_dot *= gpu_constants.moment_of_inertia_inv;

	
	double omega_squared = dot(omega, omega);

	double4 e_dot;
	e_dot.w = - 0.5 * dot(orientation, omega);
	e_dot.x = 0.5 * (orientation.w * omega.x - orientation.y * omega.z + orientation.z * omega.y);
	e_dot.y = 0.5 * (orientation.w * omega.y - orientation.z * omega.x + orientation.x * omega.z);
	e_dot.z = 0.5 * (orientation.w * omega.z - orientation.x * omega.y + orientation.y * omega.x);

	double temp = 0.5 * e_dot.w;

	double4 e_ddot;
	e_ddot.w = - 0.25 * (orientation.w * omega_squared + 2.0 * dot(orientation, omega_dot));
	e_ddot.x = temp * omega.x + 0.5 * (orientation.w * omega_dot.x - orientation.y * omega_dot.z + orientation.z * omega_dot.y);
	e_ddot.y = temp * omega.y + 0.5 * (orientation.w * omega_dot.y - orientation.z * omega_dot.x + orientation.x * omega_dot.z);
	e_ddot.z = temp * omega.z + 0.5 * (orientation.w * omega_dot.z - orientation.x * omega_dot.y + orientation.y * omega_dot.x);

	orientation += timestep * (e_dot + 0.5 * timestep * e_ddot);


	// make sure that e0^2 + .. + e3^2 = 1
	double norm_inv = rsqrt(orientation.w*orientation.w + orientation.x*orientation.x + orientation.y*orientation.y + orientation.z*orientation.z);

	orientation *= norm_inv;


	orientations[index] = orientation;
}
#endif

#ifdef GPU_TRACK_DISSIPATED_ENERGY
__global__
void calculateKineticEnergyKernel(
        double* particle_kinetic_energy,
        const double* RESTRICT velocities,
        unsigned int number_of_particles
        )
{
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= number_of_particles)
		return;

    double3 velocity = make_double3(velocities, index, number_of_particles);

    double E = 0.5*gpu_constants.mass * dot(velocity, velocity);

    particle_kinetic_energy[index] = E;
}


__global__
void calculateRollingEnergyKernel(
        double* particle_rolling_energy,
        const double* RESTRICT angular_velocities,
        unsigned int number_of_particles
        )
{
    unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

    if(index >= number_of_particles)
        return;

    double3 angular_velocity = make_double3(angular_velocities, index, number_of_particles);

    double E = 0.5*gpu_constants.moment_of_inertia * dot(angular_velocity, angular_velocity);

    particle_rolling_energy[index] = E;
}




__global__
void calculateContactEnergyKernel(
        double* RESTRICT V_normal,
        double* RESTRICT V_rolling,
        double* RESTRICT V_sliding,
        double* RESTRICT V_twisting,
        int2* RESTRICT contact_particle_ids,
        double4* RESTRICT contact_normals,
        double4* RESTRICT n_1,
        double4* RESTRICT n_2,
        const unsigned int number_of_particles
        )
{
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= 6 * number_of_particles)
		return;

	// determine particle ids
	int2 particle_ids = contact_particle_ids[index];

	if(particle_ids.x == NO_PARTICLE)
		return;

	//////////////////////////////////////////////////////////////////////////////
	// determine current contact pointers
	//////////////////////////////////////////////////////////////////////////////

	double4 n_c = contact_normals[index];
	double4 n1 = n_1[index];
	double4 n2 = n_2[index];

	//////////////////////////////////////////////////////////////////////////////
	// calculate potential energy
	//////////////////////////////////////////////////////////////////////////////

	// detrmine normal potential
	double temp = getContactRadius(2.0 * gpu_constants.particle_radius - n_c.w);
	temp /= gpu_constants.eq_contact_radius;
	
    V_normal[index] = 7.2684823713285586 * gpu_constants.F_c * gpu_constants.delta_c * (0.8 * temp*temp*temp*temp*temp - 1.3333333333333333 * temp*temp*temp*sqrt(temp) + 0.3333333333333333 * temp*temp);
	//double E_pot = 7.268482371328559 * gpu_constants.F_c * gpu_constants.delta_c * (0.8 * temp*temp*temp*temp*temp - 4.0 / 3.0 * pow(temp, 3.5) + 1.0 / 3.0 * temp*temp); // constants replaced with value
	//double E_pot = 4.0 * pow(6.0, 1.0 / 3.0) * gpu_constants.F_c * gpu_constants.delta_c * (0.8 * temp*temp*temp*temp*temp - 4.0 / 3.0 * pow(temp, 3.5) + 1.0 / 3.0 * temp*temp); // original from alex

	double3 displacement;

#ifdef GPU_ROLLING
	displacement = make_double3(n1 + n2);
    V_rolling[index] = 0.5 * gpu_constants.k_r * dot(displacement, displacement);
#endif // GPU_ROLLING

#ifdef GPU_SLIDING
    displacement = (make_double3(n1 - n2 + 2.0 * n_c));

	temp = dot(displacement, n_c);

	displacement -= temp * make_double3(n_c);

    V_sliding[index] = 0.5 * gpu_constants.k_s * dot(displacement, displacement);

#endif // GPU_SLIDING

#ifdef GPU_TWISTING
    V_twisting[index] = 0.5 * gpu_constants.k_t * n1.w * n1.w;
#endif // GPU_TWISTING

}
#endif // GPU_TRACK_DISSIPATED_ENERGY

#endif // CUDA_KERNELS_H
