 /***************************************************************************
 *   Simulation of particle collisions/agglomeration						*
 *   Copyright (C) 2009 Alexander Seizinger									*
 *	 alexander.seizinger[...]gmx.net										*
 *																			*
 *   This program is free software; you can redistribute it and/or modify	*
 *   it under the terms of the GNU General Public License as published by	*
 *   the Free Software Foundation; either version 2 of the License, or		*
 *   (at your option) any later version.									*
 *																			*
 *   This program is distributed in the hope that it will be useful,		*
 *   but WITHOUT ANYs WARRANTY; without even the implied warranty of		*
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the			*
 *   GNU General Public License for more details.							*
 *																			*
 *   You should have received a copy of the GNU General Public License		*
 *   along with this program; if not, write to the							*
 *   Free Software Foundation, Inc.											*
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.				*
 ***************************************************************************/

#include "SimulationCuda.h"
#include "CudaDefines.cuh"
#include "CudaWrapper.cuh"
#include "SimulationLib.h"
#include <vector_functions.h>
#include <string.h>
#include <helper_cuda.h>

#include "CudaWallKernelsWrapper.cuh"

#if defined(GPU_USE_CONSTANT_DAMPING) && defined(GPU_USE_VISCOELASTIC_DAMPING)
#error Both constant and viscoelastic damping on GPU enabled - select only one!
#endif

extern double particle_radius;
extern double density;
extern double surface_energy;
extern double nu;
extern double young_mod;
extern double crit_rolling_displacement;
extern double osc_damping_factor;
extern double T_vis;
extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;
extern double crit_sliding_displacement_modifier;
extern double crit_wall_sliding_displacement_modifier;
extern double delta_c;
extern double moment_of_inertia;
extern double mass;
extern double ENERGY_UNIT;
extern double contact_making_dist;

#if defined (GPU_TRACK_PARTICLE_ORIENTATION) && !defined(TRACK_PARTICLE_ORIENTATION)
#error GPU_TRACK_PARTICLE_ORIENTATION defined but TRACK_PARTICLE_ORIENTATION missing!
#endif


void SimulationCuda::copyParticleArrayFromGPU(double* host_ptr, double* device_ptr, int number_of_particles, double* host_buffer)
{

    size_t memsize = number_of_particles * 3 * sizeof(double);

    checkCudaErrors(cudaMemcpy(host_buffer, device_ptr, memsize, cudaMemcpyDeviceToHost));



    for(int p = 0; p < number_of_particles; ++p)
    {
        host_ptr[X_COORD(p)] = host_buffer[X_COORD_GPU(p, number_of_particles)];
        host_ptr[Y_COORD(p)] = host_buffer[Y_COORD_GPU(p, number_of_particles)];
        host_ptr[Z_COORD(p)] = host_buffer[Z_COORD_GPU(p, number_of_particles)];
    }

}

void SimulationCuda::copyParticleArrayToGPU(double* device_ptr, double* host_ptr, int number_of_particles, double* host_buffer)
{

    size_t memsize = number_of_particles * 3 * sizeof(double);


    for(int p = 0; p < number_of_particles; ++p)
    {
        host_buffer[X_COORD_GPU(p, number_of_particles)] = host_ptr[X_COORD(p)];
        host_buffer[Y_COORD_GPU(p, number_of_particles)] = host_ptr[Y_COORD(p)];
        host_buffer[Z_COORD_GPU(p, number_of_particles)] = host_ptr[Z_COORD(p)];
    }
    checkCudaErrors(cudaMemcpy(device_ptr, host_buffer, memsize, cudaMemcpyHostToDevice));


}



SimulationCuda::SimulationCuda(void) : Simulation()
{
	gpu_positions_old = NULL;
    gpu_positions_new = NULL;
    gpu_positions_sorted = NULL;

	gpu_velocities = NULL;
	gpu_angular_velocities = NULL;
	gpu_forces_old = NULL;
	gpu_forces_new = NULL;
	gpu_torques_old = NULL;
	gpu_torques_new = NULL;


#ifdef GPU_TRACK_PARTICLE_ORIENTATION
	gpu_particle_orientation = NULL;
#endif

	gpu_grid_particle_index = NULL;
	gpu_grid_particle_cell = NULL;
	gpu_cell_start = NULL;
	gpu_cell_end = NULL;

	gpu_particle_number_of_contacts = NULL;
	gpu_particle_particle_ids = NULL;
	gpu_particle_contact_ids = NULL;
	gpu_update_local_contact_list = NULL;

    gpu_contact_type = NULL;
	gpu_contact_particle_ids = NULL;
	gpu_contact_normals_new = NULL;
	gpu_contact_normals_old = NULL;
	gpu_contact_rot1 = NULL;
	gpu_contact_rot2 = NULL;
	gpu_contact_n1_initial = NULL;
	gpu_contact_n2_initial = NULL;
	gpu_contact_n1 = NULL;
	gpu_contact_n1 = NULL;


	gpu_free_contacts_list = NULL;
	gpu_next_free_contact = NULL;
	gpu_new_contact_particle_ids = NULL;
	gpu_number_of_new_contacts = NULL;


#ifdef GPU_TRACK_DISSIPATED_ENERGY

    gpu_dissipated_damping_energy = NULL;
    gpu_dissipated_rolling_energy = NULL;
    gpu_dissipated_contact_energy = NULL;
    gpu_dissipated_sliding_energy = NULL;
    gpu_dissipated_twisting_energy = NULL;

    gpu_kinetic_energy = NULL;
    gpu_rotation_energy = NULL;

    gpu_normal_energy = NULL;
    gpu_rolling_energy = NULL;
    gpu_sliding_energy = NULL;
    gpu_twisting_energy = NULL;

    E_tot = 0.0;
    E_kin = 0.0;
    E_rot = 0.0;
    V_tot = 0.0;
    V_normal = 0.0;
    V_rolling = 0.0;
    V_sliding = 0.0;
    V_twisting = 0.0;

#endif

    m_sphere_radius = 0.0;
    m_sphere_compaction_speed = 0.0;

    E_kin_last = 0.0;
    V_tot_last = 0.0;

	damping_factor = 0;
	check_filling_factor_counter = 0;

	num_grid_cells = 0;

	cuda_initialized = false;
	use_gpu = false;
}

SimulationCuda::~SimulationCuda(void)
{
	cleanUpCuda();
}

void SimulationCuda::cleanUpCuda()
{
	if(gpu_positions_old)
	{
        cubPlan.DestroyPlan();

        checkCudaErrors(cudaFree(gpu_positions_old));
        checkCudaErrors(cudaFree(gpu_positions_new));
        checkCudaErrors(cudaFree(gpu_positions_sorted));

        checkCudaErrors(cudaFree(gpu_velocities));
        checkCudaErrors(cudaFree(gpu_angular_velocities));
        checkCudaErrors(cudaFree(gpu_forces_old));
        checkCudaErrors(cudaFree(gpu_forces_new));
        checkCudaErrors(cudaFree(gpu_torques_old));
        checkCudaErrors(cudaFree(gpu_torques_new));

        checkCudaErrors(cudaFree(gpu_grid_particle_index));
        checkCudaErrors(cudaFree(gpu_grid_particle_cell));
        checkCudaErrors(cudaFree(gpu_cell_start));
        checkCudaErrors(cudaFree(gpu_cell_end));


        gpu_positions_old = NULL;
        gpu_positions_new = NULL;
        gpu_positions_sorted = NULL;

		gpu_velocities = NULL;
		gpu_angular_velocities = NULL;
		gpu_forces_old = NULL;
		gpu_forces_new = NULL;
		gpu_torques_old = NULL;
		gpu_torques_new = NULL;


#ifdef GPU_TRACK_PARTICLE_ORIENTATION
		cudaFree(gpu_particle_orientation);
		gpu_particle_orientation = NULL;
#endif

#ifdef GPU_TRACK_DISSIPATED_ENERGY

        cudaFree(gpu_dissipated_damping_energy);
        cudaFree(gpu_dissipated_rolling_energy);
        cudaFree(gpu_dissipated_contact_energy);
        cudaFree(gpu_dissipated_sliding_energy);
        cudaFree(gpu_dissipated_twisting_energy);

        cudaFree(gpu_kinetic_energy);
        cudaFree(gpu_rotation_energy);

        cudaFree(gpu_normal_energy);
        cudaFree(gpu_rolling_energy);
        cudaFree(gpu_sliding_energy);
        cudaFree(gpu_twisting_energy);
		
        gpu_kinetic_energy = NULL;
        gpu_rotation_energy = NULL;

        gpu_normal_energy = NULL;
        gpu_rolling_energy = NULL;
        gpu_sliding_energy = NULL;
        gpu_twisting_energy = NULL;
#endif

		gpu_grid_particle_index = NULL;
		gpu_grid_particle_cell = NULL;
		gpu_cell_start = NULL;
		gpu_cell_end = NULL;


        cudaDeviceSynchronize();
        getLastCudaError("Cleaning up CUDA particles failed!");


	}

	if(gpu_contact_particle_ids)
	{

        if(this->box)
        {
            wallInteraction.deleteParticleWallInteraction();
        }

        checkCudaErrors(cudaFree(gpu_particle_number_of_contacts));
        checkCudaErrors(cudaFree(gpu_particle_particle_ids));
        checkCudaErrors(cudaFree(gpu_particle_contact_ids));
        checkCudaErrors(cudaFree(gpu_update_local_contact_list));

        checkCudaErrors(cudaFree(gpu_contact_type));
        checkCudaErrors(cudaFree(gpu_contact_particle_ids));
        checkCudaErrors(cudaFree(gpu_contact_normals_new));
        checkCudaErrors(cudaFree(gpu_contact_normals_old));
        checkCudaErrors(cudaFree(gpu_contact_rot1));
        checkCudaErrors(cudaFree(gpu_contact_rot2));
        checkCudaErrors(cudaFree(gpu_contact_n1_initial));
        checkCudaErrors(cudaFree(gpu_contact_n2_initial));
        checkCudaErrors(cudaFree(gpu_contact_n1));
        checkCudaErrors(cudaFree(gpu_contact_n2));


        checkCudaErrors(cudaFree(gpu_free_contacts_list));
        checkCudaErrors(cudaFree(gpu_next_free_contact));
        checkCudaErrors(cudaFree(gpu_new_contact_particle_ids));
        checkCudaErrors(cudaFree(gpu_number_of_new_contacts));

		gpu_particle_number_of_contacts = NULL;
        gpu_contact_type = NULL;
        gpu_particle_particle_ids = NULL;
        gpu_particle_contact_ids = NULL;
		gpu_update_local_contact_list = NULL;

        gpu_contact_type = NULL;
		gpu_contact_particle_ids = NULL;
		gpu_contact_normals_new = NULL;
		gpu_contact_normals_old = NULL;
		gpu_contact_rot1 = NULL;
		gpu_contact_rot2 = NULL;
		gpu_contact_n1_initial = NULL;
		gpu_contact_n2_initial = NULL;
		gpu_contact_n1 = NULL;
		gpu_contact_n2 = NULL;


		gpu_free_contacts_list = NULL;
		gpu_next_free_contact = NULL;
		gpu_new_contact_particle_ids = NULL;
		gpu_number_of_new_contacts = NULL;
	}


    cudaDeviceSynchronize();
    getLastCudaError("Cleaning up CUDA contacts failed!");

	cuda_initialized = false;
	use_gpu = false;
}

void SimulationCuda::setDampingFactor(double damping_factor)
{
	this->damping_factor = damping_factor;
}

double SimulationCuda::getGPUFillingFactor()
{

    if(this->use_gpu)
    {

        if(this->number_of_walls == 6)
        {
            this->wallInteraction.get_box_pos_from_gpu(*this);
            this->copyPositionsFromGPU(false);

            return this->getBoxFillingFactor();

        }
        else if(this->number_of_walls == 2)
        {

            this->wallInteraction.get_box_pos_from_gpu(*this);

            double cross_section;
            this->copyPositionsFromGPU(false);
            SimLib::getCrossSectionNoRotation(*this, 0.1*particle_radius, cross_section);

            return (this->number_of_particles * 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius) / (cross_section * this->box->height);

        }
    }
    else
    {
        return this->getBoxFillingFactor();
    }
}



double SimulationCuda::getGPUFillingFactorNoSW()
{

    copyPositionsFromGPU();
    double cross_section = 0.0;
    SimLib::getCrossSectionNoRotation(*this, 0.2*particle_radius, cross_section);


    if(cross_section > 0 && box)
        return (number_of_particles * 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius) / (cross_section * box->height);
    else
        return 0;

}

int SimulationCuda::getGPUNumContacts()
{
    int gpu_num_contacts = 0;
    checkCudaErrors(cudaMemcpy(&gpu_num_contacts, gpu_next_free_contact, sizeof(int), cudaMemcpyDeviceToHost));

    return gpu_num_contacts;

}


bool SimulationCuda::sufficientMemory(size_t required_memory)
{
	size_t free;
	size_t total;

	cudaError_t cuda_error = cudaMemGetInfo(&free, &total);

	if(cuda_error != cudaSuccess || free < required_memory)
		return false;
	else
		return true;
}

ErrorCode SimulationCuda::initCuda(int GPU_id)
{
	check_filling_factor_counter = 0;

	/////////////////////////////////////////////////////////////////////////////////////////
	// check walls
	/////////////////////////////////////////////////////////////////////////////////////////

    if(number_of_walls != 0 && number_of_walls != 1 && number_of_walls != 5 && number_of_walls != 6 && number_of_walls != 2)
		return EC_CUDA_INVALID_BOUNDARIES;

	/////////////////////////////////////////////////////////////////////////////////////////
	// select proper device
	/////////////////////////////////////////////////////////////////////////////////////////

	if(GPU_id < 0)
	{
		int devices;
		cudaGetDeviceCount( &devices );

        int mp_count = 0;

		for(int i = 0; i < devices; i++) 
		{
			cudaGetDeviceProperties(&device_properties, i );

			// check if comput capability is sufficient (when double support is required)
#ifdef USE_DOUBLE_PRECISION
			if(device_properties.major >= 2 || device_properties.minor >= 3)
#else
			if(true)
#endif
			{
				if(device_properties.multiProcessorCount > mp_count)
				{
					mp_count = device_properties.multiProcessorCount;
					GPU_id = i;
				}
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// set device
	/////////////////////////////////////////////////////////////////////////////////////////

	if(GPU_id < 0)
		return EC_NO_CUDA_DEVICE;
	else
	{
		cudaError_t cuda_error = cudaSetDevice(GPU_id);

		if(cuda_error != cudaSuccess)
			return EC_CUDA_DEVICE_UNAVAILABLE;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// allocate memory and copy data
	/////////////////////////////////////////////////////////////////////////////////////////
    num_grid_cells = 128;

    initGPUConstants(
                this->normal_interaction,
                this->num_grid_cells
                );



    const size_t memsize = number_of_particles * sizeof(double3);
    const int max_contacts = 6 * number_of_particles;

	size_t required_memory = 0;
    required_memory += 9 * memsize;	// arrays for particle data
    required_memory += 2 * number_of_particles * sizeof(unsigned int) + 2 * num_grid_cells*num_grid_cells*num_grid_cells * sizeof(unsigned int); // grid

    required_memory += 2 * number_of_particles * sizeof(unsigned int); // local contact list
    required_memory += 2 * number_of_particles * MAX_CONTACTS * sizeof(int); // local contact list

    required_memory += max_contacts * (6 * sizeof(double4) + 2 * sizeof(double3) + 2 * sizeof(int2) + sizeof(int)); // global contact list
	required_memory += max_contacts * (sizeof(int) + sizeof(int2)) + 2 * sizeof(int);	// new contacts
	required_memory += number_of_particles * sizeof(int); // deleted contacts 
	required_memory += 2 * number_of_particles * sizeof(unsigned int); // for cub sorting algorithm 


    //printf("memory required = %d mb\n", required_memory/(1024*1024));

    if(number_of_walls > 0)
        required_memory += 2 * number_of_particles * sizeof(double); // for summing up the force the monomers exert on the top & bottom wall of the compression box

#ifdef GPU_TRACK_PARTICLE_ORIENTATION
	required_memory += number_of_particles * sizeof(double4);
#endif

#ifdef GPU_TRACK_DISSIPATED_ENERGY
    required_memory += 4*number_of_particles * sizeof(double);
    required_memory += 12*number_of_contacts * sizeof(double);
#endif

	if(!sufficientMemory(required_memory))
	{
		size_t free;
		size_t total;

		cudaError_t cuda_error = cudaMemGetInfo(&free, &total);

        printf("total / free / required_memory: %zu / %zu / %zu\n", total/(1024*1024), free/(1024*1024), required_memory/(1024*1024));
	
		return EC_INSUFFICIENT_GPU_MEMORY;
	}

	// particle data
	cudaMalloc((void**)&gpu_positions_old, memsize);
    cudaMalloc((void**)&gpu_positions_new, memsize);
    cudaMalloc((void**)&gpu_positions_sorted, memsize);

	cudaMalloc((void**)&gpu_velocities, memsize);
	cudaMalloc((void**)&gpu_angular_velocities, memsize);
	cudaMalloc((void**)&gpu_forces_old, memsize);
	cudaMalloc((void**)&gpu_forces_new, memsize);
	cudaMalloc((void**)&gpu_torques_old, memsize);
	cudaMalloc((void**)&gpu_torques_new, memsize);

#ifdef GPU_TRACK_PARTICLE_ORIENTATION
	cudaMalloc((void**)&gpu_particle_orientation, number_of_particles * sizeof(double4));
#endif

#ifdef GPU_TRACK_DISSIPATED_ENERGY

	checkCudaErrors(cudaMalloc((void**)&gpu_dissipated_damping_energy, number_of_particles * sizeof(double)));

    checkCudaErrors(cudaMalloc((void**)&gpu_dissipated_rolling_energy, max_contacts * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&gpu_dissipated_contact_energy, max_contacts * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&gpu_dissipated_sliding_energy, max_contacts * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&gpu_dissipated_twisting_energy, max_contacts * sizeof(double)));

	checkCudaErrors(cudaMemset(gpu_dissipated_damping_energy, 0, number_of_particles * sizeof(double)));

    checkCudaErrors(cudaMemset(gpu_dissipated_rolling_energy, 0, max_contacts * sizeof(double)));
    checkCudaErrors(cudaMemset(gpu_dissipated_contact_energy, 0, max_contacts * sizeof(double)));
    checkCudaErrors(cudaMemset(gpu_dissipated_sliding_energy, 0, max_contacts * sizeof(double)));
    checkCudaErrors(cudaMemset(gpu_dissipated_twisting_energy, 0, max_contacts * sizeof(double)));

	checkCudaErrors(cudaMalloc((void**)&gpu_kinetic_energy, number_of_particles * sizeof(double)));
	checkCudaErrors(cudaMalloc((void**)&gpu_rotation_energy, number_of_particles * sizeof(double)));

	checkCudaErrors(cudaMalloc((void**)&gpu_normal_energy, max_contacts * sizeof(double)));
	checkCudaErrors(cudaMalloc((void**)&gpu_rolling_energy, max_contacts * sizeof(double)));
	checkCudaErrors(cudaMalloc((void**)&gpu_sliding_energy, max_contacts * sizeof(double)));
	checkCudaErrors(cudaMalloc((void**)&gpu_twisting_energy, max_contacts * sizeof(double)));

#endif


    // grid
	checkCudaErrors(cudaMalloc((void**)&gpu_grid_particle_index, number_of_particles * sizeof(unsigned int)));
	checkCudaErrors(cudaMalloc((void**)&gpu_grid_particle_cell, number_of_particles * sizeof(unsigned int)));
	checkCudaErrors(cudaMalloc((void**)&gpu_cell_start, num_grid_cells*num_grid_cells*num_grid_cells * sizeof(unsigned int)));
	checkCudaErrors(cudaMalloc((void**)&gpu_cell_end, num_grid_cells*num_grid_cells*num_grid_cells * sizeof(unsigned int)));


    // use buffer to convert from cpu to gpu addressing scheme
	double* buffer = new double[3*number_of_particles];

    copyParticleArrayToGPU(gpu_positions_old, pos_old, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_velocities, vel, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_angular_velocities, vel_angular, number_of_particles, buffer);



    copyParticleArrayToGPU(gpu_positions_new, pos_new, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_velocities, vel, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_angular_velocities, vel_angular, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_forces_old, force_old, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_forces_new, force_new, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_torques_old, torque_old, number_of_particles, buffer);

    copyParticleArrayToGPU(gpu_torques_new, torque_new, number_of_particles, buffer);


	delete [] buffer;


#ifdef GPU_TRACK_PARTICLE_ORIENTATION
	double4 *orientation_buffer = new double4[number_of_particles];

	for(int p = 0; p < number_of_particles; ++p)
	{
		orientation_buffer[p].w = orientation[4*p];
		orientation_buffer[p].x = orientation[4*p+1];
		orientation_buffer[p].y = orientation[4*p+2];
		orientation_buffer[p].z = orientation[4*p+3];
	}

	checkCudaErrors(cudaMemcpy(gpu_particle_orientation, orientation_buffer, number_of_particles * sizeof(double4), cudaMemcpyHostToDevice));

	delete [] orientation_buffer;
#endif

	/////////////////////////////////////////////////////////////////////////////////////////
	// set up contact list
	/////////////////////////////////////////////////////////////////////////////////////////

    checkCudaErrors(cudaMalloc((void**)&gpu_particle_number_of_contacts, number_of_particles * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&gpu_particle_contact_ids, number_of_particles * MAX_CONTACTS * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&gpu_particle_particle_ids, number_of_particles * MAX_CONTACTS * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&gpu_update_local_contact_list, number_of_particles * sizeof(int)));


    checkCudaErrors(cudaMalloc((void**)&gpu_contact_type, max_contacts * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_particle_ids, max_contacts * sizeof(int2)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_n1, max_contacts * sizeof(double4)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_n2, max_contacts * sizeof(double4)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_normals_old, max_contacts * sizeof(double4)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_normals_new, max_contacts * sizeof(double4)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_n1_initial, max_contacts * sizeof(double3)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_n2_initial, max_contacts * sizeof(double3)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_rot1, max_contacts * sizeof(double4)));
    checkCudaErrors(cudaMalloc((void**)&gpu_contact_rot2, max_contacts * sizeof(double4)));

    checkCudaErrors(cudaMalloc((void**)&gpu_next_free_contact, sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&gpu_number_of_new_contacts, sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&gpu_new_contact_particle_ids, max_contacts * sizeof(int2)));
    checkCudaErrors(cudaMalloc((void**)&gpu_free_contacts_list, max_contacts * sizeof(int)));

	checkCudaErrors(cudaMemset(gpu_contact_normals_new, 0, max_contacts * sizeof(double4) ));
	checkCudaErrors(cudaMemset(gpu_number_of_new_contacts, 0, sizeof(int) ));
	checkCudaErrors(cudaMemset(gpu_update_local_contact_list, 0, number_of_particles * sizeof(int)));
	checkCudaErrors(cudaMemset(gpu_new_contact_particle_ids, NO_PARTICLE, max_contacts * sizeof(int2)));

	// allocate buffers to be copied to gpu later
	int *contact_number_buffer = new int[number_of_particles];
	int *contact_id_buffer = new int[number_of_particles * MAX_CONTACTS];
	int *particle_id_buffer = new int[number_of_particles * MAX_CONTACTS];

	int next_free_contact_id = 0;
	int *free_contacts_list_buffer = new int[max_contacts];
	int2 *contact_particle_id_buffer = new int2[max_contacts];
	double4 *n1_buffer = new double4[max_contacts];
	double4 *n2_buffer = new double4[max_contacts];
	double4 *contact_normals_old_buffer = new double4[max_contacts];
    double4 *contact_normals_new_buffer = new double4[max_contacts];
    double3 *contact_n1_initial_buffer = new double3[max_contacts];
    double3 *contact_n2_initial_buffer = new double3[max_contacts];
	double4 *contact_rot1_buffer = new double4[max_contacts];
	double4 *contact_rot2_buffer = new double4[max_contacts];


	// init empty contact list
	for(int p = 0; p < number_of_particles; ++p)
	{
		contact_number_buffer[p] = 0;

		for(int c = 0; c < MAX_CONTACTS; ++c)
		{
			contact_id_buffer[CONTACT_ID(c, p, number_of_particles)] = -1;
			particle_id_buffer[CONTACT_ID(c, p, number_of_particles)] = NO_PARTICLE;
		}
	}

	for(int c = 0; c < max_contacts; ++c)
	{
		free_contacts_list_buffer[c] = c;
		contact_particle_id_buffer[c] = make_int2(NO_PARTICLE, NO_PARTICLE);
	}


	// copy existing contact list
	ContactListEntry *cl_entry = NULL;
	vec3 n1, n2;

	for(int p = 0; p < number_of_particles; ++p)
	{

		cl_entry = contact_list[p];

		while(cl_entry)
		{


            if(cl_entry->id < 0) // ignore contact with wall
            {
                cl_entry = cl_entry->next;
                continue;
            }

            if(cl_entry->id < 0)
                printf("wall bond:  %d  %d\n", cl_entry->contact->id1, cl_entry->contact->id2);


			// determine contact normal
			double4 n_c;
			n_c.x = pos_old[X_COORD(p)] - pos_old[X_COORD(cl_entry->id)];
			n_c.y = pos_old[Y_COORD(p)] - pos_old[Y_COORD(cl_entry->id)];
			n_c.z = pos_old[Z_COORD(p)] - pos_old[Z_COORD(cl_entry->id)];
			n_c.w = sqrt(n_c.x*n_c.x + n_c.y*n_c.y + n_c.z*n_c.z);

            double tmp = 1.0 / n_c.w;

            n_c.x *= tmp;
            n_c.y *= tmp;
            n_c.z *= tmp;


			// add contact to contact list
			contact_particle_id_buffer[next_free_contact_id] = make_int2(p, cl_entry->id);

			cl_entry->contact->getCurrentN1(&n1);
			cl_entry->contact->getCurrentN2(&n2);

			n1_buffer[next_free_contact_id] = make_double4(n1[0], n1[1], n1[2], cl_entry->contact->twisting_displacement);
			n2_buffer[next_free_contact_id] = make_double4(n2[0], n2[1], n2[2], cl_entry->contact->compression_length);


            contact_normals_new_buffer[next_free_contact_id] = n_c;

            n_c.w = 2.0 * particle_radius - cl_entry->contact->compression_length; // to be consistent with the CPU version
            contact_normals_old_buffer[next_free_contact_id] = n_c;


			contact_rot1_buffer[next_free_contact_id] = make_double4(cl_entry->contact->rot1.e1, cl_entry->contact->rot1.e2, cl_entry->contact->rot1.e3, cl_entry->contact->rot1.e0);
			contact_rot2_buffer[next_free_contact_id] = make_double4(cl_entry->contact->rot2.e1, cl_entry->contact->rot2.e2, cl_entry->contact->rot2.e3, cl_entry->contact->rot2.e0);

            contact_n1_initial_buffer[next_free_contact_id] = make_double3(cl_entry->contact->n1_initial[0], cl_entry->contact->n1_initial[1], cl_entry->contact->n1_initial[2]);
            contact_n2_initial_buffer[next_free_contact_id] = make_double3(cl_entry->contact->n2_initial[0], cl_entry->contact->n2_initial[1], cl_entry->contact->n2_initial[2]);

			// add contact id to list of contacts of both particles
			for(int c = 0; c < MAX_CONTACTS; ++c)
			{
				if(contact_id_buffer[CONTACT_ID(c, p, number_of_particles)] < 0)
				{
					contact_id_buffer[CONTACT_ID(c, p, number_of_particles)] = next_free_contact_id;
					particle_id_buffer[CONTACT_ID(c, p, number_of_particles)] = cl_entry->id;
					++contact_number_buffer[p];
					break;
				}
			}

			for(int c = 0; c < MAX_CONTACTS; ++c)
			{

            if(cl_entry->id >= 0) // only contacts with other particles
            {
				if(contact_id_buffer[CONTACT_ID(c, cl_entry->id, number_of_particles)] < 0)
				{
					contact_id_buffer[CONTACT_ID(c, cl_entry->id, number_of_particles)] = next_free_contact_id;
					particle_id_buffer[CONTACT_ID(c, cl_entry->id, number_of_particles)] = p;
					++contact_number_buffer[cl_entry->id];
					break;
				}
            }
			}

            if(cl_entry->id < 0)
                printf("wall bond:  %d  %d\n", cl_entry->contact->id1, cl_entry->contact->id2);


			++next_free_contact_id;
			cl_entry = cl_entry->next;
		}
	}


	number_of_contacts = next_free_contact_id;


	checkCudaErrors(cudaMemcpy(gpu_next_free_contact, &next_free_contact_id, sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(gpu_free_contacts_list, free_contacts_list_buffer, max_contacts * sizeof(int), cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMemcpy(gpu_particle_number_of_contacts, contact_number_buffer, number_of_particles * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_particle_contact_ids, contact_id_buffer, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_particle_particle_ids, particle_id_buffer, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMemcpy(gpu_contact_type, contact_particle_id_buffer, max_contacts * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_contact_particle_ids, contact_particle_id_buffer, max_contacts * sizeof(int2), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_contact_normals_old, contact_normals_old_buffer, max_contacts * sizeof(double4), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(gpu_contact_normals_new, contact_normals_new_buffer, max_contacts * sizeof(double4), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_contact_n1, n1_buffer, max_contacts * sizeof(double4), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_contact_n2, n2_buffer, max_contacts * sizeof(double4), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_contact_rot1, contact_rot1_buffer, max_contacts * sizeof(double4), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(gpu_contact_rot2, contact_rot2_buffer, max_contacts * sizeof(double4), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(gpu_contact_n1_initial, contact_n1_initial_buffer, max_contacts * sizeof(double3), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(gpu_contact_n2_initial, contact_n2_initial_buffer, max_contacts * sizeof(double3), cudaMemcpyHostToDevice));

    delete [] free_contacts_list_buffer;
	delete [] contact_number_buffer;
	delete [] contact_id_buffer;
	delete [] particle_id_buffer;
	delete [] contact_particle_id_buffer;
	delete [] contact_normals_old_buffer;
    delete [] contact_normals_new_buffer;
	delete [] n1_buffer;
	delete [] n2_buffer;
	delete [] contact_rot1_buffer;
	delete [] contact_rot2_buffer;
	delete [] contact_n1_initial_buffer;
	delete [] contact_n2_initial_buffer;

    uint3 gridDim = make_uint3(num_grid_cells, num_grid_cells, num_grid_cells);

    cubPlan.CreateSortPlan(number_of_particles, gridDim, gpu_grid_particle_cell, gpu_grid_particle_index);
    cubPlan.CreateReducePlan(6*number_of_particles);

    if(this->box)
    {
        wallInteraction.initParticleWallInteraction(number_of_particles);
        wallInteraction.initMaterialConstants(
                    this->normal_interaction
                    );

        wallInteraction.load_box_to_gpu(*this);
    }



#ifdef GPU_TRACK_DISSIPATED_ENERGY

    double3 v1;
    double3 v2;
    cudaGetEnergy(cubPlan,
                  gpu_velocities,
                  gpu_dissipated_damping_energy,
                  gpu_dissipated_rolling_energy,
                  gpu_dissipated_contact_energy,
                  gpu_dissipated_sliding_energy,
                  gpu_dissipated_twisting_energy,
                  gpu_kinetic_energy,
                  gpu_rotation_energy,
                  gpu_normal_energy,
                  gpu_rolling_energy,
                  gpu_sliding_energy,
                  gpu_twisting_energy,
                  gpu_velocities,
                  gpu_angular_velocities,
                  gpu_contact_particle_ids,
                  gpu_contact_normals_new,
                  gpu_contact_n1,
                  gpu_contact_n2,
                  E_tot,
                  E_kin,
                  E_rot,
                  V_tot,
                  V_normal,
                  V_rolling,
                  V_sliding,
                  V_twisting,
                  dissipated_damping_energy,
                  dissipated_rolling_energy,
                  dissipated_sliding_energy,
                  dissipated_twisting_energy,
                  dissipated_contact_energy,
                  v1,v2,
                  timestep,
                  number_of_particles
                  );

#endif

	cuda_initialized = true;

	return EC_OK;
}

ErrorCode SimulationCuda::toggleGPUMode(bool use_gpu)
{
	if(use_gpu)
	{
		// enable gpu mode
		if(cuda_initialized)
		{
			this->use_gpu = true;
			return EC_OK;
		}
		else
			return EC_CUDA_NOT_INITIALIZED;
	}
	else
	{
		this->use_gpu = false;
		return EC_OK;
	}
}

ErrorCode SimulationCuda::addParticlesFromSim(Simulation &sim)
{
	cleanUpCuda();

	return Simulation::addParticlesFromSim(sim);
}

ErrorCode SimulationCuda::resizeArrays(int new_number_of_particles, int new_number_of_walls)
{
	cleanUpCuda();

	return Simulation::resizeArrays(new_number_of_particles, new_number_of_walls);
}




void SimulationCuda::GPUprintEnergies()
{
#ifdef GPU_TRACK_DISSIPATED_ENERGY

    double3 v1; // mean velocity of first half of the particles (first particle in case of collision)
    double3 v2; // mean velocity of second half of particles
    cudaGetEnergy(cubPlan,
                  gpu_velocities,
                  gpu_dissipated_damping_energy,
                  gpu_dissipated_rolling_energy,
                  gpu_dissipated_contact_energy,
                  gpu_dissipated_sliding_energy,
                  gpu_dissipated_twisting_energy,
                  gpu_kinetic_energy,
                  gpu_rotation_energy,
                  gpu_normal_energy,
                  gpu_rolling_energy,
                  gpu_sliding_energy,
                  gpu_twisting_energy,
                  gpu_velocities,
                  gpu_angular_velocities,
                  gpu_contact_particle_ids,
                  gpu_contact_normals_new,
                  gpu_contact_n1,
                  gpu_contact_n2,
                  E_tot,
                  E_kin,
                  E_rot,
                  V_tot,
                  V_normal,
                  V_rolling,
                  V_sliding,
                  V_twisting,
                  dissipated_damping_energy,
                  dissipated_rolling_energy,
                  dissipated_sliding_energy,
                  dissipated_twisting_energy,
                  dissipated_contact_energy,
                  v1,
                  v2,
                  timestep,
                  number_of_particles
                  );

    dissipated_contact_energy += (double(broken_contacts) * normal_interaction.getJKRPotentialEnergy(- delta_c) - double(created_contacts) * normal_interaction.getJKRPotentialEnergy(0.0)) * ENERGY_UNIT;
    E_tot += (double(broken_contacts) * normal_interaction.getJKRPotentialEnergy(- delta_c) - double(created_contacts) * normal_interaction.getJKRPotentialEnergy(2*particle_radius - contact_making_dist)) * ENERGY_UNIT;

    if(energies_filename)
    {
        FILE *file = fopen(energies_filename, "a+");

        if(file)
        {

            fprintf(file, "%.12g    %.12g    %.12g    %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g   %.12g\n",
                    current_time,
                    2.0*(double)this->getGPUNumContacts()/(double)this->number_of_particles,
                    E_tot,
                    E_kin,
                    E_rot,

                    V_tot,
                    V_normal,
                    V_rolling,
                    V_sliding,
                    V_twisting,

                    dissipated_rolling_energy,
                    dissipated_sliding_energy,
                    dissipated_twisting_energy,
                    dissipated_damping_energy,
                    dissipated_contact_energy,
                    v1.x,
                    v1.y,
                    v1.z,
                    v2.x,
                    v2.y,
                    v2.z
                    );

            fclose(file);
        }
    }
    else
    {


        //printf("#V_tot, V_normal, V_roll, V_slide, V_twist\n");
        printf("%.2g %.2g %.2g %.2g %.2g %.2g\n",
                V_tot,
                V_normal,
                V_rolling,
                V_sliding,
                V_twisting,
                dissipated_damping_energy
                );


    }
#endif
    return;

}



void SimulationCuda::update()
{
	if(use_gpu)
	{

		// update position and calculate new grid cell
        cudaPosPredictor(
                    gpu_positions_new,
                    gpu_positions_old,
                    gpu_velocities,
                    gpu_forces_old,
                    timestep,
                    gpu_grid_particle_cell,
                    gpu_grid_particle_index,
                    number_of_particles
                    );


		// build new grid
        cubPlan.sort(gpu_grid_particle_cell, gpu_grid_particle_index, number_of_particles);


        cudaFindCellStart(
			gpu_cell_start,
			gpu_cell_end,
			gpu_grid_particle_cell,
			gpu_grid_particle_index,
			gpu_positions_new,
            gpu_positions_sorted,
			number_of_particles,
			num_grid_cells
                    );



		// update contact list and calculate contact normals
		cudaUpdateSticking(
            created_contacts,
            broken_contacts,
#ifdef GPU_TRACK_DISSIPATED_ENERGY
            gpu_dissipated_contact_energy,
            gpu_contact_n1,
            gpu_contact_n2,
#endif // GPU_TRACK_DISSIPATED_ENERGY
			gpu_update_local_contact_list,
			gpu_next_free_contact,
			gpu_free_contacts_list,
			gpu_new_contact_particle_ids,
			gpu_number_of_new_contacts,
			gpu_particle_number_of_contacts,
			gpu_particle_particle_ids,
			gpu_particle_contact_ids,
			gpu_contact_particle_ids,
            gpu_cell_start,
            gpu_cell_end,
            gpu_grid_particle_index,
			gpu_contact_normals_new,
            gpu_positions_new,
            gpu_positions_sorted,
			number_of_particles
            );

		// check if new contacts have to be added
		int number_of_new_contacts;
		cudaMemcpy(&number_of_new_contacts, gpu_number_of_new_contacts, sizeof(int), cudaMemcpyDeviceToHost);


		if(number_of_new_contacts > 0)
		{
            cudaAddNewContacts(
                        gpu_particle_number_of_contacts,
                        gpu_particle_particle_ids,
                        gpu_particle_contact_ids,
                        gpu_next_free_contact,
                        gpu_free_contacts_list,
                        gpu_new_contact_particle_ids,
                        number_of_new_contacts,
                        gpu_contact_particle_ids,
                        gpu_contact_normals_new,
                        gpu_contact_normals_old,
                        gpu_contact_rot1,
                        gpu_contact_rot2,
                        gpu_contact_n1_initial,
                        gpu_contact_n2_initial,
                        gpu_contact_n1,
                        gpu_positions_new,
                        number_of_particles
                        );

			// reset counter
			cudaMemset(gpu_number_of_new_contacts, 0, sizeof(int));
		}
		

        cudaUpdateContacts(
                    gpu_contact_particle_ids,
                    gpu_contact_rot1,
                    gpu_contact_rot2,
                    gpu_contact_n1_initial,
                    gpu_contact_n2_initial,
                    gpu_contact_n1,
                    gpu_contact_n2,
                    gpu_contact_normals_new,
                    gpu_contact_normals_old,
                    gpu_positions_new,
                    gpu_angular_velocities,
                    gpu_torques_old,
                    timestep,
                    number_of_particles
                #ifdef GPU_TRACK_DISSIPATED_ENERGY
                    ,
                    gpu_dissipated_rolling_energy,
                    gpu_dissipated_sliding_energy,
                    gpu_dissipated_twisting_energy
                #endif // GPU_TRACK_DISSIPATED_ENERGY
                    );




		// calculate particle interaction
        cudaUpdateInteraction(
                    gpu_forces_new,
                    gpu_torques_new,
                    gpu_particle_number_of_contacts,
                    gpu_particle_contact_ids,
                    gpu_contact_particle_ids,
                    gpu_contact_n1,
                    gpu_contact_n2,
                    gpu_contact_normals_new,
                    number_of_particles
            #ifdef GPU_TRACK_DISSIPATED_ENERGY
                    ,
                    gpu_dissipated_damping_energy
            #endif
                    );


        //
        if(box)
        {
            //updateBox();
            if(number_of_walls == 6)
            {
                wallInteraction.updateWallSticking(gpu_positions_new,
                    gpu_angular_velocities, gpu_torques_old, timestep);

                wallInteraction.updateWallContactNormals(gpu_positions_new);

                wallInteraction.updateWallContacts(gpu_forces_new, gpu_torques_new, gpu_velocities, gpu_angular_velocities, gpu_torques_old, timestep);

            }
            else if(number_of_walls == 2)
            {
                wallInteraction.updateWallStickingNoSw(gpu_positions_new,
                    gpu_angular_velocities, gpu_torques_old, timestep);

                wallInteraction.updateWallContactNormals(gpu_positions_new);

                wallInteraction.updateWallContacts(gpu_forces_new, gpu_torques_new, gpu_velocities, gpu_angular_velocities, gpu_torques_old, timestep);

            }
        }
        else if (m_sphere_radius > 0.0)
        {
            m_sphere_radius -= m_sphere_compaction_speed*timestep;
            cudaUpdateEnclosingSphereWallInteractionKernel(
                        gpu_forces_new,
                        gpu_positions_new,
                        m_sphere_radius,
                        number_of_particles);
        }

		if(damping_factor == 0)
            cudaCorrector(
                        gpu_velocities,
                        gpu_angular_velocities,
                        gpu_forces_old,
                        gpu_forces_new,
                        gpu_torques_old,
                        gpu_torques_new,
                        timestep,
                        number_of_particles
                        );
		else
            cudaDampingCorrector(
                        gpu_velocities,
                        gpu_angular_velocities,
                        gpu_forces_old,
                        gpu_forces_new,
                        gpu_torques_old,
                        gpu_torques_new,
                        damping_factor,
                        timestep,
                        number_of_particles
                        );


#ifdef GPU_TRACK_PARTICLE_ORIENTATION
		cudaUpdateParticleOrientation(gpu_particle_orientation, gpu_angular_velocities, gpu_torques_old, timestep, number_of_particles);
#endif


		// check filling factor
        if(sim_info.info_storage[0] > 0 && wallInteraction.is_initialized())
		{

            // stop sim if walls are too close
            if(wallInteraction.get_height() < 4.0 * particle_radius)
            {
                printf("Wall too low\n");
                stop_simulation = true;
            }

            if(sim_info.info_storage[7] > 0 && wallInteraction.get_height() < sim_info.info_storage[7])
            {
                printf("Stop height reached\n");
                stop_simulation = true;
            }

            /*
            if(getGPUFillingFactor() > sim_info.info_storage[0])
            {
                printf("Stop filling factor reached %e > %e\n", getGPUFillingFactor(), sim_info.info_storage[0]);
                stop_simulation = true;
            }
            */

		}



        current_time += timestep;


#ifdef GPU_TRACK_DISSIPATED_ENERGY

        // log data if necessary
        if(print_energies_interval > 0)
        {

            if(print_energies_counter >= print_energies_interval - 1)
            {

                print_energies_counter = -1;

                GPUprintEnergies();

            }

            ++print_energies_counter;
        }


#endif


        // check changes of potential energy
        if(check_potential_variation_interval > 0)
        {
            ++check_potential_variation_counter;

            if(check_potential_variation_counter > check_potential_variation_interval)
            {
                check_potential_variation_counter = 0;
			printf("%d	%d	||	%e	%e\n", broken_contacts, created_contacts, E_kin_last, V_tot_last);


            if(broken_contacts > 0)
            {

                if((double(broken_contacts) - E_kin_last + double(created_contacts) - V_tot_last) < 0.5 && current_time > min_end_time)
                {
                    printf("Stopping simulation because no more contacts are broken or created.\n");
                    stop_simulation = true;
                }

                //last_dissipated_energy = dissipated_energy;
                E_kin_last = double(broken_contacts);
                V_tot_last = double(created_contacts);

            //printf("E_diss: %g    Delta E_diss: %g    Threshold: %g    t: %g\n", dissipated_energy, dissipated_energy - last_dissipated_energy, potential_variation_stop_threshold, current_time);

            }
            /*
                if(broken_contacts > 0 && created_contacts > 0)
                {

                    if((double(broken_contacts + created_contacts) - (E_kin_last + V_tot_last))/(E_kin_last + V_tot_last + 1.e-10) < 1.e-4 && current_time > min_end_time)
                    {
                        printf("Stopping simulation because no more contacts are broken or created.\n");
                        stop_simulation = true;
                    }

                    //last_dissipated_energy = dissipated_energy;
                    E_kin_last = double(broken_contacts);
                    V_tot_last = double(created_contacts);

                //printf("E_diss: %g    Delta E_diss: %g    Threshold: %g    t: %g\n", dissipated_energy, dissipated_energy - last_dissipated_energy, potential_variation_stop_threshold, current_time);

                }
                */
            }
        }

        cudaDeviceSynchronize();


        // switch pointers
        double *temp;

        temp = gpu_positions_new;
        gpu_positions_new = gpu_positions_old;
        gpu_positions_old = temp;

        temp = gpu_forces_new;
        gpu_forces_new = gpu_forces_old;
        gpu_forces_old = temp;

        temp = gpu_torques_new;
        gpu_torques_new = gpu_torques_old;
        gpu_torques_old = temp;

        double4 *temp2;
        temp2 = gpu_contact_normals_new;
        gpu_contact_normals_new = gpu_contact_normals_old;
        gpu_contact_normals_old = temp2;

        // check time
        if(current_time > end_time)
        {
            printf("Simulation end_time reached %e > %e\n", current_time, end_time);
            stop_simulation = true;
        }


	}
	else
    {
		Simulation::update();
    }
}

void SimulationCuda::copyPositionsFromGPU(bool new_pos)
{
	size_t memsize = number_of_particles * 3 * sizeof(double);

	// use buffer to convert from cpu to gpu addressing scheme
	double* buffer = new double[3*number_of_particles];

	checkCudaErrors(cudaMemcpy(buffer, gpu_positions_old, memsize, cudaMemcpyDeviceToHost));

	for(int p = 0; p < number_of_particles; ++p)
	{
		pos_old[X_COORD(p)] = (double)buffer[X_COORD_GPU(p, number_of_particles)];
		pos_old[Y_COORD(p)] = (double)buffer[Y_COORD_GPU(p, number_of_particles)]; 
		pos_old[Z_COORD(p)] = (double)buffer[Z_COORD_GPU(p, number_of_particles)];
	}

    if(new_pos)
    {
        checkCudaErrors(cudaMemcpy(buffer, gpu_positions_new, memsize, cudaMemcpyDeviceToHost));
        for(int p = 0; p < number_of_particles; ++p)
        {
            pos_new[X_COORD(p)] = (double)buffer[X_COORD_GPU(p, number_of_particles)];
            pos_new[Y_COORD(p)] = (double)buffer[Y_COORD_GPU(p, number_of_particles)];
            pos_new[Z_COORD(p)] = (double)buffer[Z_COORD_GPU(p, number_of_particles)];
        }
    }
    delete [] buffer;

#ifdef GPU_TRACK_PARTICLE_ORIENTATION
	double4 *orientation_buffer = new double4[number_of_particles];
	checkCudaErrors(cudaMemcpy(orientation_buffer, gpu_particle_orientation, number_of_particles * sizeof(double4), cudaMemcpyDeviceToHost));

	for(int p = 0; p < number_of_particles; ++p)
	{
		orientation[4*p] = orientation_buffer[p].w;
		orientation[4*p+1] = orientation_buffer[p].x;
		orientation[4*p+2] = orientation_buffer[p].y;
		orientation[4*p+3] = orientation_buffer[p].z;
	}

	delete [] orientation_buffer;
#endif
}

void SimulationCuda::copySimDataFromGPU()
{
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// copy particle data
	//////////////////////////////////////////////////////////////////////////////////////////////////

	// use buffer to convert from cpu to gpu addressing scheme
	double* buffer = new double[3*number_of_particles];

    copyParticleArrayFromGPU(pos_old, gpu_positions_old, number_of_particles, buffer);

    copyParticleArrayFromGPU(pos_new, gpu_positions_new, number_of_particles, buffer);

    copyParticleArrayFromGPU(vel, gpu_velocities, number_of_particles, buffer);

    copyParticleArrayFromGPU(vel_angular, gpu_angular_velocities, number_of_particles, buffer);

    copyParticleArrayFromGPU(force_old, gpu_forces_old, number_of_particles, buffer);

    copyParticleArrayFromGPU(force_new, gpu_forces_new, number_of_particles, buffer);

    copyParticleArrayFromGPU(torque_old, gpu_torques_old, number_of_particles, buffer);

    copyParticleArrayFromGPU(torque_new, gpu_torques_new, number_of_particles, buffer);








#ifdef GPU_TRACK_DISSIPATED_ENERGY
#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE


    for(int p = 0; p < number_of_particles; ++p)
    {
        dissipated_energy_of_particle[p] = 0.0;
    }

    checkCudaErrors(cudaMemcpy(buffer, gpu_dissipated_damping_energy, number_of_particles * sizeof(double), cudaMemcpyDeviceToHost));
    for(int p = 0; p < number_of_particles; ++p)
    {
        dissipated_energy_of_particle[p] += buffer[p];
    }

    checkCudaErrors(cudaMemcpy(buffer, gpu_dissipated_rolling_energy, number_of_particles * sizeof(double), cudaMemcpyDeviceToHost));
    for(int p = 0; p < number_of_particles; ++p)
    {
        dissipated_energy_of_particle[p] += buffer[p];
    }

    checkCudaErrors(cudaMemcpy(buffer, gpu_dissipated_contact_energy, number_of_particles * sizeof(double), cudaMemcpyDeviceToHost));
    for(int p = 0; p < number_of_particles; ++p)
    {
        dissipated_energy_of_particle[p] += buffer[p];
    }

    checkCudaErrors(cudaMemcpy(buffer, gpu_dissipated_sliding_energy, number_of_particles * sizeof(double), cudaMemcpyDeviceToHost));
    for(int p = 0; p < number_of_particles; ++p)
    {
        dissipated_energy_of_particle[p] += buffer[p];
    }

    checkCudaErrors(cudaMemcpy(buffer, gpu_dissipated_twisting_energy, number_of_particles * sizeof(double), cudaMemcpyDeviceToHost));
    for(int p = 0; p < number_of_particles; ++p)
    {
        dissipated_energy_of_particle[p] += buffer[p];
    }


#endif
#endif





	delete [] buffer;

#ifdef GPU_TRACK_PARTICLE_ORIENTATION
	double4 *orientation_buffer = new double4[number_of_particles];
	cudaMemcpy(orientation_buffer, gpu_particle_orientation, number_of_particles * sizeof(double4), cudaMemcpyDeviceToHost);

	for(int p = 0; p < number_of_particles; ++p)
	{
		orientation[4*p] = orientation_buffer[p].w;
		orientation[4*p+1] = orientation_buffer[p].x;
		orientation[4*p+2] = orientation_buffer[p].y;
		orientation[4*p+3] = orientation_buffer[p].z;
	}

	delete [] orientation_buffer;
#endif

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// copy contact list
	//////////////////////////////////////////////////////////////////////////////////////////////////

	// delete old contact list (if existing)
	deleteContacts();

	int max_contacts = 6 * number_of_particles;

	int *contact_id_buffer = new int[MAX_CONTACTS * number_of_particles];
	int2 *contact_particle_id_buffer = new int2[max_contacts];
	double4 *contact_n1_buffer = new double4[max_contacts];
	double4 *contact_n2_buffer = new double4[max_contacts];
	double4 *contact_rot1_buffer = new double4[max_contacts];
	double4 *contact_rot2_buffer = new double4[max_contacts];
	double3 *contact_n1_initial_buffer = new double3[max_contacts];
	double3 *contact_n2_initial_buffer = new double3[max_contacts];
	double4 *contact_normals_old_buffer = new double4[max_contacts];

	cudaMemcpy(contact_id_buffer, gpu_particle_contact_ids, MAX_CONTACTS * number_of_particles * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_particle_id_buffer, gpu_contact_particle_ids, max_contacts * sizeof(int2), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_n1_buffer, gpu_contact_n1, max_contacts * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_n2_buffer, gpu_contact_n2, max_contacts * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_rot1_buffer, gpu_contact_rot1, max_contacts * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_rot2_buffer, gpu_contact_rot2, max_contacts * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_n1_initial_buffer, gpu_contact_n1_initial, max_contacts * sizeof(double3), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_n2_initial_buffer, gpu_contact_n2_initial, max_contacts * sizeof(double3), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_normals_old_buffer, gpu_contact_normals_old, max_contacts * sizeof(double4), cudaMemcpyDeviceToHost);

	for(int p = 0; p < number_of_particles; ++p)
	{
		//entry = contact_list[p];
		for(int c = 0; c < MAX_CONTACTS; ++c)
		{
			int contact_id = contact_id_buffer[CONTACT_ID(c, p, number_of_particles)];

			if(contact_id < 0)
				break;

			int2 particle_id = contact_particle_id_buffer[contact_id];

			if(particle_id.x == p)
			{
				ContactListEntry *cl_entry = getContactInsertPos(p, particle_id.y);

				if(cl_entry)
				{
					// create new contact
					cl_entry->contact = new Contact;

					cl_entry->contact->id1 = particle_id.x;
					cl_entry->contact->id2 = particle_id.y;

					cl_entry->contact->n1_initial[0] = contact_n1_initial_buffer[contact_id].x;
					cl_entry->contact->n1_initial[1] = contact_n1_initial_buffer[contact_id].y;
					cl_entry->contact->n1_initial[2] = contact_n1_initial_buffer[contact_id].z;
					cl_entry->contact->n2_initial[0] = contact_n2_initial_buffer[contact_id].x;
					cl_entry->contact->n2_initial[1] = contact_n2_initial_buffer[contact_id].y;
					cl_entry->contact->n2_initial[2] = contact_n2_initial_buffer[contact_id].z;

					cl_entry->contact->rot1.e0 = contact_rot1_buffer[contact_id].w;
					cl_entry->contact->rot1.e1 = contact_rot1_buffer[contact_id].x;
					cl_entry->contact->rot1.e2 = contact_rot1_buffer[contact_id].y;
					cl_entry->contact->rot1.e3 = contact_rot1_buffer[contact_id].z;
					cl_entry->contact->rot2.e0 = contact_rot2_buffer[contact_id].w;
					cl_entry->contact->rot2.e1 = contact_rot2_buffer[contact_id].x;
					cl_entry->contact->rot2.e2 = contact_rot2_buffer[contact_id].y;
					cl_entry->contact->rot2.e3 = contact_rot2_buffer[contact_id].z;

					cl_entry->contact->twisting_displacement = contact_n1_buffer[contact_id].w;


                    double compression_length = 2.0 * particle_radius - contact_normals_old_buffer[contact_id].w;

                    cl_entry->contact->compression_length = compression_length;

					cl_entry->contact->old_contact_normal[0] = contact_normals_old_buffer[contact_id].x;
					cl_entry->contact->old_contact_normal[1] = contact_normals_old_buffer[contact_id].y;
					cl_entry->contact->old_contact_normal[2] = contact_normals_old_buffer[contact_id].z;

					++number_of_contacts;
				}
			}
		}
	}

	delete [] contact_id_buffer;
	delete [] contact_particle_id_buffer;
	delete [] contact_n1_buffer;
	delete [] contact_n2_buffer;
	delete [] contact_rot1_buffer;
	delete [] contact_rot2_buffer;
	delete [] contact_n1_initial_buffer;
	delete [] contact_n2_initial_buffer;
	delete [] contact_normals_old_buffer;


    if(box)
        wallInteraction.load_box_from_gpu(*this);


}

void SimulationCuda::printGPUContactList(const char *filename)
{
	FILE *file = fopen(filename, "w+");

	if(!file)
		return;

	int *particle_particle_ids_buffer = new int[number_of_particles * MAX_CONTACTS];
	int *particle_contact_ids_buffer = new int[number_of_particles * MAX_CONTACTS];
	cudaMemcpy(particle_particle_ids_buffer, gpu_particle_particle_ids, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(particle_contact_ids_buffer, gpu_particle_contact_ids, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyDeviceToHost);

	for(int p = 0; p < number_of_particles; ++p)
	{
		for(int c = 0; c < MAX_CONTACTS; ++c)
		{
			int p2 = particle_particle_ids_buffer[CONTACT_ID(c, p, number_of_particles)];

			//if(p2 != NO_PARTICLE && p2 < p)
			fprintf(file, "%i <-> %i: %i\n", p, p2, particle_contact_ids_buffer[CONTACT_ID(c, p, number_of_particles)]);
		}

		fprintf(file, "\n");
	}

	delete [] particle_particle_ids_buffer;
	delete [] particle_contact_ids_buffer;

	fprintf(file, "\n\n");

	int2 *contact_particle_ids_buffer = new int2[number_of_particles * 6];
	cudaMemcpy(contact_particle_ids_buffer, gpu_contact_particle_ids, number_of_particles * 6 * sizeof(int2), cudaMemcpyDeviceToHost);

	for(int c = 0; c < number_of_particles * 6; ++c)
		fprintf(file, "%i: %i <-> %i\n", c, contact_particle_ids_buffer[c].x, contact_particle_ids_buffer[c].y);

	delete [] contact_particle_ids_buffer;

	fprintf(file, "\n\n");

	int next_free_id;
	cudaMemcpy(&next_free_id, gpu_next_free_contact, sizeof(int), cudaMemcpyDeviceToHost);

	fprintf(file, "next_id: %i\n", next_free_id);

	int *free_contacts_list_buffer = new int[number_of_particles * 6];
	cudaMemcpy(free_contacts_list_buffer, gpu_free_contacts_list, number_of_particles * 6 * sizeof(int), cudaMemcpyDeviceToHost);

	for(int c = 0; c < number_of_particles * 6; ++c)
		fprintf(file, "%-3i", free_contacts_list_buffer[c]);
	fprintf(file, "\n");

	delete [] free_contacts_list_buffer;

	fclose(file);
}

void SimulationCuda::printGPUContactNormals(const char *filename)
{
	int max_contacts = 6 * number_of_particles;

	int *contact_id_buffer = new int[MAX_CONTACTS * number_of_particles];
	int2 *contact_particle_id_buffer = new int2[max_contacts];
	double4 *contact_normals_new = new double4[max_contacts];
	double4 *contact_normals_old = new double4[max_contacts];

	cudaMemcpy(contact_id_buffer, gpu_particle_contact_ids, MAX_CONTACTS * number_of_particles * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_particle_id_buffer, gpu_contact_particle_ids, max_contacts * sizeof(int2), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_normals_new, gpu_contact_normals_new, max_contacts * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(contact_normals_old, gpu_contact_normals_old, max_contacts * sizeof(double4), cudaMemcpyDeviceToHost);

	FILE *file = fopen(filename, "w+");

	if(file)
	{
		for(int p = 0; p < number_of_particles; ++p)
		{
			for(int c = 0; c < MAX_CONTACTS; ++c)
			{
				int contact_id = contact_id_buffer[CONTACT_ID(c, p, number_of_particles)];

				if(contact_id < 0)
					break;

				int2 particle_id = contact_particle_id_buffer[contact_id];

				if(p == particle_id.x)
				{
					fprintf(file, "%i %i: n_c_new = ( %.16lf , %.16lf , %.16lf ), dist_new = %.16lg   n_c_old = ( %.16lf , %.16lf , %.16lf ), dist_old = %.16lg\n", particle_id.x, particle_id.y, 
						contact_normals_new[contact_id].x, contact_normals_new[contact_id].y, contact_normals_new[contact_id].z, contact_normals_new[contact_id].w,
						contact_normals_old[contact_id].x, contact_normals_old[contact_id].y, contact_normals_old[contact_id].z, contact_normals_old[contact_id].w);
				}
			}
		}

		fclose(file);
	}

	delete [] contact_id_buffer;
	delete [] contact_particle_id_buffer;
	delete [] contact_normals_new;
	delete [] contact_normals_old;
}

void SimulationCuda::printGPUContacts(const char *filename)
{
	double4 *n1 = new double4[number_of_particles * 6];
	double4 *n2 = new double4[number_of_particles * 6];
	double3 *n1_initial = new double3[number_of_particles * 6];
	double3 *n2_initial = new double3[number_of_particles * 6];
	double4 *rot1 = new double4[number_of_particles * 6];
	double4 *rot2 = new double4[number_of_particles * 6];
	double4 *normals_new = new double4[number_of_particles * 6];
	double4 *normals_old = new double4[number_of_particles * 6];

	int *particle_contact_ids = new int[number_of_particles*MAX_CONTACTS];
	int *particle_particle_ids = new int[number_of_particles*MAX_CONTACTS];

	cudaMemcpy(n1, gpu_contact_n1, number_of_particles * 6 * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(n2, gpu_contact_n2, number_of_particles * 6 * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(n1_initial, gpu_contact_n1_initial, number_of_particles * 6 * sizeof(double3), cudaMemcpyDeviceToHost);
	cudaMemcpy(n2_initial, gpu_contact_n2_initial, number_of_particles * 6 * sizeof(double3), cudaMemcpyDeviceToHost);
	cudaMemcpy(rot1, gpu_contact_rot1, number_of_particles * 6 * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(rot2, gpu_contact_rot2, number_of_particles * 6 * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(normals_old, gpu_contact_normals_old, number_of_particles * 6 * sizeof(double4), cudaMemcpyDeviceToHost);
	cudaMemcpy(normals_new, gpu_contact_normals_new, number_of_particles * 6 * sizeof(double4), cudaMemcpyDeviceToHost);
	
	cudaMemcpy(particle_contact_ids, gpu_particle_contact_ids, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(particle_particle_ids, gpu_particle_particle_ids, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyDeviceToHost);


	FILE *file = fopen(filename, "w+");

	if(file)
	{
		for(int p = 0; p < number_of_particles; ++p)
		{
			for(int c = 0; c < MAX_CONTACTS; ++c)
			{	
				int p2 = particle_particle_ids[CONTACT_ID(c, p, number_of_particles)];

				if(p2 != NO_PARTICLE && p > p2)
				{
					int contact_id = particle_contact_ids[CONTACT_ID(c, p, number_of_particles)];

					fprintf(file, "%i %i: n1 = ( %.18lg, %.18lg, %.18lg ), n2 = ( %.18lg, %.18lg, %.18lg ), n_c_new = ( %.16lf , %.16lf , %.16lf ), dist = %.18lg\n", 
						p, particle_particle_ids[CONTACT_ID(c,p, number_of_particles)],
						n1[contact_id].x, n1[contact_id].y, n1[contact_id].z,
						n2[contact_id].x, n2[contact_id].y, n2[contact_id].z,
						normals_old[contact_id].x, normals_old[contact_id].y, normals_old[contact_id].z, normals_old[contact_id].w);

					/*fprintf(file, "%i %i: n1 = ( %.16lg, %.16lg, %.16lg ), n2 = ( %.16lg, %.16lg, %.16lg ), n1_i = ( %.16lg, %.16lg, %.16lg ), n2_i= ( %.16lg, %.16lg, %.16lg ), rot1 = ( %.16lf, %.16lf, %.16lf, %.16lf ), rot2 = ( %.16lf, %.16lf, %.16lf, %.16lf ), Phi = %.16lg,    n_c_old = ( %.16lf , %.16lf , %.16lf ), n_c_new = ( %.16lf , %.16lf , %.16lf )\n", 
						p, contact_ids[CONTACT_ID(c, p, number_of_particles)], 
						n1[CONTACT_ID(c, p, number_of_particles)].x, n1[CONTACT_ID(c, p, number_of_particles)].y, n1[CONTACT_ID(c, p, number_of_particles)].z,
						n2[CONTACT_ID(c, p, number_of_particles)].x, n2[CONTACT_ID(c, p, number_of_particles)].y, n2[CONTACT_ID(c, p, number_of_particles)].z,
						n1_initial[CONTACT_ID(c, p, number_of_particles)].x, n1_initial[CONTACT_ID(c, p, number_of_particles)].y, n1_initial[CONTACT_ID(c, p, number_of_particles)].z,
						n2_initial[CONTACT_ID(c, p, number_of_particles)].x, n2_initial[CONTACT_ID(c, p, number_of_particles)].y, n2_initial[CONTACT_ID(c, p, number_of_particles)].z,
						rot1[CONTACT_ID(c, p, number_of_particles)].w, rot1[CONTACT_ID(c, p, number_of_particles)].x, rot1[CONTACT_ID(c, p, number_of_particles)].y, rot1[CONTACT_ID(c, p, number_of_particles)].z,
						rot2[CONTACT_ID(c, p, number_of_particles)].w, rot2[CONTACT_ID(c, p, number_of_particles)].x, rot2[CONTACT_ID(c, p, number_of_particles)].y, rot2[CONTACT_ID(c, p, number_of_particles)].z,
						n1_initial[CONTACT_ID(c, p, number_of_particles)].w,
						normals_old[CONTACT_ID(c, p, number_of_particles)].x, normals_old[CONTACT_ID(c, p, number_of_particles)].y, normals_old[CONTACT_ID(c, p, number_of_particles)].z,
						normals_new[CONTACT_ID(c, p, number_of_particles)].x, normals_new[CONTACT_ID(c, p, number_of_particles)].y, normals_new[CONTACT_ID(c, p, number_of_particles)].z);
				*/
				}
			}
		}

		fclose(file);
	}

	delete [] normals_new;
	delete [] normals_old;
	delete [] n1;
	delete [] n2;
	delete [] n1_initial;
	delete [] n2_initial;
	delete [] rot1;
	delete [] rot2;

	delete [] particle_contact_ids;
	delete [] particle_particle_ids;
}

void SimulationCuda::printGPUForces(const char *filename)
{
	double *force_buffer = new double[3*number_of_particles];
	double *torque_buffer = new double[3*number_of_particles];
	cudaMemcpy(force_buffer, gpu_velocities, number_of_particles * 3 * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(torque_buffer, gpu_torques_old, number_of_particles * 3 * sizeof(double), cudaMemcpyDeviceToHost);

	FILE *file = fopen(filename, "w+");

	if(file)
	{
		fprintf(file, "# Interaction:");

	#ifdef GPU_ROLLING
		fprintf(file, "Rolling ");
	#endif
	#ifdef GPU_SLIDING
		fprintf(file, "Sliding ");
	#endif
	#ifdef GPU_TWISTING
		fprintf(file, "Twisting ");
	#endif
		fprintf(file, "\n");

		for(int p = 0; p < number_of_particles; ++p)
			fprintf(file, "%i %g %g %g %g %g %g\n", p, force_buffer[X_COORD_GPU(p, number_of_particles)], force_buffer[Y_COORD_GPU(p, number_of_particles)], force_buffer[Z_COORD_GPU(p, number_of_particles)], torque_buffer[X_COORD_GPU(p, number_of_particles)] * moment_of_inertia, torque_buffer[Y_COORD_GPU(p, number_of_particles)] * moment_of_inertia, torque_buffer[Z_COORD_GPU(p, number_of_particles)]  * moment_of_inertia);	

		fclose(file);
	}

	delete [] force_buffer;
	delete [] torque_buffer;
}





ErrorCode SimulationCuda::printGPUsimToFile(const char* filename)
{
    // try to open specified file
    FILE *file = fopen(filename, "wb+");

    if(!file)
        return EC_FILE_NOT_FOUND;





    const size_t memsize = number_of_particles * sizeof(double3);
    const int max_contacts = 6 * number_of_particles;


    double3* buffer = new double3[number_of_particles];



    // write file version, simulation type, etc.
    fwrite(&number_of_particles, sizeof(int), 1, file);
    fwrite(&number_of_walls, sizeof(int), 1, file);


    fwrite(&sim_info, sizeof(SimInfo), 1, file);
    fwrite(&timestep, sizeof(double), 1, file);


    fwrite(&current_time, sizeof(double), 1, file);
    fwrite(&end_time, sizeof(double), 1, file);
    fwrite(&stop_simulation, sizeof(bool), 1, file);
    fwrite(&use_gpu, sizeof(bool), 1, file);



    GPUConstants d_gpu_constants;
    getGPUConstants(d_gpu_constants);
    fwrite(&d_gpu_constants, sizeof(GPUConstants), 1, file);



    // particle data
    checkCudaErrors(cudaMemcpy(buffer, gpu_positions_old, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_positions_new, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_positions_sorted, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_velocities, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_angular_velocities, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_forces_old, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_forces_new, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_torques_old, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(buffer, gpu_torques_new, memsize, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double3), number_of_particles, file);

    delete[] buffer;


    // contacts data
    int* int_buffer = new int[number_of_particles * MAX_CONTACTS];

    checkCudaErrors(cudaMemcpy(int_buffer, gpu_particle_number_of_contacts, number_of_particles * sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(int_buffer, sizeof(int), number_of_particles, file);

    checkCudaErrors(cudaMemcpy(int_buffer, gpu_particle_contact_ids, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(int_buffer, sizeof(int), number_of_particles * MAX_CONTACTS, file);

    checkCudaErrors(cudaMemcpy(int_buffer, gpu_particle_particle_ids, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(int_buffer, sizeof(int), number_of_particles * MAX_CONTACTS, file);

    checkCudaErrors(cudaMemcpy(int_buffer, gpu_update_local_contact_list, number_of_particles * sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(int_buffer, sizeof(int), number_of_particles, file);

    int2* int2_buffer = new int2[max_contacts];
    checkCudaErrors(cudaMemcpy(int2_buffer, gpu_contact_particle_ids, max_contacts*sizeof(int2), cudaMemcpyDeviceToHost));
    fwrite(int2_buffer, sizeof(int2), max_contacts, file);

    checkCudaErrors(cudaMemcpy(int_buffer, gpu_contact_type, max_contacts*sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(int_buffer, sizeof(int), max_contacts, file);



    double4* double4_buffer = new double4[max_contacts];

    checkCudaErrors(cudaMemcpy(double4_buffer, gpu_contact_n1, max_contacts*sizeof(double4), cudaMemcpyDeviceToHost));
    fwrite(double4_buffer, sizeof(double4), max_contacts, file);

    checkCudaErrors(cudaMemcpy(double4_buffer, gpu_contact_n2, max_contacts*sizeof(double4), cudaMemcpyDeviceToHost));
    fwrite(double4_buffer, sizeof(double4), max_contacts, file);

    checkCudaErrors(cudaMemcpy(double4_buffer, gpu_contact_normals_old, max_contacts*sizeof(double4), cudaMemcpyDeviceToHost));
    fwrite(double4_buffer, sizeof(double4), max_contacts, file);

    checkCudaErrors(cudaMemcpy(double4_buffer, gpu_contact_normals_new, max_contacts*sizeof(double4), cudaMemcpyDeviceToHost));
    fwrite(double4_buffer, sizeof(double4), max_contacts, file);

    checkCudaErrors(cudaMemcpy(double4_buffer, gpu_contact_rot1, max_contacts*sizeof(double4), cudaMemcpyDeviceToHost));
    fwrite(double4_buffer, sizeof(double4), max_contacts, file);

    checkCudaErrors(cudaMemcpy(double4_buffer, gpu_contact_rot2, max_contacts*sizeof(double4), cudaMemcpyDeviceToHost));
    fwrite(double4_buffer, sizeof(double4), max_contacts, file);


    double3* double3_buffer = new double3[max_contacts];

    checkCudaErrors(cudaMemcpy(double3_buffer, gpu_contact_n1_initial, max_contacts*sizeof(double3), cudaMemcpyDeviceToHost));
    fwrite(double3_buffer, sizeof(double3), max_contacts, file);

    checkCudaErrors(cudaMemcpy(double3_buffer, gpu_contact_n2_initial, max_contacts*sizeof(double3), cudaMemcpyDeviceToHost));
    fwrite(double3_buffer, sizeof(double3), max_contacts, file);



    int* integer = new int[1];

    checkCudaErrors(cudaMemcpy(integer, gpu_next_free_contact, sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(integer, sizeof(int), 1, file);

    checkCudaErrors(cudaMemcpy(integer, gpu_number_of_new_contacts, sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(integer, sizeof(int), 1, file);


    checkCudaErrors(cudaMemcpy(int2_buffer, gpu_new_contact_particle_ids, max_contacts*sizeof(int2), cudaMemcpyDeviceToHost));
    fwrite(int2_buffer, sizeof(int2), max_contacts, file);

    checkCudaErrors(cudaMemcpy(int_buffer, gpu_free_contacts_list, max_contacts * sizeof(int), cudaMemcpyDeviceToHost));
    fwrite(int_buffer, sizeof(int), max_contacts, file);


    delete[] double4_buffer;
    delete[] double3_buffer;
    delete[] integer;
    delete[] int2_buffer;
    delete[] int_buffer;


    if(this->wallInteraction.is_initialized())
    {
        fwrite(box, sizeof(WallBox), 1, file);

        for(int w = 0; w < number_of_walls; ++w)
        {
            fwrite(&(walls[w]),  sizeof(Wall), 1, file);
        }

        this->wallInteraction.storeToFile(file);
    }


#ifdef GPU_TRACK_DISSIPATED_ENERGY
    const size_t memsize_double = number_of_particles * sizeof(double);

    double* double_buffer = new double[number_of_particles];

    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_damping_energy, memsize_double, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_rolling_energy, memsize_double, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_contact_energy, memsize_double, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_sliding_energy, memsize_double, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_twisting_energy, memsize_double, cudaMemcpyDeviceToHost));
    fwrite(buffer, sizeof(double), number_of_particles, file);

#endif


    fclose(file);


    return EC_OK;
}



ErrorCode SimulationCuda::loadGPUsimFromFile(const char* filename, const int GPU_id)
{
    // try to open specified file
    FILE *file = fopen(filename, "rb");

    if(!file)
        return EC_FILE_NOT_FOUND;




    // write file version, simulation type, etc.
    fread(&number_of_particles, sizeof(int), 1, file);
    fread(&number_of_walls, sizeof(int), 1, file);

    printf("%d  %d\n", number_of_particles, number_of_walls);


    resizeArrays(number_of_particles, number_of_walls);
    this->initCuda(GPU_id);
    this->toggleGPUMode(true);



    fread(&sim_info, sizeof(SimInfo), 1, file);
    fread(&timestep, sizeof(double), 1, file);

    fread(&current_time, sizeof(double), 1, file);
    fread(&end_time, sizeof(double), 1, file);
    fread(&stop_simulation, sizeof(bool), 1, file);
    fread(&use_gpu, sizeof(bool), 1, file);

    printf("%e  %e  %e\n", timestep, current_time, end_time);

    const size_t memsize = number_of_particles * sizeof(double3);
    const int max_contacts = 6 * number_of_particles;



    GPUConstants d_gpu_constants;
    fread(&d_gpu_constants, sizeof(GPUConstants), 1, file);
    setGPUConstants(d_gpu_constants);


    // particle data
    double3* buffer = new double3[number_of_particles];

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_positions_old, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_positions_new, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_positions_sorted, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_velocities, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_angular_velocities, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_forces_old, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_forces_new, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_torques_old, buffer, memsize, cudaMemcpyHostToDevice));

    fread(buffer, sizeof(double3), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_torques_new, buffer, memsize, cudaMemcpyHostToDevice));

    delete[] buffer;


    // contacts data
    int* int_buffer = new int[number_of_particles * MAX_CONTACTS];

    fread(int_buffer, sizeof(int), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_particle_number_of_contacts, int_buffer, number_of_particles * sizeof(int), cudaMemcpyHostToDevice));

    fread(int_buffer, sizeof(int), number_of_particles * MAX_CONTACTS, file);
    checkCudaErrors(cudaMemcpy(gpu_particle_contact_ids, int_buffer, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyHostToDevice));

    fread(int_buffer, sizeof(int), number_of_particles * MAX_CONTACTS, file);
    checkCudaErrors(cudaMemcpy(gpu_particle_particle_ids, int_buffer, number_of_particles * MAX_CONTACTS * sizeof(int), cudaMemcpyHostToDevice));

    fread(int_buffer, sizeof(int), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(gpu_update_local_contact_list, int_buffer, number_of_particles * sizeof(int), cudaMemcpyHostToDevice));


    int2* int2_buffer = new int2[max_contacts];

    fread(int2_buffer, sizeof(int2), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_particle_ids, int2_buffer, max_contacts*sizeof(int2), cudaMemcpyHostToDevice));

    fread(int_buffer, sizeof(int), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_type, int_buffer, max_contacts*sizeof(int), cudaMemcpyHostToDevice));



    double4* double4_buffer = new double4[max_contacts];

    fread(double4_buffer, sizeof(double4), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_n1, double4_buffer, max_contacts*sizeof(double4), cudaMemcpyHostToDevice));

    fread(double4_buffer, sizeof(double4), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_n2, double4_buffer, max_contacts*sizeof(double4), cudaMemcpyHostToDevice));

    fread(double4_buffer, sizeof(double4), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_normals_old, double4_buffer, max_contacts*sizeof(double4), cudaMemcpyHostToDevice));

    fread(double4_buffer, sizeof(double4), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_normals_new, double4_buffer, max_contacts*sizeof(double4), cudaMemcpyHostToDevice));

    fread(double4_buffer, sizeof(double4), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_rot1, double4_buffer, max_contacts*sizeof(double4), cudaMemcpyHostToDevice));

    fread(double4_buffer, sizeof(double4), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_rot2, double4_buffer, max_contacts*sizeof(double4), cudaMemcpyHostToDevice));



    double3* double3_buffer = new double3[max_contacts];

    fread(double3_buffer, sizeof(double3), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_n1_initial, double3_buffer, max_contacts*sizeof(double3), cudaMemcpyHostToDevice));

    fread(double3_buffer, sizeof(double3), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_contact_n2_initial, double3_buffer, max_contacts*sizeof(double3), cudaMemcpyHostToDevice));




    int* integer = new int[1];

    fread(integer, sizeof(int), 1, file);
    checkCudaErrors(cudaMemcpy(gpu_next_free_contact, integer, sizeof(int), cudaMemcpyHostToDevice));

    fread(integer, sizeof(int), 1, file);
    checkCudaErrors(cudaMemcpy(gpu_number_of_new_contacts, integer, sizeof(int), cudaMemcpyHostToDevice));

    fread(int2_buffer, sizeof(int2), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_new_contact_particle_ids, int2_buffer, max_contacts*sizeof(int2), cudaMemcpyHostToDevice));

    fread(int_buffer, sizeof(int), max_contacts, file);
    checkCudaErrors(cudaMemcpy(gpu_free_contacts_list, int_buffer, max_contacts * sizeof(int), cudaMemcpyHostToDevice));



    delete[] double4_buffer;
    delete[] double3_buffer;
    delete[] integer;
    delete[] int2_buffer;
    delete[] int_buffer;


    if(number_of_walls > 0)
    {


        WallBox* temp_box = new WallBox;
        fread(temp_box, sizeof(WallBox), 1, file);
        box = temp_box;

        for(int w = 0; w < number_of_walls; ++w)
        {
            fread(&(walls[w]), sizeof(Wall), 1, file);
        }

        this->wallInteraction.initParticleWallInteraction(number_of_particles);
        this->wallInteraction.loadFromFile(file);
    }



#ifdef GPU_TRACK_DISSIPATED_ENERGY
    const size_t memsize_double = number_of_particles * sizeof(double);

    double* double_buffer = new double[number_of_particles];
    /*
    if(!feof(file))
    {
        delete[] double_buffer;

        printf("WARNING: tried to read energies from empty file!\n");
        fclose(file);


        return EC_OK;
    }
    */
    fread(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_damping_energy, memsize_double, cudaMemcpyDeviceToHost));

    /*
    if(!feof(file))
    {
        delete[] double_buffer;

        printf("WARNING: tried to read energies from empty file!\n");
        fclose(file);


        return EC_OK;
    }
    */

    fread(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_rolling_energy, memsize_double, cudaMemcpyDeviceToHost));

    fread(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_contact_energy, memsize_double, cudaMemcpyDeviceToHost));

    fread(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_sliding_energy, memsize_double, cudaMemcpyDeviceToHost));

    fread(buffer, sizeof(double), number_of_particles, file);
    checkCudaErrors(cudaMemcpy(double_buffer, gpu_dissipated_twisting_energy, memsize_double, cudaMemcpyDeviceToHost));

    delete[] double_buffer;

#endif



    fclose(file);


    return EC_OK;
}




