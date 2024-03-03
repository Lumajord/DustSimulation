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
 *   Free Software Foundation, Inc.,										*
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.				*
 ***************************************************************************/

#ifndef SIMULATIONCUDA_H
#define SIMULATIONCUDA_H

#include "cuda.h"
#include "cuda_runtime_api.h"

#include "Simulation.h"
#include "CudaDefines.cuh"
#include "myCubWrapper.h"
#include "CudaWallInteraction.h"

#include "helper_timer.h"


class SimulationCuda : public Simulation
{
public:
	SimulationCuda(void);
	~SimulationCuda(void);


    particleWallInteraction wallInteraction;
    CubPlan cubPlan;

	/////////////////////////////////////////////////////////////////////////////////////////
	// cuda specific functions
	/////////////////////////////////////////////////////////////////////////////////////////

	ErrorCode initCuda(int GPU_id = -1);

	void cleanUpCuda();

	bool sufficientMemory(size_t required_memory);

	ErrorCode toggleGPUMode(bool use_gpu);

	void setDampingFactor(double damping_factor);

	double getGPUFillingFactor();
    double getGPUFillingFactorNoSW();

    void copyPositionsFromGPU(bool new_pos = true);

	void copySimDataFromGPU();

	void printGPUContactList(const char *filename);

	void printGPUContactNormals(const char *filename);

	void printGPUContacts(const char *filename);

	void printGPUForces(const char *filename);

    void GPUprintEnergies();

    bool get_use_gpu() {return this->use_gpu;}
    void set_use_gpu(bool use_gpu) {this->use_gpu = use_gpu;}

    int getGPUNumContacts();

    ErrorCode printGPUsimToFile(const char* filename);
    ErrorCode loadGPUsimFromFile(const char* filename, const int GPU_id);


	/////////////////////////////////////////////////////////////////////////////////////////
	// extended functions
	/////////////////////////////////////////////////////////////////////////////////////////

	void update();

	ErrorCode addParticlesFromSim(Simulation &sim);

	ErrorCode resizeArrays(int new_number_of_particles, int new_number_of_walls);

	/////////////////////////////////////////////////////////////////////////////////////////
	// variables
	/////////////////////////////////////////////////////////////////////////////////////////

    double m_sphere_radius;
    double m_sphere_compaction_speed;

private:


    void copyParticleArrayFromGPU(double* host_ptr, double* device_ptr, int number_of_particles, double* host_buffer);
    void copyParticleArrayToGPU(double* device_ptr, double* host_ptr, int number_of_particles, double* host_buffer);

    cudaDeviceProp device_properties;

	// particle data
    double* gpu_positions_old;
    double* gpu_positions_new;
    double* gpu_positions_sorted;

    double* gpu_velocities;
    double* gpu_angular_velocities;
    double* gpu_forces_old;
    double* gpu_forces_new;
    double* gpu_torques_old;
    double* gpu_torques_new;

#ifdef GPU_TRACK_PARTICLE_ORIENTATION
	double4 *gpu_particle_orientation;
#endif

#ifdef GPU_TRACK_DISSIPATED_ENERGY


    double *gpu_dissipated_damping_energy;
    double *gpu_dissipated_rolling_energy;
    double *gpu_dissipated_contact_energy;
    double *gpu_dissipated_sliding_energy;
    double *gpu_dissipated_twisting_energy;

    double *gpu_kinetic_energy;
    double *gpu_rotation_energy;

    double *gpu_normal_energy;
    double *gpu_rolling_energy;
    double *gpu_sliding_energy;
    double *gpu_twisting_energy;


public:
	double E_tot;
    double E_kin;
    double E_rot;
    double V_tot;
    double V_normal;
    double V_rolling;
    double V_sliding;
    double V_twisting;

private:
#endif

    double E_kin_last;
    double V_tot_last;

	// grid specifications
	unsigned int num_grid_cells;

	// data for spatial grid
    unsigned int* gpu_grid_particle_cell;	// grid hash value for each particle
	unsigned int* gpu_grid_particle_index;	// particle index for each particle
	unsigned int* gpu_cell_start;			// index of start of each cell in sorted list
	unsigned int* gpu_cell_end;				// index of end of cell


	// for every particle
    int* gpu_particle_number_of_contacts;	// stores the number of contacts of a certain particle
    int* gpu_particle_particle_ids;			// stores the ids of the other particles a certain particle is in contact with
    int* gpu_particle_contact_ids;			// stores the id of the contacts of a particle (to look it up in the contact list)
    int* gpu_update_local_contact_list;		// > 0 if local contact list of a specific particle needs to be updated

	// for every contact
    int * gpu_contact_type;
    int2* gpu_contact_particle_ids;			// stores the ids of the two particles of the contact
    double4* gpu_contact_normals_new;		// stores contact normal and distance between particles (in .w component)
    double4* gpu_contact_normals_old;
    double3* gpu_contact_n1_initial;		// initial contact pointer
    double3* gpu_contact_n2_initial;		// initial contact pointer
    double4* gpu_contact_rot1;				// rotation state
    double4* gpu_contact_rot2;				// rotation state
    double4* gpu_contact_n1;				// current n1 (.w component stores twisting displacement)
    double4* gpu_contact_n2;				// current n2 (.w component stores relative velocity)

    int* gpu_free_contacts_list;			// stores ids of all free contacs
    int* gpu_next_free_contact;				// index of the next free element in gpu_free_contacts_list array
    int* gpu_number_of_new_contacts;		// number of contacts that have been created in one update step, additionaly reused as the number of broken contacts in one update step
    int2* gpu_new_contact_particle_ids;		// particles ids of the newly created contacts that have yet to be added to the contact lists

	// true if simulation is prepared to be executed on the gpu
	bool cuda_initialized;

	// true  -> the simulation will be executed on the gpu
	// false -> ""    ""                 ""     on th cpu
	bool use_gpu;

	double damping_factor;

	static const int check_filling_factor_interval = 500;
	int check_filling_factor_counter;

};

#endif
