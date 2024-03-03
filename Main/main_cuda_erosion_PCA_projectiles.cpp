#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
//#include <cutil.h>

#include "CudaDefines.cuh"
#include "SimulationCuda.h"
#include "SimulationLib.h"

extern double delta_0;

double atanh (double x)
{
   return (log(1.0+x) - log(1.0-x))/2.0;
}

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double impact_speed;

	int number_of_projectiles;
	int projectile_min_size;
	int projectile_max_size;
	double injection_min_distance;
	double injection_max_distance;

	int sample_file_index;
	int result_file_index;
	int material_file_index;

	int seed;
	int GPU_id = 0;

	unsigned int snapshot_interval = 0;
	int snapshot_path_index;

	if(argc == 14)
	{
		timestep = atof(argv[1]);
		impact_speed = atof(argv[2]);
		number_of_projectiles = atoi(argv[3]);
		projectile_min_size = atoi(argv[4]);
		projectile_max_size = atoi(argv[5]);
		injection_min_distance = atof(argv[6]);
		injection_max_distance = atof(argv[7]);
		sample_file_index = 8;
		result_file_index = 9;
		material_file_index = 10;
		seed = atoi(argv[11]);
		snapshot_interval = atoi(argv[12]);
		snapshot_path_index = 13;
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -impact_speed -num_projectiles -min_projectile_size -max_projectile_size -min_injection_distance -max_injection_distance -sample_filename -result_filename -material_filename -seed\n");
		return EXIT_SUCCESS;
	}

	srand(seed);

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup sim sim
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

	SimulationCuda sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	double outer_radius1, outer_radius2;
	error_code = sim.loadFromFile(argv[sample_file_index]);

	if(error_code != EC_OK)
		return error_code;

	SimLib::centerCMS(&sim);
	SimLib::getSize(sim, NULL, &outer_radius1);

	SimLib::rotateSimRandomly(&sim);

	Simulation sim2;
	memcpy(&(sim2.sim_info), &(sim.sim_info), sizeof(SimInfo));

	/////////////////////////////////////////////////////////////////////////////////////////////
	// add projectiles
	/////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> projectile_sizes(number_of_projectiles);
	int number_of_added_particles = 0;

	// determine size of individual projectiles in order to know how many particles will have to be added to the sim
	for(int projectile = 0; projectile < number_of_projectiles; ++projectile)
	{
		double x = (double)rand() / ((double)RAND_MAX+1.0);
		//int projectile_size = (int) exp( log(x * (-pow(projectile_min_size,0.2)+pow(projectile_max_size,0.2))+pow(projectile_min_size,0.2)) * 5.0);
		int projectile_size = (int)(  exp(log( (double)(projectile_max_size-projectile_min_size) ) * x ) * (double)projectile_min_size);

		projectile_sizes[projectile] = projectile_size;
		number_of_added_particles += projectile_size;
	}

	int target_size = sim.number_of_particles;
	int next_particle_id = sim.number_of_particles;
	sim.addParticles(number_of_added_particles);

	for(int projectile = 0; projectile < number_of_projectiles; ++projectile)
	{
		sim2.cleanUp();
		SimLib::initBAMAggregate(&sim2, NULL, projectile_sizes[projectile], 0.0, 0.0, BAM_SELECT_RANDOM);
		SimLib::centerCMS(&sim2);
		SimLib::getSize(sim2, NULL, &outer_radius2);

		vec3 dir;
		vec3 pos, particle_pos;
		vec3 n;

		bool projectile_added = false;

		while(!projectile_added)
		{
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// find trajectory that hits the target
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// determine random direction
			double phi = 2.0 * M_PI * ( (double)rand() / ((double)RAND_MAX));
			double r = 0.15 * ( (double)rand() / ((double)RAND_MAX));

			dir[0] = 1.0;
			dir[1] = r * sin(phi);
			dir[2] = r * cos(phi);
			normalize(&dir);

			SimLib::getOrthoVector(&n, dir);

			phi = 2.0 * M_PI * ( (double)rand() / ((double)RAND_MAX));
			
			double x = (double)rand() / ((double)RAND_MAX+1.0);
			r = 4.0 * outer_radius1* exp( log(x * pow(1.0,1.3)) / 1.3);

			// determine distance
			double dist = 1e-4 * injection_min_distance;

			x = (double)rand() / ((double)RAND_MAX+1.0);
			dist += 1e-4 * (injection_max_distance - injection_min_distance) * (0.5 - atanh(1 - 2.0 * x)  / 6.0 );

			/*if(x < 0.1)
				dist += 1e-4 * (injection_max_distance - injection_min_distance) * exp(log( x*pow(0.2,2) )/2.0);
			else if(x < 0.95)
				dist += 1e-4 * (injection_max_distance - injection_min_distance) * (0.2 + 0.7/0.85 * (x-0.1));
			else
				dist += 1e-4 * (injection_max_distance - injection_min_distance) * ( 0.9 + (x -0.95) * 2.0);*/

			pos[0] = -(outer_radius1+outer_radius2) - dist;
			pos[1] = r * sin(phi);
			pos[2] = r * cos(phi);

			bool can_add_projectile = true;

			for(int p = 0; p < sim2.number_of_particles; ++p)
			{
				particle_pos[0] = pos[0] + sim2.pos_old[X_COORD(p)];
				particle_pos[1] = pos[1] + sim2.pos_old[Y_COORD(p)];
				particle_pos[2] = pos[2] + sim2.pos_old[Z_COORD(p)];

				if(!sim.grid.canAddParticleAt(particle_pos, sim.pos_old))
				{
					can_add_projectile = false;
					break;
				}
			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// add aggregate
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			if(can_add_projectile)
			{
				int id_offset = next_particle_id;

				for(int p = 0; p < sim2.number_of_particles; ++p)
				{
					// add particles
					particle_pos[0] = pos[0] + sim2.pos_old[X_COORD(p)];
					particle_pos[1] = pos[1] + sim2.pos_old[Y_COORD(p)];
					particle_pos[2] = pos[2] + sim2.pos_old[Z_COORD(p)];

					sim.grid.addParticle(particle_pos, next_particle_id);
					sim.pos_old[X_COORD(next_particle_id)] = particle_pos[0];
					sim.pos_old[Y_COORD(next_particle_id)] = particle_pos[1];
					sim.pos_old[Z_COORD(next_particle_id)] = particle_pos[2];
					sim.vel[X_COORD(next_particle_id)] = impact_speed * dir[0];
					sim.vel[Y_COORD(next_particle_id)] = impact_speed * dir[1];
					sim.vel[Z_COORD(next_particle_id)] = impact_speed * dir[2];

					// add contacts
					sim.contact_list[next_particle_id] = NULL;

					ContactListEntry *new_cl_entry = NULL;
					ContactListEntry *cl_entry = sim2.contact_list[p];

					while(cl_entry)
					{
						Contact *new_contact = new Contact();
						*new_contact = *(cl_entry->contact);

						new_contact->id1 += id_offset;
						new_contact->id2 += id_offset;

						if(sim.contact_list[p + id_offset] == NULL) // particle has no other contacts yet
						{
							// create contact list entry
							sim.contact_list[p + id_offset] = new ContactListEntry;
							sim.contact_list[p + id_offset]->next = NULL;
							sim.contact_list[p + id_offset]->id = cl_entry->id + id_offset;
							sim.contact_list[p + id_offset]->contact = new_contact;
							new_cl_entry = sim.contact_list[p + id_offset];
						}
						else
						{
							new_cl_entry->next = new ContactListEntry;
							new_cl_entry->next->next = NULL;
							new_cl_entry->next->id = cl_entry->id + id_offset;
							new_cl_entry->next->contact = new_contact;
							new_cl_entry = new_cl_entry->next;
						}

						cl_entry = cl_entry->next;
					}

					++next_particle_id;
				}

				projectile_added = true;
			}
		}
	}

	sim.initSimData(SIM_TYPE_COLLISION, impact_speed, 0, 0, 0, 0);

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup cuda
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Setting up CUDA...\n");
	fflush(stdout);
	error_code = sim.initCuda(GPU_id);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up CUDA the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	error_code = sim.toggleGPUMode(true);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
		return EXIT_SUCCESS;
	}

	double sim_time = (1e-4 * 1.5 * injection_max_distance) / impact_speed + 2e-4;

	error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation!\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////

	if(snapshot_interval > 0)
	{
		char filename[300];
		char buf[30];
		strcpy(filename, argv[snapshot_path_index]);
		sprintf(buf, "positions_0.dat");
		strcat(filename, buf);

		sim.copyPositionsFromGPU();
		SimLib::printPositions(sim, filename);

		printf("Total number of frames: %i\n", (int) ( (sim_time / timestep) / (double)snapshot_interval) );
	}

	printf("Running simulation on GPU - Simulation time: %g µs...", sim_time * 1e6);
	fflush(stdout);

	unsigned int snapshot_counter = 0;
	unsigned int snapshot_id_counter = 0;

	while(!sim.stop_simulation)
	{
		sim.update();

		++snapshot_counter;

		if(snapshot_interval > 0 && snapshot_counter >= snapshot_interval)
		{
			snapshot_counter = 0;
			++snapshot_id_counter;

			char filename[300];
			char buf[30];

			strcpy(filename, argv[snapshot_path_index]);
			sprintf(buf, "positions_%i.dat", snapshot_id_counter);
			strcat(filename, buf);

			sim.copyPositionsFromGPU();
			SimLib::printPositions(sim, filename);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// finalize
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("...done. Saving result to %s.\n", argv[result_file_index]);
	fflush(stdout);

	sim.copySimDataFromGPU();
	sim.saveToFile(argv[result_file_index], true);

	return EXIT_SUCCESS;
}