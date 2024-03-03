#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
//#include <cutil.h>

#include "CudaDefines.cuh"
#include "SimulationCuda.h"
#include "SimulationLib.h"

extern double particle_radius;

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double impact_speed;
	double impact_distance;
	double impact_angle;
	double injection_edge_distance = 10;

	int projectile_count;
	int projectile_size;

	int target_file_index;
	int log_file_index;
	int result_file_index;
	int material_file_index;
	
	unsigned int snapshot_interval = 0;
	int snapshot_path_index;

	int GPU_id = 0;
	int seed;

	if(argc == 12 || argc == 14)
	{
		timestep = atof(argv[1]);
		impact_speed = atof(argv[2]);
		impact_distance = atof(argv[3]);
		impact_angle = atof(argv[4]);
		projectile_count = atoi(argv[5]);
		projectile_size = atoi(argv[6]);
		target_file_index = 7;
		result_file_index = 8;
		material_file_index = 9;
		log_file_index = 10;
		seed = atoi(argv[11]);
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -impact_speed -impact_distance (µm) -impact_angle -projectile_count -projectile_size -target_filename -result_filename -material_filename -log_filename -seed\n");
		return EXIT_SUCCESS;
	}

	if(argc == 14)
	{
		snapshot_interval = atoi(argv[12]);
		snapshot_path_index = 13;
	}

	srand(seed);

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup sim
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

	error_code = sim.loadFromFile(argv[target_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile loading the target cake the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	SimLib::centerCMS(&sim);
	SimLib::initBox(&sim, NULL, false, 10.0, false, 0, 0);

	// determine size
	vec3 lower_box, upper_box, lower_target, upper_target;
	sim.getEnclosingBox(&lower_box, &upper_box);

	lower_box[0] += particle_radius;
	lower_box[2] += particle_radius;
	upper_box[0] -= particle_radius;
	upper_box[2] -= particle_radius;

	// dont seed particles to close to the edges
	lower_target[0] = lower_box[0] + 1e-4 * injection_edge_distance;
	lower_target[1] = lower_box[1];
	lower_target[2] = lower_box[2] + 1e-4 * injection_edge_distance;
	upper_target[0] = upper_box[0] - 1e-4 * injection_edge_distance;
	upper_target[1] = upper_box[1];
	upper_target[2] = upper_box[2] - 1e-4 * injection_edge_distance;

	/////////////////////////////////////////////////////////////////////////////////////////////
	// add projectiles
	/////////////////////////////////////////////////////////////////////////////////////////////

	int target_size = sim.number_of_particles;
	int next_particle_id = sim.number_of_particles;

	int number_of_added_particles = projectile_count * projectile_size;
	sim.addParticles(number_of_added_particles);

	Simulation sim2;
	memcpy(&(sim2.sim_info), &(sim.sim_info), sizeof(SimInfo));

	vec3 dir = {sin(impact_angle * (M_PI / 180.0)), - cos(impact_angle * (M_PI / 180.0)), 0.0};

	for(int projectile = 0; projectile < projectile_count; ++projectile)
	{
		sim2.cleanUp();
		SimLib::initBAMAggregate(&sim2, NULL, projectile_size, 0.0, 0.0, BAM_SELECT_RANDOM);
		SimLib::centerCMS(&sim2);

		bool projectile_added = false;

		while(!projectile_added)
		{
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// determine injection pos
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			vec3 pos, particle_pos;

			pos[0] = lower_target[0] + (upper_target[0] - lower_target[0]) * (double)rand() / ((double)RAND_MAX+1.0) - dir[0] * 1e-4 * impact_distance;
			pos[1] = upper_target[1] - dir[1] * 1e-4 * impact_distance;
			pos[2] = lower_target[2] + (upper_target[2] - lower_target[2]) * (double)rand() / ((double)RAND_MAX+1.0);

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

	sim.saveToFile("debug_setup.dat", true);

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

	double sim_time = 1.2e-4 + 1e-4 * impact_distance / impact_speed;
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

	printf("Running simulation on GPU - Simulation time: %g µs...", sim_time * 1e6);
	fflush(stdout);

	unsigned int snapshot_counter = 0;
	unsigned int snapshot_id_counter = 0;
	unsigned int log_counter = 0;

	while(!sim.stop_simulation)
	{
		sim.update();

		++snapshot_counter;
		++log_counter;

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

	sim.removeWalls();

	std::vector<int> fragment_ids;							// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;						// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment

	SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

	FILE *file = fopen(argv[log_file_index], "w+");

	if(file)
	{
		fprintf(file, "# Number of particles / fragments: %i / %i\n# Impact velocity / angle: %g / %g\n", sim.number_of_particles, size_of_fragment.size(), impact_speed, impact_angle);
		fprintf(file, "# fragment size   speed   vel_x   vel_y   vel_z\n");
	}

	for(int fragment = 0; fragment < size_of_fragment.size(); ++fragment)
	{
		if(size_of_fragment[fragment] < 1000)
		{
			// determine velocity/angle
			vec3 v = {0,0,0};
			int particles = 0;

			for(std::list<int>::iterator p = particles_of_fragment[fragment].begin(); p != particles_of_fragment[fragment].end(); ++p)
			{
				++particles;
				v[0] += sim.vel[X_COORD(*p)];
				v[1] += sim.vel[Y_COORD(*p)];
				v[2] += sim.vel[Z_COORD(*p)];
			}

			if(particles > 0)
			{
				v[0] /= (double)particles;
				v[1] /= (double)particles;
				v[2] /= (double)particles;

				double velocity = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

				if(file)
					fprintf(file, "%i %g %g %g %g\n", size_of_fragment[fragment], velocity, v[0] / velocity, v[1] / velocity, v[2] / velocity);
			}
		}
	}

	fclose(file);

	return EXIT_SUCCESS;
}