#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
//#include <cutil.h>

#include "CudaDefines.cuh"
#include "SimulationCuda.h"
#include "SimulationLib.h"

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double impact_speed;

	unsigned int number_of_projectiles;
	double injection_edge_distance;
	double injection_height;

	int sample_file_index;
	int result_file_index;
	int material_file_index;
	
	int seed;
	int GPU_id = 0;

	unsigned int snapshot_interval = 0;
	int snapshot_path_index;

	bool target_cake = false;

	if(argc == 9 || argc == 11)
	{
		target_cake = false;

		timestep = atof(argv[1]);
		impact_speed = atof(argv[2]);
		number_of_projectiles = atoi(argv[3]);
		injection_height = atof(argv[4]);
		sample_file_index = 5;
		result_file_index = 6;
		material_file_index = 7;
		seed = atoi(argv[8]);
	}
	else if(argc == 10 || argc == 12)
	{
		timestep = atof(argv[1]);
		impact_speed = atof(argv[2]);
		number_of_projectiles = atoi(argv[3]);
		injection_edge_distance = atof(argv[4]);
		injection_height = atof(argv[5]);
		sample_file_index = 6;
		result_file_index = 7;
		material_file_index = 8;
		seed = atoi(argv[9]);
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -impact_speed -num_projectiles -injection_edge_distance (µm) -injection_height (µm) -sample_filename -result_filename -material_filename -seed\n");
		return EXIT_SUCCESS;
	}

	if(argc == 11)
	{
		snapshot_interval = atoi(argv[9]);
		snapshot_path_index = 10;
	}

	if(argc == 12)
	{
		snapshot_interval = atoi(argv[10]);
		snapshot_path_index = 11;
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

	error_code = sim.loadFromFile(argv[sample_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up the collision the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	SimLib::filterFragments(&sim);
	SimLib::centerCMS(&sim);

	double outer_radius;
	SimLib::getSize(sim, NULL, &outer_radius);

	// shoot projectiles on cake
	if(target_cake)
	{
		SimLib::centerCMS(&sim);

		int next_particle_id = sim.number_of_particles;
		sim.addParticles(number_of_projectiles);

		// determine size
		vec3 lower, upper;
		sim.getEnclosingBox(&lower, &upper);

		// dont seed particles to close to the edges
		lower[0] += 1e-4 * injection_edge_distance;
		lower[2] += 1e-4 * injection_edge_distance;
		upper[0] -= 1e-4 * injection_edge_distance;
		upper[2] -= 1e-4 * injection_edge_distance;

		while(next_particle_id < sim.number_of_particles)
		{
			vec3 pos;
			pos[0] = lower[0] + (double)rand() / ((double)RAND_MAX+1.0) * (upper[0] - lower[0]);
			pos[1] = upper[1] + 1e-4 * injection_height * (double)rand() / ((double)RAND_MAX+1.0);
			pos[2] = lower[2] + (double)rand() / ((double)RAND_MAX+1.0) * (upper[2] - lower[2]);

			if(sim.grid.canAddParticleAt(pos, sim.pos_old))
			{
				sim.grid.addParticle(pos, next_particle_id);
				sim.pos_old[X_COORD(next_particle_id)] = pos[0];
				sim.pos_old[Y_COORD(next_particle_id)] = pos[1];
				sim.pos_old[Z_COORD(next_particle_id)] = pos[2];
				sim.vel[X_COORD(next_particle_id)] = 0;
				sim.vel[Y_COORD(next_particle_id)] = - impact_speed;
				sim.vel[Z_COORD(next_particle_id)] = 0;
				++next_particle_id;
			}
		}

		// init box and move top wall up
		SimLib::initBox(&sim, NULL, true, 1.0, true, 0, 0);
		sim.walls[sim.box->top_wall_id].pos[1] += 1e-3 * injection_height;
	}
	// shoot particles on aggregate from random directions
	else
	{
		ErrorCode error_code = SimLib::sandblastAggregate(&sim, NULL, number_of_projectiles, impact_speed, injection_height);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("Adding projectiles failed: %s\n", message);
			return EXIT_SUCCESS;
		}
	}

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

	double sim_time;

	if(target_cake)
		sim_time = 1e-4 * injection_height / impact_speed + 1.2e-4;
	else
		sim_time = (1e-4 * injection_height + outer_radius) / impact_speed + 1.5e-4;

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

	int initial_particles = sim.number_of_particles - number_of_projectiles;
	sim.removeWalls();

	std::vector<int> fragment_ids;							// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;						// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment

	SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

	if(target_cake)
	{
		int max_fragment_size = 0;

		// determine biggest agglomerate
		for(unsigned int agg = 0; agg < size_of_fragment.size(); ++agg)
		{
			if(size_of_fragment[agg] > max_fragment_size)
				max_fragment_size = size_of_fragment[agg];
		}

		double erosion_efficiency = (double)(initial_particles - max_fragment_size) / (double)(number_of_projectiles);

		FILE *file = fopen("analysis.txt", "a+");
		fprintf(file, "%i %g\n", (int) impact_speed, erosion_efficiency);
		fclose(file);
	}
	else
	{
		int total_particles = 0;
		int max_fragment_size = 0;

		// fragments with more than 1000 monomers will not be considered fragments
		for(unsigned int agg = 0; agg < size_of_fragment.size(); ++agg)
		{
			if(size_of_fragment[agg] > max_fragment_size)
				max_fragment_size = size_of_fragment[agg];

			if(size_of_fragment[agg] > 9)
				total_particles += size_of_fragment[agg];
		}

		double erosion_efficiency_max = (double)(initial_particles - max_fragment_size) / (double)(number_of_projectiles);
		double erosion_efficiency = (double)(initial_particles - total_particles) / (double)(number_of_projectiles);

		FILE *file = fopen("analysis.txt", "a+");
		fprintf(file, "%i %g %g\n", (int) impact_speed, erosion_efficiency_max, erosion_efficiency);
		fclose(file);
	}

	return EXIT_SUCCESS;
}