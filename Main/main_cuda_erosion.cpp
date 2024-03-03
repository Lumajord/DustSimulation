#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
//#include <cutil.h>

#if !defined(GPU_OPEN_BOX)
	#error GPU_OPEN_BOX must be defined!
#endif

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
	double injection_distance;

	int sample_file_index;
	int result_file_index;
	int material_file_index;
	
	int seed;
	int GPU_id = 0;

	unsigned int snapshot_interval = 0;
	int snapshot_path_index;

	bool target_cake = true;

	if(argc == 9 || argc == 11)
	{
		target_cake = false;

		timestep = atof(argv[1]);
		impact_speed = atof(argv[2]);
		number_of_projectiles = atoi(argv[3]);
		injection_distance = atof(argv[4]);
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
		injection_distance = atof(argv[5]);
		sample_file_index = 6;
		result_file_index = 7;
		material_file_index = 8;
		seed = atoi(argv[9]);
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -impact_speed -num_projectiles -injection_edge_distance (µm) -injection_distance (µm) -sample_filename -result_filename -material_filename -seed\n");
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

		//SimLib::initOpenBox(&sim, NULL, true, 1.0, 1.0, 1.0, 1.0);
		SimLib::initBox(&sim, NULL, true, 10.0, true, 0, 0);

		int next_particle_id = sim.number_of_particles;
		sim.addParticles(number_of_projectiles);

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

		while(next_particle_id < sim.number_of_particles)
		{
			//vec3 dir = {0.0, -1.0, 0.0};

			// determine random direction
			vec3 dir;

			while(true)
			{
				double phi = 2.0 * M_PI * ( (double)rand() / ((double)RAND_MAX));
				double theta = asin( 2.0 * (double)rand() / ((double)RAND_MAX) - 1.0 ) + 0.5 * M_PI;

				dir[0] = sin(theta) * cos(phi);
				dir[1] = sin(theta) * sin(phi);
				dir[2] = cos(theta);

				if(dir[1] > 0)
					dir[1] *= -1.0;

				// determine angle
				if( fabs(dir[1]) > 0.17)
					break;
			}

			//
			double dist = ( (double)rand() / (double)RAND_MAX) * 1e-4 * injection_distance;

			vec3 pos;
			pos[0] = lower_target[0] + (double)rand() / ((double)RAND_MAX+1.0) * (upper_target[0] - lower_target[0]) - dist * dir[0];
			pos[1] = upper_target[1] - dist * dir[1];
			pos[2] = lower_target[2] + (double)rand() / ((double)RAND_MAX+1.0) * (upper_target[2] - lower_target[2]) - dist * dir[2];

			// make sure particles spwan inside of the box
			if(pos[0] < lower_box[0])
			{
				double k = (lower_box[0] - pos[0]) / dir[0];
				pos[0] += dir[0] * k;
				pos[1] += dir[1] * k;
				pos[2] += dir[2] * k;
			}

			if(pos[2] < lower_box[2])
			{
				double k = (lower_box[2] - pos[2]) / dir[2];
				pos[0] += dir[0] * k;
				pos[1] += dir[1] * k;
				pos[2] += dir[2] * k;
			}

			if(pos[0] > upper_box[0])
			{
				double k = (upper_box[0] - pos[0]) / dir[0];
				pos[0] += dir[0] * k;
				pos[1] += dir[1] * k;
				pos[2] += dir[2] * k;
			}	
				
			if(pos[2] > upper_box[2])
			{
				double k = (upper_box[2] - pos[2]) / dir[2];
				pos[0] += dir[0] * k;
				pos[1] += dir[1] * k;
				pos[2] += dir[2] * k;
			}

			if(sim.grid.canAddParticleAt(pos, sim.pos_old))
			{
				sim.grid.addParticle(pos, next_particle_id);
				sim.pos_old[X_COORD(next_particle_id)] = pos[0];
				sim.pos_old[Y_COORD(next_particle_id)] = pos[1];
				sim.pos_old[Z_COORD(next_particle_id)] = pos[2];
				sim.vel[X_COORD(next_particle_id)] = impact_speed * dir[0];
				sim.vel[Y_COORD(next_particle_id)] = impact_speed * dir[1];
				sim.vel[Z_COORD(next_particle_id)] = impact_speed * dir[2];
				++next_particle_id;
			}
		}
	}
	// shoot particles on aggregate from random directions
	else
	{
		ErrorCode error_code = SimLib::sandblastAggregate(&sim, NULL, number_of_projectiles, impact_speed, injection_distance);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("Adding projectiles failed: %s\n", message);
			return EXIT_SUCCESS;
		}
	}

	sim.saveToFile("erosion_setup.dat");

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

	double sim_time = 3e-4;

	if(target_cake)
		sim_time += 1e-4 * injection_distance / impact_speed;
	else
		sim_time += (1e-4 * injection_distance + outer_radius) / impact_speed + 0.3e-4;

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

	printf("...done.\nSaving result to %s.\n", argv[result_file_index]);
	fflush(stdout);

	sim.copySimDataFromGPU();

	printf("Data copied from GPU\n");
	fflush(stdout);

	sim.saveToFile(argv[result_file_index], true);

	printf("Saved to file\n");
	fflush(stdout);

	sim.removeWalls();

	printf("Walls removed\n");
	fflush(stdout);

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