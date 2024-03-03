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

	/*unsigned int number_of_projectile_particles;
	unsigned int number_of_cake_particles;
	double cake_edge_length;
	double cake_top_slice_factor;
	double migration_rate;
	int seed*/

	int projectile_file_index;
	int target_file_index;
	int log_file_index;
	int result_file_index;
	int material_file_index;
	
	unsigned int log_interval = 0;
	unsigned int snapshot_interval = 0;
	int snapshot_path_index;

	if(argc == 10 || argc == 12)
	{
		timestep = atof(argv[1]);
		impact_speed = atof(argv[2]);
		impact_distance = atof(argv[3]);
		projectile_file_index = 4;
		target_file_index = 5;
		result_file_index = 6;
		material_file_index = 7;
		log_interval = atoi(argv[8]);
		log_file_index = 9;
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -impact_speed -impact_distance (µm) -projectile_filename -target_filename -result_filename -material_filename -log_interval -log_file\n");
		return EXIT_SUCCESS;
	}

	if(argc == 12)
	{
		snapshot_interval = atoi(argv[10]);
		snapshot_path_index = 11;
	}

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

	sim.loadFromFile(argv[target_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile loading the target the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	Simulation sim2;
	error_code = sim2.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	sim2.loadFromFile(argv[projectile_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile loading the target the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	unsigned int first_projectile_id = sim.number_of_particles;
	unsigned int last_projectile_id = sim.number_of_particles + sim2.number_of_particles - 1;

	// determine size
	vec3 lower, upper;
	sim.getEnclosingBox(&lower, &upper);

	SimLib::centerCMS(&sim2);
	double outer_radius;
	SimLib::getSize(sim2, NULL, &outer_radius);

	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
		sim2.pos_old[Y_COORD(p)] += outer_radius + upper[1] + impact_distance * 1e-4;

		sim2.vel[X_COORD(p)] = 0;
		sim2.vel[Y_COORD(p)] = - impact_speed;
		sim2.vel[Z_COORD(p)] = 0;
	}

	sim.addParticlesFromSim(sim2);

	SimLib::initBox(&sim, NULL, true, 1.0, true, 0, 0);
	sim.walls[sim.box->top_wall_id].pos[1] += 10.0 * outer_radius;

	sim.initSimData(SIM_TYPE_COLLISION, impact_speed, 0, 1.0, 0, 0);

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup cuda
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Setting up CUDA...\n");
	error_code = sim.initCuda();

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

	double sim_time = (1e-4 * impact_distance + 1.5 * outer_radius ) / impact_speed + 1.3e-4;

	error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation!\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// prepare logging
	/////////////////////////////////////////////////////////////////////////////////////////

	sim.saveToFile("DebugSetup.dat");

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

	FILE *log_file = fopen(argv[log_file_index], "w+");

	if(log_file)
	{
		fprintf(log_file, "# Target / projectile particles: %u / %u, Projectile diameter: %g, Impact speed: %g\n# time  cms x / y / z\n", first_projectile_id, last_projectile_id - first_projectile_id + 1, outer_radius, impact_speed);
		fclose(log_file);
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

		if(log_counter >= log_interval)
		{
			log_counter = 0;

			log_file = fopen(argv[log_file_index], "a");

			if(log_file)
			{
				sim.copyPositionsFromGPU();

				vec3 cms;
				SimLib::getCenterOfMass(&cms, sim, first_projectile_id, last_projectile_id);

				fprintf(log_file, "%g %g %g %g\n", sim.current_time, cms[0], cms[1], cms[2]);
				fclose(log_file);
			}
		}

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