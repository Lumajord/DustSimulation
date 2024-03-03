#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double ENERGY_UNIT;
extern double gravity_modifier;

#ifndef ENABLE_GRAVITY
	#error Gravity not enabled!
#endif

int main(int argc, char **argv)
{
	double duration;
	double timestep;
	double wall_speed_slow;
	double wall_speed_fast;
	double relaxation_time;
	int sample_index;
	int result_index;
	int material_index;

	double ff_slow_wall[21] = {0.158, 0.178, 0.198, 0.218, 0.237, 0.257, 0.277, 0.296, 0.316, 0.335, 0.355, 0.375, 0.395, 0.415, 0.435, 0.455, 0.475, 0.495, 0.515, 0.535, 0.555};
	double ff_stop_wall[21] = {0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56};
	int steps = 21;

	if(argc == 7)
	{
		timestep = atof(argv[1]);
		duration = atof(argv[2]);
		gravity_modifier = atof(argv[3]);
		sample_index = 4;
		result_index = 5;
		material_index = 6;
	}
	else if(argc == 9)
	{
		timestep = atof(argv[1]);
		wall_speed_fast = atof(argv[2]);
		wall_speed_slow = atof(argv[3]);
		relaxation_time = atof(argv[4]);
		gravity_modifier = atof(argv[5]);
		sample_index = 6;
		result_index = 7;
		material_index = 8;
	}
	else
	{
		printf("ERROR: Wrong number of arguments!\nUse: -timestep -duration -gravity_strength_modifier -sample_filename -result_filename -material_filename\n");
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	ErrorCode error_code;
	Simulation sim;

	error_code = sim.loadMaterial(argv[material_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_index], message);
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(argc == 7)
	{
		error_code = SimLib::initCompactionBox(&sim, argv[sample_index], false);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR:\nFailed to setup compaction:\n%s\n", message);
			return EXIT_SUCCESS;
		}

		error_code = sim.startSimulation(duration, timestep, 0, 0, NULL, NULL, true);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR:\nFailed to start compaction:\n%s\n", message);
			return EXIT_SUCCESS;
		}
	}
	else if(argc == 9)
	{
		error_code = SimLib::initCompressionBox(&sim, argv[sample_index], false, true, false, wall_speed_fast, ff_slow_wall[0], 1.0, 0.0, 0.0); 

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR:\nFailed to setup compaction:\n%s\n", message);
			return EXIT_SUCCESS;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compact
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Starting compaction...");

	if(argc == 7)
	{
		while(!sim.stop_simulation)
			sim.update();

		// damp velocities
		sim.startSimulation(1e-6, timestep, 0, 0, NULL, NULL, false);

		while(!sim.stop_simulation)
		{
			sim.update();
			sim.dampVelocities(0.995);
		}

		memset(sim.vel, 0, 3 * sim.number_of_particles * sizeof(double));
		memset(sim.vel_angular, 0, 3 * sim.number_of_particles * sizeof(double));

		sim.removeWalls();

		// save file with contacts in initial state
		sim.saveToFile(argv[5], false);
		sim.loadFromFile(argv[5]);
		memcpy(sim.pos_new, sim.pos_old, sim.number_of_particles * 3 * sizeof(double));
		sim.updateSticking();
		sim.saveToFile(argv[result_index]);
	}
	else if(argc == 9)
	{
		for(int next_ff_index = 0; next_ff_index < steps; ++next_ff_index)
		{
			error_code = SimLib::initCompressionBox(&sim, NULL, false, true, false, wall_speed_fast, ff_slow_wall[next_ff_index], 1.0, 0.0, 0.0);
			error_code = sim.startSimulation(1.0, timestep, 0, 0, NULL, NULL, true);

			while(!sim.stop_simulation)
				sim.update();

			error_code = SimLib::initCompressionBox(&sim, NULL, false, true, false, wall_speed_slow, ff_stop_wall[next_ff_index]-0.001, 1.0, 0.0, 0.0);
			error_code = sim.startSimulation(1.0, timestep, 0, 0, NULL, NULL, false);

			while(!sim.stop_simulation)
				sim.update();

			error_code = SimLib::initCompressionBox(&sim, NULL, false, true, false, 0, ff_stop_wall[next_ff_index], 1.0, 0.0, 0.0);
			error_code = sim.startSimulation(relaxation_time, timestep, 0, 0, NULL, NULL, false);

			while(!sim.stop_simulation)
				sim.update();

			sim.removeWalls();

			char filename[200];
			sprintf(filename, "%s_%.2f.dat", argv[result_index], ff_stop_wall[next_ff_index]);
			sim.saveToFile(filename);

			if(ff_slow_wall[next_ff_index] > 0.35)
				gravity_modifier += 25.0;
		}
	}

    return EXIT_SUCCESS;
}
