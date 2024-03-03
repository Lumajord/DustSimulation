#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double ENERGY_UNIT;
extern double gravity_modifier;

int main(int argc, char **argv)
{
	double timestep;
	double wall_speed;
	int sample_index;
	int result_index;
	int material_index;

	double ff_stop_wall[23] = {0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60};
	int steps = 23;

	if(argc == 6)
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		sample_index = 3;
		result_index = 4;
		material_index = 5;
	}
	else
	{
		printf("ERROR: Wrong number of arguments!\nUse: -timestep -wall_speed -sample_filename -result_filename -material_filename\n");
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
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	error_code = SimLib::initCompactionBox(&sim, argv[sample_index], 0, 0.8);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to setup compaction:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compact
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef VERBOSE
	printf("Starting compaction...");
	printf("Initial Filling factor: %g\n", sim.getBoxFillingFactor());
#endif

	for(int next_ff_index = 0; next_ff_index < steps; ++next_ff_index)
	{
#ifdef VERBOSE
		printf("Next Stop FF: %g\n", ff_stop_wall[next_ff_index]);
#endif

		error_code = SimLib::initCompactionBox(&sim, NULL, wall_speed, ff_stop_wall[next_ff_index]);
		error_code = sim.startSimulation(1.0, timestep);

		while(!sim.stop_simulation)
			sim.update();

		sim.removeWalls();

		char filename[200];
		sprintf(filename, "%s_%.2f.dat", argv[result_index], ff_stop_wall[next_ff_index]);

		sim.saveToFile(filename);

#ifdef VERBOSE
		printf("Result written to %s\n\n", filename);
#endif
	}

	return EXIT_SUCCESS;
}
