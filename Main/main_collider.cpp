#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SimulationLib.h"

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double sim_time;
	double collision_speed;
	double impact_parameter;

	int material_file_index;
	int result_file_index;
	int input_file1_index;
	int input_file2_index;

	if(argc == 9 || argc == 10)
	{
		timestep = atof(argv[1]);
		sim_time = atof(argv[2]);
		collision_speed = atof(argv[3]);
		impact_parameter = atof(argv[4]);
		input_file1_index = 5;
		input_file2_index = 6;
		result_file_index = 7;
		material_file_index = 8;
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -sim_time -collision_speed -impact_parameter -filename1 -filename2 -result_filename -material\nor\n-timestep -sim_time -collision_speed -impact_parameter -filename1 -filename2 -result_filename -material -seed\n");
		return EXIT_SUCCESS;
	}

	if(argc == 10)
		srand(atoi(argv[9]));
	else
	{
		srand(time(NULL));
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup sim sim
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	error_code = SimLib::collideTwoAgglomerates(&sim, argv[input_file1_index], argv[input_file2_index], collision_speed, impact_parameter, 0.0, true, true, true, true);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up the collision the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// start simulation
	/////////////////////////////////////////////////////////////////////////////////////////

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

	while(!sim.stop_simulation)
		sim.update();

	/////////////////////////////////////////////////////////////////////////////////////////
	// finalize
	/////////////////////////////////////////////////////////////////////////////////////////

	sim.saveToFile(argv[result_file_index]);

	return EXIT_SUCCESS;
}