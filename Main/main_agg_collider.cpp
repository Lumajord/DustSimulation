#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SimulationLib.h"

#include <sstream>
#include <string>

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double sim_time;
	double collision_speed;
	double impact_parameter;
    double impact_distance;

	int material_file_index;
    int result_file_index;
	int input_file1_index;
	int input_file2_index;
    int seed = 0;

    if(argc == 11)
	{
		timestep = atof(argv[1]);
		sim_time = atof(argv[2]);
		collision_speed = atof(argv[3]);
		impact_parameter = atof(argv[4]);
        impact_distance = atof(argv[5]);
        input_file1_index = 6;
        input_file2_index = 7;
        result_file_index = 8;
        material_file_index = 9;
        seed = atoi(argv[10]);
	}
	else
	{
        printf("Wrong number of arguments! Use:\n-timestep -sim_time -collision_speed -impact_parameter -filename1 -filename2 -aggregate_result_filename -data_result_filename -material\nor\n-timestep -sim_time -collision_speed -impact_parameter -filename1 -filename2 -result_filename -material -seed\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup sim sim
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

	Simulation sim;
    sim.rand_generator.seed(seed);
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);


	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}


    //error_code = SimLib::collideTwoAgglomerates(&sim, argv[input_file1_index], argv[input_file2_index], collision_speed, impact_parameter, impact_distance, true, true, true, true, seed);
    error_code = sim.loadFromFile(argv[input_file1_index]);


	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up the collision the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}



    char energy_file[1024];
    sprintf(energy_file, "../energiesCPU_");
    strcat(energy_file, argv[result_file_index]+3);

	/////////////////////////////////////////////////////////////////////////////////////////
	// start simulation
	/////////////////////////////////////////////////////////////////////////////////////////
    printf("starting Simulation...\n\n");
    error_code = sim.startSimulation(sim_time, timestep, 100, 0, energy_file);

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
	{


		sim.update();


	}

    /////////////////////////////////////////////////////////////////////////////////////////
    // finalize
    /////////////////////////////////////////////////////////////////////////////////////////

    printf("...done\n");
    fflush(stdout);

    sim.saveToFile(argv[result_file_index], true);

	return EXIT_SUCCESS;
}
