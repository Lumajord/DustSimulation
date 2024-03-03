#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double ENERGY_UNIT;
extern double gravity_modifier;
extern double particle_radius;
extern double azimuthal_acceleration;

int main(int argc, char **argv)
{
	double timestep;
	double sim_time;
	int log_interval = 0;

	int sample_filename_index = 4;
	int material_filename_index = 5;
	int log_filename_index = 7;

	if(argc == 8)
	{
		timestep = atof(argv[1]);
		azimuthal_acceleration = atof(argv[2]);
		sim_time = atof(argv[3]);
		log_interval = atoi(argv[6]);
	}
	else
	{
		printf("ERROR: Wrong number of arguments! Use:\n-timestep -azimuthal_acceleration (in 1/s^2) -sim_time -sample_filename -material_filename -log_interval -result_filename\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	ErrorCode error_code;
	Simulation sim;

	error_code = sim.loadMaterial(argv[material_filename_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_filename_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load sim
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	error_code = sim.loadFromFile(argv[sample_filename_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to load simulation - %s\n", message);
		return EXIT_SUCCESS;
	}

	error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation - %s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run sim
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	FILE *log_file = fopen(argv[log_filename_index], "w+");

	if(log_file)
	{
		fprintf(log_file, "# time    omega    actual mean omega    sigma omega\n");
		fclose(log_file);
	}

	int log_counter = 0;

	while(!sim.stop_simulation)
	{
		sim.update();

		++log_counter;

		if(log_interval > 0 && log_counter >= log_interval)
		{
			log_counter = 0;

			double mean, sigma;
			sim.calculateRotation(&mean, &sigma);

			FILE *log_file = fopen(argv[log_filename_index], "a");

			if(log_file)
			{
				fprintf(log_file, "%g %g %g %g\n", sim.current_time, sim.current_time * azimuthal_acceleration, mean, sigma);
				fclose(log_file);
			}
		}
	}

    return EXIT_SUCCESS;
}
