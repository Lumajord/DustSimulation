#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <cutil.h>

#include "CudaDefines.cuh"
#include "SimulationCuda.h"
#include "SimulationLib.h"

extern double moment_of_inertia;
extern double mass;

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double wall_speed;
	double stop_filling_factor;
	int material_file_index = -1;
	int result_file_index = -1;
	int input_file_index = -1;
	int log_interval = 0;

	double ff_stop_wall[23] = {0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60};
	int steps = 23;


    /*
     * How to use
     * timestep = d_c / 10
     * wall rolling modifier = simulation sample size / real sample size = 1/1000
     * filling_factor_interval = particle diameter
     * stop_filling_factor =< 0.6
     * wall speed in cm/s ??
     */


	if(argc == 6)
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		input_file_index = 3;
		result_file_index = 4;
		material_file_index = 5;
	}
	else if(argc == 7 || argc == 8)
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		stop_filling_factor = atof(argv[3]);
		input_file_index = 4;
		result_file_index = 5;
		material_file_index = 6;

		if(argc == 8)
			log_interval = atoi(argv[7]);
	}
	else
	{
        printf("Wrong number of arguments! Use:\n"
               "-timestep -wall_speed -input_filename -result_filename -material_filename\n"
               "-timestep -wall_speed -stop_filling_factor -input_filename -result_filename -material_filename\n"
               "or\n"
               "-timestep -wall_speed -stop_filling_factor -input_filename -result_filename -material_filename\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// load sim
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

	SimulationCuda sim;

	ErrorCode error_code;
	error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

    error_code = SimLib::initCompressionBox(&sim, argv[input_file_index], true, true, wall_speed, stop_filling_factor, 1.0, 0.001, 0.001);
	
	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up the compression box the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

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

	/////////////////////////////////////////////////////////////////////////////////////////
	// start simulation
	/////////////////////////////////////////////////////////////////////////////////////////

	error_code = sim.startSimulation(1, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation!\n%s\n", message);
		return EXIT_SUCCESS;
	}

	printf("Running simulation with %i particles on GPU...\n", sim.number_of_particles);

	if(argc == 6)
	{
		for(int next_ff_index = 0; next_ff_index < steps; ++next_ff_index)
		{
			printf("Next Stop FF: %g\n", ff_stop_wall[next_ff_index]);

			sim.sim_info.info_storage[0] = ff_stop_wall[next_ff_index];
			error_code = sim.startSimulation(1, timestep);

			while(!sim.stop_simulation)
				sim.update();

			char filename[200];
			sprintf(filename, "%s_%.2f.dat", argv[result_file_index], ff_stop_wall[next_ff_index]);

			sim.copySimDataFromGPU();
			sim.saveToFile(filename);

			printf("Result written to %s\n\n", filename);
		}
	}
	else
	{
#ifdef MEASURE_SPEED
		unsigned int timer;
		cutCreateTimer(&timer);
		cutStartTimer(timer);
#endif

		if(log_interval <= 0)
		{
			while(!sim.stop_simulation)
				sim.update();
		}
		else
		{
			int print_positions_counter = 0;
			int print_positions_file_counter = 0;

			while(!sim.stop_simulation)
			{
				sim.update();

				++print_positions_counter;

				if(print_positions_counter >= log_interval)
				{
					print_positions_counter = 0;
					++print_positions_file_counter;

					char filename[300];
					char buf[30];

					strcpy(filename, "positions/");
					sprintf(buf, "positions_%i.dat", print_positions_file_counter);
					strcat(filename, buf);

					sim.copyPositionsFromGPU();
					SimLib::printPositions(sim, filename);
				}
			}
		}

#ifdef MEASURE_SPEED
		cutStopTimer(timer);
#endif

		printf(" ... finished\n");

		sim.copySimDataFromGPU();
		sim.saveToFile(argv[result_file_index]);

#ifdef GPU_TRACK_NUMBER_OF_CONTACTS
		int broken_contacts, created_contacts;
		cudaMemcpy(&broken_contacts, sim.gpu_contacts_broken, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&created_contacts, sim.gpu_contacts_created, sizeof(int), cudaMemcpyDeviceToHost);
		printf("Broken contacts: %i\nCreated contacts: %i\n", broken_contacts, created_contacts);
#endif

#ifdef MEASURE_SPEED
		printf("Cuda execution time: %f ms\nTime per update step: %f ms\n\n", cutGetTimerValue(timer), cutGetTimerValue(timer) / (float)update_steps);
		cutDeleteTimer(timer);
#endif
	}




#ifdef _WIN32
	system("pause");
#endif

	return EXIT_SUCCESS;
}
