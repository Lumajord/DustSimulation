#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cutil.h>

#include "CudaDefines.cuh"
#include "SimulationCuda.h"
#include "SimulationLib.h"

extern double moment_of_inertia;
extern double mass;

bool use_gpu = false;
bool use_cpu = true;

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	int update_steps = 1;
	int log_interval = 0;
	double timestep;
	int material_file_index = -1;
	int result_file_index = -1;
	int input_file_index = -1;

	if(argc == 5)
	{
		timestep = atof(argv[1]);
		input_file_index = 2;
		material_file_index = 3;
		update_steps = atoi(argv[4]);
	}
	else if(argc == 6)
	{
		timestep = atof(argv[1]);
		input_file_index = 2;
		material_file_index = 3;
		update_steps = atoi(argv[4]);
		log_interval = atoi(argv[5]);
	}
	else if(argc == 7)
	{
		timestep = atof(argv[1]);
		input_file_index = 2;
		result_file_index = 3;
		material_file_index = 4;
		update_steps = atoi(argv[5]);
		log_interval = atoi(argv[6]);
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -filename -material -integration_steps -log_interval\nor\n-timestep -filename -result_filename -material -integration_steps -log_interval\n");
		return EXIT_SUCCESS;
	}

	double simulation_time = timestep * (update_steps+1);

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
		return 0;
	}

	error_code = sim.loadFromFile(argv[input_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nDuring simulation setup the following error occurred:\n%s\n", message);
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup cuda
	/////////////////////////////////////////////////////////////////////////////////////////

	if(use_gpu)
	{
		printf("Setting up CUDA...\n");
		error_code = sim.initCuda();

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR:\nWhile setting up CUDA the following error occurred:\n%s\n", message);
			return 0;
		}

		error_code = sim.toggleGPUMode(use_gpu);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
			return 0;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// start simulation
	/////////////////////////////////////////////////////////////////////////////////////////

	error_code = sim.startSimulation(simulation_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation!\n%s\n", message);
		return 0;
	}

#ifdef MEASURE_SPEED
	cudaEvent_t start, stop;
#endif

	if(use_gpu)
	{
		printf("Running simulation with %i particles on GPU...\n", sim.number_of_particles);

#ifdef MEASURE_SPEED
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);
#endif

		if(log_interval <= 0)
		{
			for(int i = 0; i < update_steps; ++i)
			{
				sim.update();

				if(sim.stop_simulation)
					break;
			}
		}
		else
		{
			int print_positions_counter = 0;
			int print_positions_file_counter = 0;

			for(int i = 0; i < update_steps; ++i)
			{
				sim.update();

				if(sim.stop_simulation)
					break;

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
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
#endif
		printf(" ... finished\n");

		if(argc == 7)
		{
			sim.copySimDataFromGPU();
			sim.saveToFile(argv[result_file_index]);
		}

#ifdef GPU_TRACK_NUMBER_OF_CONTACTS
		int broken_contacts, created_contacts;
		cudaMemcpy(&broken_contacts, sim.gpu_contacts_broken, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&created_contacts, sim.gpu_contacts_created, sizeof(int), cudaMemcpyDeviceToHost);
		printf("Broken contacts: %i\nCreated contacts: %i\n", broken_contacts, created_contacts);
#endif

#ifdef MEASURE_SPEED
		float elapsed_time;
		cudaEventElapsedTime(&elapsed_time, start, stop);
		printf("Cuda execution time: %f ms\nTime per update step: %f ms\n\n", elapsed_time, elapsed_time / (float)update_steps);
#endif
	}

	if(use_cpu)
	{
		error_code = sim.toggleGPUMode(false);

		printf("Running simulation on CPU...");
		error_code = sim.loadFromFile(argv[input_file_index]);

		error_code = sim.startSimulation(simulation_time, timestep);

#ifdef MEASURE_SPEED
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);
#endif

		for(int i = 0; i < update_steps; ++i)
			sim.update();

#ifdef MEASURE_SPEED
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
#endif

		printf(" ... finished\n");
		printf("Broken contacts: %i\nCreated contacts: %i\n", sim.broken_contacts, sim.created_contacts);

#ifdef MEASURE_SPEED
		float elapsed_time;
		cudaEventElapsedTime(&elapsed_time, start, stop);
		printf("CPU execution time: %f ms\nTime per update step: %f ms\n", elapsed_time, elapsed_time / (float)update_steps);
#endif
	}

#ifdef _WIN32
	system("pause");
#endif

	return EXIT_SUCCESS;
}