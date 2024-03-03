#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cutil.h>

#include "CudaDefines.cuh"
#include "SimulationCuda.h"
#include "SimulationLib.h"

//#define MEASURE_SPEED

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double sim_time;
	int material_file_index = -1;
	int result_file_index = -1;
	int input_file_index = -1;

	if(argc == 6)
	{
		timestep = atof(argv[1]);
		sim_time = atof(argv[2]);
		input_file_index = 3;
		result_file_index = 4;
		material_file_index = 5;
	}
	else
	{
		printf("Wrong number of arguments! Use:\n-timestep -sim_time -input_filename -result_filename -material_filename\n");
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

	error_code = sim.loadFromFile(argv[input_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nDuring simulation setup the following error occurred:\n%s\n", message);
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

	error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation!\n%s\n", message);
		return EXIT_SUCCESS;
	}

#ifdef MEASURE_SPEED
	unsigned int timer;
#endif

	printf("Running simulation with %i particles on GPU...\n", sim.number_of_particles);

#ifdef MEASURE_SPEED
	cutCreateTimer(&timer);
	cutStartTimer(timer);
#endif
	while(!sim.stop_simulation)
			sim.update();

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

	return EXIT_SUCCESS;
}