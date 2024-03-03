#include <stdlib.h>
#include <stdio.h>

#include "Simulation.h"
#include "SimulationLib.h"

const char *material_filename = "Silicate.mat";
const char *data_filename = "TestCollision2.dat";
const char *log_filename = "Data_dump.txt";

int main(int argc, char **argv)
{
	Simulation sim;

	ErrorCode error_code = sim.loadMaterial(material_filename);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", material_filename, message);
		return EXIT_SUCCESS;
	}

	error_code = sim.loadFromFile(data_filename);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile loading data the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	double duration = 1e-4;
	double timestep = 3e-10;

	sim.startSimulation(duration, timestep);

	// sim.stop_simulation will be set to true when max sim time is reached
	//while(!sim.stop_simulation)
	for(int i = 0; i < 10; ++i)
	{
		sim.update();
	}

	sim.debugDump(log_filename);

	return EXIT_SUCCESS;
}
