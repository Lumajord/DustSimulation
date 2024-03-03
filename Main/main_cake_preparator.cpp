#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern bool damping_enabled;
extern double damping_factor;

int main(int argc, char **argv)
{
	double timestep;
	double top_slice_factor;
	double bottom_slice_factor;
	double disturbance;
	double damping_time;
	int sample_index;
	int result_index;
	int material_index;

	if(argc == 9)
	{
		timestep = atof(argv[1]);
		top_slice_factor = atof(argv[2]);
		bottom_slice_factor = atof(argv[3]);
		disturbance = atof(argv[4]);
		damping_time = atof(argv[5]);
		sample_index = 6;
		result_index = 7;
		material_index = 8;
	}
	else
	{
		printf("ERROR: Wrong number of arguments!\nUse: -timestep -top_slice_factor -bottom_slice_factor -disturbance -damping_time -sample_filename -result_filename -material_filename\n");
		return EXIT_SUCCESS;
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

	sim.loadFromFile(argv[sample_index]);

	if(top_slice_factor > 0)
		error_code = SimLib::sliceTop(&sim, top_slice_factor);

	if(error_code == EC_OK && bottom_slice_factor > 0)
		error_code = SimLib::sliceBottom(&sim, bottom_slice_factor);

	printf("Finished slicing - new number of particles: %i\n", sim.number_of_particles);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to slice aggregate:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	error_code = SimLib::initCompressionBox(&sim, NULL, false, true, true, 0, 0.75, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0);


	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to init box:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	error_code = SimLib::resetContacts(&sim);

	if(error_code == EC_OK)
		error_code = SimLib::disturbParticles(&sim, disturbance * particle_radius);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to disturb box:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// relax
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Starting relaxation\n");
	error_code = sim.startSimulation(damping_time, timestep, 0, 0, NULL, NULL, true);



	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start relaxation:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	while(!sim.stop_simulation)
		sim.update();

	// damp velocities
	printf("Starting damping of velocities\n");
	SimLib::resetContacts(&sim);
	error_code = sim.startSimulation(8000.0 * timestep, timestep, 0, 0, NULL, NULL, false);
	damping_enabled = true;
	damping_factor = 0.999;

	while(!sim.stop_simulation)
			sim.update();

	printf("Saving result to %s\n", argv[result_index]);
	fflush(stdout);

	sim.saveToFile(argv[result_index]);

    return EXIT_SUCCESS;
}
