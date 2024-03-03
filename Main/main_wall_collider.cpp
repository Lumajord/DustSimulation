#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double ENERGY_UNIT;

extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;

extern double particle_radius;
extern double density;
extern double surface_energy;
extern double nu;
extern double young_mod;
extern double crit_rolling_displacement;
extern double yield_strength;
extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	#undef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
#endif

#ifdef TRACK_DISSIPATED_ENERGY
	#undef TRACK_DISSIPATED_ENERGY
#endif

int main(int argc, char **argv)
{
	srand(time(NULL));

	// for main simulation
	double timestep;
	double impact_velocity;
	double impact_angle;

	int index_aggregate_filename;
	int index_result_filename;
	int index_material_filename;

	if(argc == 7)
	{
		timestep = atof(argv[1]);
		impact_velocity = atof(argv[2]);
		impact_angle = atoi(argv[3]);
		index_aggregate_filename = 4;
		index_result_filename = 5;
		index_material_filename = 6;
	}
	else
	{
		printf("Incorrect arguments! Use:\n -timestep -impact_velocity -impact_angle -aggregate_filename -result_filename -material_filename\n");
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[index_material_filename]);

	if(error_code != EC_OK)
	{
		char message[300];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %200s, the following error occurred:\n%s\n", argv[index_material_filename], message);
		return EXIT_SUCCESS;
	}

	error_code = sim.loadFromFile(argv[index_aggregate_filename]);

	if(error_code != EC_OK)
	{
		char message[300];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to load %200s:\n%s\n", argv[index_aggregate_filename], message);
		return EXIT_SUCCESS;
	}

	// get outer radius (to calculate sim time)
	double outer_radius;
	SimLib::getSize(sim, NULL, &outer_radius);

	error_code = SimLib::collideAgglomerateWithWall(&sim, NULL, impact_velocity, impact_angle, 0.01, true, 0, 0, 0);

	if(error_code != EC_OK)
	{
		char message[300];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to setup the collision, the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	// determine sim time
	double sim_time = outer_radius / impact_velocity;

	if(sim_time < 1e-5)
		sim_time = 1e-5;

	sim.startSimulation(sim_time, timestep);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	while(!sim.stop_simulation)
		sim.update();

	sim.saveToFile(argv[index_result_filename]);

    return EXIT_SUCCESS;
}
