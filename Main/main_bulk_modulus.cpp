#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <cstring>
#include <vector>
//#include "cutil.h"

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

extern double gravity_modifier;

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	#undef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
#endif

#ifdef TRACK_DISSIPATED_ENERGY
	#undef TRACK_DISSIPATED_ENERGY
#endif

bool doesFileExist(const char *filename) 
{
  struct stat stFileInfo;
  int intStat;

  // attempt to get the file attributes
  intStat = stat(filename, &stFileInfo);
  
  if(intStat == 0)
	  return true;
  else 
	  return false;
}

int main(int argc, char **argv)
{
	srand(time(NULL));

	double timestep;
	double max_sim_time;
	double applied_pressure;

	char sample_filename[200];
	char result_filename[200];
	char material_filename[200];

	double side_wall_rolling_modifier = 1.0;
	double side_wall_sliding_modifier = 1.0;

	if(argc == 7 || argc == 9)
	{
		timestep = atof(argv[1]);
		max_sim_time = atof(argv[2]);
		applied_pressure = atof(argv[3]);
		sprintf(sample_filename, "%s", argv[4]);
		sprintf(result_filename, "%s", argv[5]);
		sprintf(material_filename, "%s", argv[6]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n");
		printf("\n -timestep -max_sim_time -applied_pressure -sample_filename -result_filename -material_filename -side_wall_rolling_modifier -side_wall_sliding_modifier\n");
		return 0;
	}

	if(argc == 9)
	{
		side_wall_rolling_modifier = atof(argv[7]);
		side_wall_sliding_modifier = atof(argv[8]);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Performing aggregate compression using particle simulation core v%s\n", CORE_VERSION);

	// load material file
	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(material_filename);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", material_filename, message);
		return 0;
	}

	// calculate the wall thickness for the desired pressure
	double wall_thickness = 10 * applied_pressure / (density * norm(GRAVITY));

	// setup box
	error_code = SimLib::initDynamicCompressionBox(&sim, sample_filename, false, false, 0, wall_thickness, side_wall_rolling_modifier, side_wall_sliding_modifier);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to setup simulation - %s\n", message);
		return 0;
	}

	// start sim
	error_code = sim.startSimulation(max_sim_time, timestep, 0, 0, NULL, NULL, true);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation - %s\n", message);
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double initial_height = sim.box->height;
	double initial_filling_factor = sim.getBoxFillingFactor();

	while(!sim.stop_simulation)
	{
		sim.update();
	}

	// calculate bulk modulus
	double delta_height = initial_height - sim.box->height;

	double delta_V = sim.box->base * delta_height;

	double bulk_modulus = delta_V / applied_pressure;

	// print result
	FILE *file = fopen(result_filename, "a+");

	if(file)
	{
		fprintf(file, "%g %g %g\n", initial_filling_factor, bulk_modulus, applied_pressure);
		fclose(file);
	}

    return EXIT_SUCCESS;
}
