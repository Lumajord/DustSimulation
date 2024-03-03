#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <cstring>
#include <vector>

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
	double wall_speed;

	double start_filling_factor;
	double stop_filling_factor;
	double filling_factor_interval;
	double filling_factor_intervals;
	double filling_factor_profile_interval = 0.05;
	double relaxation_time;

	char sample_filename[200];
	char result_filename[200];
	char material_filename[200];

	double side_wall_rolling_modifier = 1.0;
	double side_wall_sliding_modifier = 1.0;


	if(argc == 10 || argc == 12)
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		start_filling_factor = atof(argv[3]);
		filling_factor_interval = atof(argv[4]);
		filling_factor_intervals = atoi(argv[5]);
		relaxation_time = atof(argv[6]);
		sprintf(sample_filename, "%s", argv[7]);
		sprintf(result_filename, "%s", argv[8]);
		sprintf(material_filename, "%s", argv[9]);
	}
	else if(argc == 9 || argc == 11)
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		start_filling_factor = atof(argv[3]);
		stop_filling_factor = atof(argv[4]);
		filling_factor_profile_interval = atof(argv[5]);
		sprintf(sample_filename, "%s", argv[6]);
		sprintf(result_filename, "%s", argv[7]);
		sprintf(material_filename, "%s", argv[8]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n");
		printf("\n -timestep -wall_speed -start_filling_factor -filling_factor_interval -filling_factor_intervals -relaxation_time -sample_filename -result_filename -material_filename\n");
		printf("\n -timestep -wall_speed -start_filling_factor -filling_factor_interval -filling_factor_intervals -relaxation_time -sample_filename -result_filename -material_filename -side_wall_rolling -side_wall_sliding\n");
		return 0;
	}

	if(argc == 12)
	{
		side_wall_rolling_modifier = atof(argv[10]);
		side_wall_sliding_modifier = atof(argv[11]);
	}

	if(argc == 11)
	{
		side_wall_rolling_modifier = atof(argv[9]);
		side_wall_sliding_modifier = atof(argv[10]);
	}

	// determine filling factors where sample is stored
	std::vector<double> log_filling_factors(filling_factor_intervals);

	for(int i = 0; i < filling_factor_intervals; ++i)
		log_filling_factors[i] = start_filling_factor + (double)i * filling_factor_interval;

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




	// setup box
	if(argc == 9 || argc == 11)
		error_code = SimLib::initCompressionBox(&sim, sample_filename, false, false, true, wall_speed, stop_filling_factor, 1.0, side_wall_rolling_modifier, side_wall_sliding_modifier);
	else
		error_code = SimLib::initCompressionBox(&sim, sample_filename, false, false, true, wall_speed, start_filling_factor + (double)filling_factor_intervals * filling_factor_interval, 1.0, side_wall_rolling_modifier, side_wall_sliding_modifier);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to setup simulation - %s\n", message);
		return 0;
	}

	// start sim
	double sim_time = sim.box->height / wall_speed;
	error_code = sim.startSimulation(sim_time, timestep);

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

	double filling_factor;
	int log_interval = 10;
	int log_counter = 1;
	int filling_factor_counter = 0;

	double next_profile_filling_factor = start_filling_factor;

	while(!sim.stop_simulation)
	{
		sim.update();

		++log_counter;

		if(log_counter >= log_interval)
		{
			log_counter = 1;

			filling_factor = sim.getBoxFillingFactor();

			// density profile
			if(argc == 9 || argc == 11)
			{
				if(filling_factor > next_profile_filling_factor)
				{
					next_profile_filling_factor += filling_factor_profile_interval;

					char filename[300];
					sprintf(filename, "%s_profile_%.0lf_%.2lf.dat", result_filename, wall_speed, filling_factor);
					SimLib::printFillingFactorProfile(sim, filename, 2.0 * particle_radius);
				}
			}
			else
			{
				// store file
				for(int i = 0; i < filling_factor_intervals; ++i)
				{
					if(i >= filling_factor_counter && filling_factor > log_filling_factors[i])
					{
						++filling_factor_counter;

						// determine filename
						char filename[300];
						sprintf(filename, "%s_%.2lf.dat", result_filename, filling_factor);

						sim.saveToFile(filename, true);

						Simulation save_sim;
						save_sim.loadFromFile(filename);
						save_sim.startSimulation(relaxation_time, timestep);

						while(!save_sim.stop_simulation)
							save_sim.update();

						save_sim.removeWalls();
						save_sim.saveToFile(filename, true);
					}
				}
			}
		}
	}

    return EXIT_SUCCESS;
}
