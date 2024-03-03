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
extern double T_vis;
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

void logPressure(Simulation* sim, const char* log_filename, std::vector<double> &top_forces, std::vector<double> &bottom_forces)
{
	// determine pressure on top&bottom wall
	double force_top = 0;
	double force_bottom = 0;

	for(unsigned int i = 0; i < top_forces.size(); ++i)
	{
		force_top += top_forces[i];
		force_bottom += bottom_forces[i];
	}

	force_top /= (double)top_forces.size();
	force_bottom /= (double)top_forces.size();

	double pressure_top = 0.1 * force_top/sim->box->base; // /10 to convert from CGS to SI
	double pressure_bottom = 0.1 * force_bottom/sim->box->base; // /10 to convert from CGS to SI

	// write to log file
	FILE* log_file = fopen(log_filename, "a");

	if(log_file)
	{	
		fprintf(log_file, "%g %g %g %g\n", sim->current_time, sim->getBoxFillingFactor(), pressure_top, pressure_bottom);
		fclose(log_file);
	}
}

void logKineticEnergy(Simulation* sim, const char* log_filename, std::vector<double> &top_forces, std::vector<double> &bottom_forces)
{
	// determine pressure on top&bottom wall
	double force_top = 0;
	double force_bottom = 0;

	for(unsigned int i = 0; i < top_forces.size(); ++i)
	{
		force_top += top_forces[i];
		force_bottom += bottom_forces[i];
	}

	force_top /= (double)top_forces.size();
	force_bottom /= (double)top_forces.size();

	double pressure_top = 0.1 * force_top/sim->box->base; // /10 to convert from CGS to SI
	double pressure_bottom = 0.1 * force_bottom/sim->box->base; // /10 to convert from CGS to SI

	// write to log file
	FILE* log_file = fopen(log_filename, "a");

	if(log_file)
	{
		fprintf(log_file, "%g %g %g %g\n", sim->current_time, pressure_top, pressure_bottom, SimLib::getKineticEnergy(*sim) * ENERGY_UNIT);
		fclose(log_file);
	}
}

int main(int argc, char **argv)
{
	srand(time(NULL));

	// for main simulation
	double timestep;
	double wall_speed;
	double stop_filling_factor;

	double wall_compression_modifier = 1.0; 
	double wall_rolling_modifier = 1.0;
	double wall_sliding_modifier = 1.0;

	int log_interval;
	int sample_filename_index;
	int log_filename_index;
	int material_filename_index;

	/*
	// after the has moved a certain distance (in particle radii), wall is stopped for a certain amount of time to allow relaxation of the aggregate
	double relaxation_distance;
	double relaxation_time;

	// relaxation stops when the dissipated kinetic energy is below this factor
	double stop_energy_dissipation_factor;
	double initial_kinetic_energy = 0;
	*/

	if(argc == 6) // continue simulation
	{
		timestep = atof(argv[1]);
		sample_filename_index = 2;
		material_filename_index = 3;
		log_filename_index = 4;
		log_interval = atoi(argv[5]);
	}
	else if(argc == 8 || argc == 11)	// main mode -> measure pressure while wall is moving downwards
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		stop_filling_factor = atof(argv[3]);
		sample_filename_index = 4;
		material_filename_index = 5;
		log_filename_index = 6;
		log_interval = atoi(argv[7]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n -timestep -wall_speed -stop_filling_factor -sample_filename -material_filename -log_filename -log_interval\n");
		return EXIT_SUCCESS;
	}

	if(argc == 11)
	{
		wall_compression_modifier = atof(argv[8]);
		wall_rolling_modifier = atof(argv[9]);
		wall_sliding_modifier = atof(argv[10]);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Performing aggregate compression using particle simulation core v%s\n", CORE_VERSION);

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_filename_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_filename_index], message);
		return EXIT_SUCCESS;
	}

	double sim_time; 

	// setup sim
	if(argc == 6) // continue simulation
	{
		error_code = sim.loadFromFile(argv[sample_filename_index]);
		
		sim_time = sim.box->height / wall_speed;
		stop_filling_factor = sim.sim_info.info_storage[0];
		wall_speed = norm(sim.walls[sim.box->top_wall_id].velocity);
	}
	else
	{
		error_code = SimLib::initCompressionBox(&sim, argv[sample_filename_index], true, false, wall_speed, stop_filling_factor, wall_compression_modifier, wall_rolling_modifier, wall_sliding_modifier);
		sim_time = sim.box->height / wall_speed;
	}

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to setup simulation - %s\n", message);
		return EXIT_SUCCESS;
	}

	// start sim
	error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation - %s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare log file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *log_file;
	log_file = fopen(argv[log_filename_index], "w+");
	fprintf(log_file, "# number of particles: %i\n", sim.number_of_particles);
	fprintf(log_file, "# box base size (in m^2): %g\n", 1e-4 * sim.box->base);
	fprintf(log_file, "# box height (in m): %g \n", 0.01 * sim.box->height);
	fprintf(log_file, "# wall speed (in cm/s): %g \n", wall_speed);
	fprintf(log_file, "# initial/stop filling factor: %g %g\n", sim.getBoxFillingFactor(), stop_filling_factor);
	fprintf(log_file, "# T_vis: %g\n", T_vis);


	if(argc == 11)
		fprintf(log_file, "# compression/rolling/sliding modifier: %g / %g / %g\n", wall_compression_modifier, wall_rolling_modifier, wall_sliding_modifier);

	//if(sim.use_gravity)
	//	fprintf(log_file, "# Gravity enabled! - gravity strength modifier: %g\n", gravity_modifier);

	fprintf(log_file, "#\n# time      filling factor     top pressure (in Pa)     bottom pressuren (in Pa)\n");
	fclose(log_file);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int log_counter = 0;
	std::vector<double> top_forces(log_interval, 0);
	std::vector<double> bottom_forces(log_interval, 0);

	while(!sim.stop_simulation)
	{
		sim.update();

		top_forces[log_counter] = sim.walls[sim.box->top_wall_id].getWallForce();
		bottom_forces[log_counter] = sim.walls[sim.box->bottom_wall_id].getWallForce();
		++log_counter;

		// write data to log file
		if(log_counter == log_interval-1)
		{
			log_counter = 0;
			logPressure(&sim, argv[log_filename_index], top_forces, bottom_forces);
		}
	}

    return EXIT_SUCCESS;
}
