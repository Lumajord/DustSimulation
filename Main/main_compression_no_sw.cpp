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

void logPressure(Simulation* sim, const char* log_filename, std::vector<double> &top_forces, std::vector<double> &top_box_forces, std::vector<double> &bottom_forces, std::vector<double> &bottom_box_forces)
{
	// determine pressure on top&bottom wall
	double force_top = 0;
	double force_bottom = 0;
	double force_box_top = 0;
	double force_box_bottom = 0;
	double max_force_top = 0;

	for(unsigned int i = 0; i < top_forces.size(); ++i)
	{
		force_top += top_forces[i];
		force_box_top += top_box_forces[i];
		force_bottom += bottom_forces[i];
		force_box_bottom += bottom_box_forces[i];

		if(top_forces[i] > max_force_top)
			max_force_top = top_forces[i];
	}

	force_top /= (double)top_forces.size();
	force_box_top /= (double)top_forces.size();
	force_bottom /= (double)top_forces.size();
	force_box_bottom /= (double)top_forces.size();

	double cross_section = sim->getCrossSection();

	double pressure_top = 0.1 * force_top/cross_section; // /10 to convert from CGS to SI
	double pressure_bottom = 0.1 * force_bottom/cross_section; // /10 to convert from CGS to SI
	double pressure_box_top = 0.1 * force_box_top/sim->box->base;
	double pressure_box_bottom = 0.1 * force_box_bottom/sim->box->base;
	double max_pressure_top = 0.1 * max_force_top/cross_section;
	
	// write to log file
	FILE* log_file = fopen(log_filename, "a");

	if(log_file)
	{	
		fprintf(log_file, "%g %g %g %g %g %g %g %g %g\n", sim->current_time, sim->getBoxFillingFactor(), sim->getFillingFactor(cross_section), cross_section / 1e4, pressure_box_top, pressure_top, pressure_box_bottom, pressure_bottom, max_pressure_top);
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
	int log_interval;
	char sample_filename[200];
	char log_filename[200];
	char material_filename[200];

	double wall_compression_modifier = 1.0;
	double wall_rolling_modifier = 1.0;
	double wall_sliding_modifier = 1.0;

	if(argc == 8 || argc == 11)
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		stop_filling_factor = atof(argv[3]);
		sprintf(sample_filename, "%s", argv[4]);
		sprintf(log_filename, "%s", argv[5]);
		log_interval = atoi(argv[6]);
		sprintf(material_filename, "%s", argv[7]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n");
		printf("-timestep -wall_speed -stop_filling_factor -sample_filename -log_filename -log_interval -material_filename\n");
		return 0;
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

	// load material file
	ErrorCode error_code;
	Simulation sim;

	error_code = sim.loadMaterial(material_filename);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", material_filename, message);
		return 0;
	}

	// setup box
	error_code = SimLib::initCompressionBox(&sim, sample_filename, false, false, false, wall_speed, stop_filling_factor, wall_compression_modifier, wall_rolling_modifier, wall_sliding_modifier);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nDuring simulation setup the following error occurred:\n%s\n", message);
		return 0;
	}

	// set wall modifiers
	sim.walls[0].rolling_modifier = wall_rolling_modifier;
	sim.walls[0].sliding_modifier = wall_sliding_modifier;
	sim.walls[1].rolling_modifier = wall_rolling_modifier;
	sim.walls[1].sliding_modifier = wall_sliding_modifier;

	// start sim
	error_code = sim.startSimulation(sim.box->height / wall_speed, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhen trying to start the simulation the following error occurred:\n%s\n", message);
		return 0;
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare log file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *log_file;
	log_file = fopen(log_filename, "w+");
	fprintf(log_file, "# number of particles: %i\n", sim.number_of_particles);
	fprintf(log_file, "# box base size (in m^2): %g\n", 1e-4 * sim.box->base);
	fprintf(log_file, "# box height (in m): %g \n", 0.01 * sim.box->height);
	fprintf(log_file, "# wall speed (in cm/s): %lf \n", wall_speed);
	fprintf(log_file, "# initial/stop filling factor: %lf %lf\n", sim.getBoxFillingFactor(), stop_filling_factor);
	//fprintf(log_file, "#\n");
	fprintf(log_file, "# wall compression/rolling/sliding modifier: %g / %g / %g\n", wall_compression_modifier, wall_rolling_modifier, wall_sliding_modifier);
	//if(sim.use_gravity)
	//	fprintf(log_file, "# Gravity enabled! - gravity strength modifier: %g\n", gravity_modifier);
	fprintf(log_file, "# time (in s)    box filling factor     cross section filling factor     cross section (in m^2)     top box pressure (in Pa)     top pressure (in Pa)     bottom box pressure (in Pa)     bottom pressure (in Pa)     max pressure (in Pa)\n");
	fclose(log_file);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int log_counter = 0;
	std::vector<double> top_forces(log_interval, 0);
	std::vector<double> top_box_forces(log_interval, 0);
	std::vector<double> bottom_forces(log_interval, 0);
	std::vector<double> bottom_box_forces(log_interval, 0);

	while(!sim.stop_simulation)
	{
		sim.update();

		top_forces[log_counter] = sim.walls[sim.box->top_wall_id].getWallForce();
		top_box_forces[log_counter] = sim.walls[sim.box->top_wall_id].getBoxVolumeForce();
		bottom_forces[log_counter] = sim.walls[sim.box->bottom_wall_id].getWallForce();
		bottom_box_forces[log_counter] = sim.walls[sim.box->bottom_wall_id].getBoxVolumeForce();
		++log_counter;

		// write data to log file
		if(log_counter == log_interval-1)
		{
			log_counter = 0;
			logPressure(&sim, log_filename, top_forces, top_box_forces, bottom_forces, bottom_box_forces);
		}
	}

    return 0;
}
