#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <cstring>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double particle_radius;
extern double ENERGY_UNIT;

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	#undef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
#endif

#define TRACK_DISSIPATED_ENERGY


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

void trackPressureWavePropagation(Simulation &sim, double wall_speed, double penetration_depth, int log_interval, const char *log_filename)
{
	int log_counter = 0;
	int save_counter = 0;
	double top_pressure, bottom_pressure;
	double top_force, bottom_force;
	double kinetic_energy;

	// prepare force averaging
	double *top_forces = new double[log_interval];
	memset(top_forces, 0, sizeof(double) * log_interval);

	double *bottom_forces = new double[log_interval];
	memset(bottom_forces, 0, sizeof(double) * log_interval);

	// shock wave
	double a = 0.5 * wall_speed*wall_speed / (penetration_depth * particle_radius);
	double T = wall_speed / a;

	while( sim.current_time < T && !sim.stop_simulation)
	{
		sim.update();
		sim.walls[sim.box->top_wall_id].velocity[1] -= a * sim.timestep;
	}

	// stop top wall
	sim.walls[sim.box->top_wall_id].velocity[0] = 0;
	sim.walls[sim.box->top_wall_id].velocity[1] = 0;
	sim.walls[sim.box->top_wall_id].velocity[2] = 0;

	while(!sim.stop_simulation)
	{
		sim.update();

		top_forces[log_counter] = sim.walls[sim.box->top_wall_id].getWallForce();
		bottom_forces[log_counter] = sim.walls[sim.box->bottom_wall_id].getWallForce();
		++log_counter;
		++save_counter;

		if(log_counter == log_interval-1)
		{
			log_counter = 0;

			// calc average force/pressure
			top_force = 0;
			bottom_force = 0;

			for(int i = 0; i < log_interval; ++i)
			{
				top_force += top_forces[i];
				bottom_force += bottom_forces[i];
			}
			top_force /= (double)log_interval;
			bottom_force /= (double)log_interval;

			top_pressure = 0.1 * top_force/sim.box->base; // /10 to convert from CGS to SI
			bottom_pressure = 0.1 * bottom_force/sim.box->base; // /10 to convert from CGS to SI
			kinetic_energy = SimLib::getKineticEnergy(sim);

			FILE *log_file = fopen(log_filename, "a");
			if(log_file)
			{
				fprintf(log_file, "%g %g %g %g %g %g\n", sim.current_time, top_force/F_c, bottom_force/F_c, top_pressure, bottom_pressure, kinetic_energy*ENERGY_UNIT);
				fclose(log_file);
			}
		}
	}

	delete [] top_forces;
	delete [] bottom_forces;
}

int main(int argc, char **argv)
{
	srand(time(NULL));

	double timestep = 1e-10;
	double prep_top_slice_factor = 0;	// percentage of the upper part of the agglomerate that will be cut off
	double shock_min_speed = 200;		// (in cm/s)
	double shock_speed_interval = 50;
	int shock_speed_intervals = 1;
	double perturbation_layer_thickness = 2.0; // in particle radii
	double force_avg_interval = 20;

	const size_t max_str_size = 301;
	int log_interval = 10;
	char sample_filename[max_str_size];
	char log_filename[max_str_size];
	char material_filename[max_str_size];

	double wall_rolling_modifier = 0.0;
	double wall_sliding_modifier = 0.0;

	// sound speed is calculated based on the time it takes until this threshold is exceeded for the first time
	const int number_of_measuring_points = 4;
	double pressure_thresholds[number_of_measuring_points] = {0.1, 1, 10, 25};	// in Pa
	double timestamps[number_of_measuring_points] = {0.0, 0.0, 0.0, 0.0};
	int next_threshold = 0;

	if(argc == 10 || argc == 12)
	{
		timestep = atof(argv[1]);
		shock_min_speed = atof(argv[2]);
		shock_speed_interval = atof(argv[3]);
		shock_speed_intervals = atoi(argv[4]);
		perturbation_layer_thickness = atof(argv[5]);
		sprintf(sample_filename, "%s", argv[6]);
		sprintf(log_filename, "%s", argv[7]);
		log_interval = atoi(argv[8]);
		sprintf(material_filename, "%s", argv[9]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n");
		printf("-timestep -shock_min_speed -shock_speed_interval -shock_speed_intervals -perturbation_layer_thickness -sample_filename -log_filename -log_interval -material_filename\n");
		return 0;
	}

	if(argc == 12)
	{
		wall_rolling_modifier = atof(argv[10]);
		wall_sliding_modifier = atof(argv[11]);
	}

	printf("Performing aggregate compaction using particle simulation core v%s\n\n", CORE_VERSION);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(material_filename);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", material_filename, message);
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare log file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	error_code = SimLib::initShockwaveBox(&sim, sample_filename, false, 100, perturbation_layer_thickness * particle_radius, 0, wall_rolling_modifier, wall_sliding_modifier);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR: Could not load from sample file: %s\n", message);
		return 0;
	}

	FILE *log_file;
	log_file = fopen(log_filename, "w+");
	fprintf(log_file, "# number of particles: %i\n", sim.number_of_particles);
	fprintf(log_file, "# filling factor: %lf\n", sim.getBoxFillingFactor());
	fprintf(log_file, "# box base size (in m^2): %g\n", 1e-4 * sim.box->base);
	fprintf(log_file, "# box height (in m): %g\n", 0.01 * sim.box->height);
	fprintf(log_file, "# perturbation layer thickness (in m): %g\n\n", 0.01 * particle_radius * perturbation_layer_thickness);
	fprintf(log_file, "# bash wall speed  (m/s)      soundspeed (m/s) (threshold in Pa)");

	for(int i = 0; i < number_of_measuring_points; ++i)
		fprintf(log_file, "   %lf", pressure_thresholds[i]);

	fprintf(log_file, "\n");

	fclose(log_file);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// measure sound speed
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool stop_simulation;
	double pressure;
	double force;

	// prepare force averaging
	double *forces = new double[log_interval];
	int log_counter = 0;

	printf("Determination of sound speed has been started\n");

	for(int run = 0; run < shock_speed_intervals; ++run)
	{
		// reset counters
		log_counter = 0;
		next_threshold = 0;
		stop_simulation = false;
		memset(forces, 0, sizeof(double) * log_interval);
		memset(timestamps, 0, sizeof(double) * number_of_measuring_points);

		double shock_speed = shock_min_speed + (double)run * shock_speed_interval;
		error_code = SimLib::initShockwaveBox(&sim, sample_filename, false, shock_speed, perturbation_layer_thickness * particle_radius, 0, wall_rolling_modifier, wall_sliding_modifier);
		error_code = sim.startSimulation(sim.box->height / 500, timestep); // sound speed should be above 5 m/s

		FILE *log_file = fopen(log_filename, "a");

		if(log_file)
		{	
			fprintf(log_file, "%g", 0.01 * shock_speed);
			fclose(log_file);
		}
		else
			printf("ERROR: Failed to open %s\n", log_filename);

		while(!sim.stop_simulation && !stop_simulation)
		{
			sim.update();

			forces[log_counter] = sim.walls[sim.box->bottom_wall_id].getWallForce();
			++log_counter;

			if(log_counter == log_interval-1)
			{
				log_counter = 0;

				// calc average force/pressure
				force = 0;
				for(int i = 0; i < log_interval; ++i)
					force += forces[i];

				force /= (double)log_interval;
				pressure = 0.1 * force/sim.box->base; // /10 to convert from CGS to SI

				/*if(pressure > 0.0001)
				{
					FILE *log_file = fopen(log_filename, "a");
					fprintf(log_file, "%g %g %g\n", sim.current_time - relaxation_time, pressure, 0.01 * init_height / (sim.current_time - relaxation_time));
					fclose(log_file);
				}*/

				for(int i = 0; i < number_of_measuring_points; ++i)
				{
					if(i >= next_threshold && pressure > pressure_thresholds[i])
					{
						timestamps[i] = sim.current_time;
						++next_threshold;

						log_file = fopen(log_filename, "a");

						if(log_file)
						{	
							fprintf(log_file, " %g", 0.01 * sim.box->height / timestamps[i]);
							fclose(log_file);
						}	
					}
				}

				// stop measurement when last threshhold has been exceeded
				if(next_threshold >= number_of_measuring_points)
					stop_simulation = true;
			}
		}

		// log result
		log_file = fopen(log_filename, "a");

		if(log_file)
		{	
			fprintf(log_file, "\n");
			fclose(log_file);
		}

		printf("Run: %i finished\n", run+1);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// clean up
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	delete [] forces;

    return 0;
}
