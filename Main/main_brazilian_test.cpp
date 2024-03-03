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

int main(int argc, char **argv)
{
	srand(time(NULL));

	double timestep;
	double wall_speed;
	double max_compression_factor;
	int log_interval;
	int sample_file_index = 4;
	int log_file_index = 6;
	int material_file_index = 7;

	if(argc == 8)
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		max_compression_factor = atof(argv[3]);
		log_interval = atoi(argv[5]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n");
		printf("-timestep -wall_speed -max_compression_dist -sample_filename -log_interval -log_file -material_filename\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Performing aggregate compression using particle simulation core v%s\n", CORE_VERSION);

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[300];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	// setup box
	error_code = SimLib::initCompressionBox(&sim, argv[sample_file_index], false, false, true, wall_speed, 0.7);

	if(error_code != EC_OK)
	{
		char message[300];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up simulation using %s the following error occurred:\n%s\n", argv[sample_file_index], message);
		return EXIT_SUCCESS;
	}

	// determine box to determine ceter filling factor
	vec3 lower_center_pos, upper_center_pos;
	sim.getEnclosingBox(&lower_center_pos, &upper_center_pos);

	double width = upper_center_pos[0] - lower_center_pos[0];
	lower_center_pos[0] += 0.4 * width;
	upper_center_pos[0] -= 0.4 * width;

	// start sim
	error_code = sim.startSimulation(sim.box->height * max_compression_factor / wall_speed, timestep);

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

	FILE* log_file = fopen(argv[log_file_index], "w+");
	fprintf(log_file, "# time (s)   sample height (m)    force top (Pa)   force bottom (Pa)   cross section (m^2)   center filling factor   center particles   center contacts\n");
	fclose(log_file);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int log_counter = 0;
	int avg_counter = 0;
	int avg_interval = 100;

	std::vector<double> top_forces(avg_interval, 0);
	std::vector<double> bottom_forces(avg_interval, 0);

	while(!sim.stop_simulation)
	{
		sim.update();

		top_forces[avg_counter] = sim.walls[sim.box->top_wall_id].getWallForce();
		bottom_forces[avg_counter] = sim.walls[sim.box->bottom_wall_id].getWallForce();

		++log_counter;
		++avg_counter;

		if(avg_counter >= avg_interval)
			avg_counter = 0;

		// write data to log file
		if(log_counter >= log_interval)
		{
			log_counter = 0;

			// determine pressure on top&bottom wall
			double force_top = 0;
			double force_bottom = 0;

			for(unsigned int i = 0; i < avg_interval; ++i)
			{
				force_top += top_forces[i];
				force_bottom += bottom_forces[i];
			}

			force_top /= (double)avg_interval;
			force_bottom /= (double)avg_interval;

			// convert from dyn to N
			force_top *= 1e-5;
			force_bottom *= 1e-5;

			// determine center filling factor
			lower_center_pos[1] = sim.walls[sim.box->bottom_wall_id].pos[1];
			upper_center_pos[1] = sim.walls[sim.box->top_wall_id].pos[1];
			double filling_factor = SimLib::getFillingFactorOfBox(sim, lower_center_pos, upper_center_pos);

			int particles, contacts;
			SimLib::getParticlesInSlice(sim, lower_center_pos, upper_center_pos, &particles, &contacts);

			// write to log file
			FILE* log_file = fopen(argv[log_file_index], "a");

			if(log_file)
			{	
				fprintf(log_file, "%g %g %g %g %g %g %i %i\n", sim.current_time, sim.getCurrentBoxHeight() * 1e-2, force_top, force_bottom, sim.getCrossSection() * 1e-4, filling_factor, particles, contacts);
				fclose(log_file);
			}
		}
	}

	char res_filename[300];
	strcpy(res_filename, argv[log_file_index]);
	strcat(res_filename, ".dat");
	sim.saveToFile(res_filename);
	
	return EXIT_SUCCESS;
}
