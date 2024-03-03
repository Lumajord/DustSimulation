#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

#if !defined(ENABLE_FIXED_PARTICLES)
#error ENABLE_FIXED_PARTICLES must be defined!
#endif

extern double F_c;
extern double ENERGY_UNIT;
extern double gravity_modifier;
extern double particle_radius;

const int force_avg_interval = 20;

int main(int argc, char **argv)
{
	double timestep;
	double pull_distance;
	double pull_speed;
	int log_interval;
	int snapshot_interval = 0;

	int log_filename_index;
	int sample_filename_index;
	int material_filename_index;
	int result_filename_index = 0;
	int snapshot_path_index = 0;

	if(argc == 8 || argc == 9 || argc == 11)
	{
		timestep = atof(argv[1]);
		pull_distance = 1e-4 * atof(argv[2]);
		pull_speed = atof(argv[3]);
		log_interval = atoi(argv[4]);
		log_filename_index = 5;
		sample_filename_index = 6;
		material_filename_index = 7;
	}
	else
	{
		printf("ERROR: Wrong number of arguments! Use:\n-timestep (in s) -push_distance (in µm) -pull_distance (in µm) -pull_speed (in cm/s) -log_interval -log_filename -sample_filename -material_filename\n");
		return EXIT_SUCCESS;
	}

	if(argc == 9)
		result_filename_index = 8;

	if(argc == 11)
	{
		snapshot_interval = atoi(argv[9]);
		snapshot_path_index = 10;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;

	ErrorCode error_code = sim.loadMaterial(argv[material_filename_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_filename_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare setup
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef VERBOSE
	printf("Loading data...");
	fflush(stdout);
#endif

	error_code = sim.loadFromFile(argv[sample_filename_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to load simulation - %s\n", message);
		return EXIT_SUCCESS;
	}

	sim.removeWalls();
	sim.initSimData(SIM_TYPE_GENERAL);

	// determine upper/lower particles
	vec3 lower_pos, upper_pos;
	sim.getEnclosingBox(&lower_pos, &upper_pos);

	double base = (upper_pos[0] - lower_pos[0]) * (upper_pos[2] - lower_pos[2]);

	std::list<int> top_particles;
	std::list<int> bottom_particles;

	for(int p = 0; p < sim.number_of_particles; ++p)
	{
		if(sim.pos_old[Y_COORD(p)] > upper_pos[1] - 1.5 * particle_radius)
			top_particles.push_back(p);

		if(sim.pos_old[Y_COORD(p)] < lower_pos[1] + 1.5 * particle_radius)
			bottom_particles.push_back(p);
	}

	// fix bottom particles
	sim.initFixedPaticles(false);
	
	for(std::list<int>::iterator p = top_particles.begin(); p != top_particles.end(); ++p)
	{
		sim.fixateParticle(*p, true);
		sim.vel[X_COORD(*p)] = 0;
		sim.vel[Y_COORD(*p)] = pull_speed;
		sim.vel[Z_COORD(*p)] = 0;
	}
	
	for(std::list<int>::iterator p = bottom_particles.begin(); p != bottom_particles.end(); ++p)
		sim.fixateParticle(*p, true);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// start pull strength test
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef VERBOSE
	sim.contacts_broken = 0;
	sim.contacts_created = 0;
	sim.dissipated_contact_energy = 0;
	sim.dissipated_damping_energy = 0;
	sim.dissipated_rolling_energy = 0;
	sim.dissipated_sliding_energy = 0;
	sim.dissipated_twisting_energy = 0;
	sim.dissipated_wall_energy = 0;
#endif

	// determine sim time
	double sim_time = pull_distance / pull_speed;
	error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation - %s\n", message);
		return EXIT_SUCCESS;
	}

	// prepare log file
	FILE *log_file = fopen(argv[log_filename_index], "w+");

	if(!log_file)
	{
		printf("ERROR: Cannot open log file %s\n", argv[log_filename_index]);
		return EXIT_SUCCESS;
	}
	else
	{
#ifdef VERBOSE
		fprintf(log_file, "# pull distance (µm)   top pressure (Pa)   bottom pressure (Pa)   contacts broken   contacts created   E_diss_tot   E_diss_contact   E_diss_roll   E_diss_slide\n");
#else
		fprintf(log_file, "# pull distance (µm)   top pressure (Pa)   bottom pressure (Pa)\n");
#endif

#ifdef ENABLE_WALL_GLUE
		fprintf(log_file, "# Wall glue distance / factor: %g / %g\n", wall_glue_distance, wall_glue_strength);
#endif

#ifdef USE_NORMAL_INTERACTION_MODIFIER
		fprintf(log_file, "# Normal interaction modifier: %g\n", NORMAL_INTERACTION_MODIFIER);
#endif

		fclose(log_file);
	}

	bool stop_simulation = false;
	int log_counter = 0;
	int snapshot_id = 0;
	int snapshot_counter = 0;
	int force_counter = 0;
	double *top_forces = new double[force_avg_interval];
	double *bottom_forces = new double[force_avg_interval];

	for(int i = 0; i < force_avg_interval; ++i)
	{
		top_forces[i] = 0;
		bottom_forces[i] = 0;
	}

	while(!stop_simulation)
	{
		sim.update();

		double displacement = sim.current_time * pull_speed;

		if(displacement > pull_distance)
			stop_simulation = true;

		++log_counter;
		++snapshot_counter;
		++force_counter;

		if(force_counter >= force_avg_interval)
			force_counter = 0;

		for(std::list<int>::iterator p = top_particles.begin(); p != top_particles.end(); ++p)
			top_forces[force_counter] += sim.force_old[Y_COORD(*p)];

		for(std::list<int>::iterator p = bottom_particles.begin(); p != bottom_particles.end(); ++p)
			bottom_forces[force_counter] += sim.force_old[Y_COORD(*p)];

		if(log_counter >= log_interval)
		{
			log_counter = 0;

			// determine average forces
			double avg_top_force = 0;
			double avg_bottom_force = 0; 

			for(int i = 0; i < force_avg_interval; ++i)
			{
				avg_top_force += top_forces[i];
				avg_bottom_force += bottom_forces[i];
			}

			avg_top_force /= (double)force_avg_interval;
			avg_bottom_force /= (double)force_avg_interval;

			double top_pressure = 0.1 * avg_top_force / base;
			double bottom_pressure = 0.1 * avg_bottom_force / base;

			log_file = fopen(argv[log_filename_index], "a");

			if(log_file)
			{
#ifdef VERBOSE
				double E_diss = sim.dissipated_contact_energy + sim.dissipated_rolling_energy + sim.dissipated_sliding_energy + sim.dissipated_twisting_energy + sim.dissipated_wall_energy + sim.dissipated_damping_energy;
				fprintf(log_file, "%g %g %g %i %i %g %g %g %g\n", 1e4 * displacement, top_pressure, bottom_pressure, sim.contacts_broken, sim.contacts_created, E_diss, sim.dissipated_contact_energy, sim.dissipated_rolling_energy, sim.dissipated_sliding_energy);
#else
				fprintf(log_file, "%g %g %g\n", 1e4 * displacement, top_pressure, bottom_pressure);
#endif
				fclose(log_file);
			}
		}

		if(snapshot_interval > 0 && snapshot_counter >= snapshot_interval)
		{
			snapshot_counter = 0;
			++snapshot_id;

			char filename[300];
			sprintf(filename, "%spositions_%i.dat", argv[snapshot_path_index], snapshot_id);

			SimLib::printPositions(sim, filename, false, true);
			//sim.saveToFile(filename);
		}
	}

	delete [] top_forces;
	delete [] bottom_forces;

#ifdef VERBOSE
	printf(" done.\n");
	fflush(stdout);
#endif

	if(argc == 10 || argc == 12)
		sim.saveToFile(argv[result_filename_index]);

    return EXIT_SUCCESS;
}
