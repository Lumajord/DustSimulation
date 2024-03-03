#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double ENERGY_UNIT;
extern double gravity_modifier;
extern double particle_radius;

const int force_avg_interval = 20;

int main(int argc, char **argv)
{
	double timestep;
	double push_distance;
	double push_speed = 1.0;
	double pull_distance;
	double pull_speed;
	int log_interval;
	int snapshot_interval = 0;

	int log_filename_index;
	int sample_filename_index;
	int material_filename_index;
	int result_filename_index = 0;
	int snapshot_path_index = 0;

	if(argc == 9 || argc == 10 || argc == 12)
	{
		timestep = atof(argv[1]);
		push_distance = 1e-4 * atof(argv[2]);
		pull_distance = 1e-4 * atof(argv[3]);
		pull_speed = atof(argv[4]);
		log_interval = atoi(argv[5]);
		log_filename_index = 6;
		sample_filename_index = 7;
		material_filename_index = 8;
	}
	else
	{
		printf("ERROR: Wrong number of arguments! Use:\n-timestep (in s) -push_distance (in µm) -pull_distance (in µm) -pull_speed (in cm/s) -log_interval -log_filename -sample_filename -material_filename\n");
		return EXIT_SUCCESS;
	}

	if(argc == 10)
		result_filename_index = 9;

	if(argc == 12)
	{
		snapshot_interval = atoi(argv[10]);
		snapshot_path_index = 11;
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

	error_code = SimLib::initCompressionBox(&sim, argv[sample_filename_index], true, true, push_speed, 0.9, 1.0, 0.0, 0.0);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to load data - %s\n", message);
		return EXIT_SUCCESS;
	}

	error_code = sim.startSimulation(push_distance / push_speed, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start precompression - %s\n", message);
		return EXIT_SUCCESS;
	}

#ifdef VERBOSE
	printf(" done.\nStarting compaction to glue top/bottom wall to sample...");
	fflush(stdout);
#endif

	while(!sim.stop_simulation)
		sim.update();

#ifdef VERBOSE
	printf(" done.\nRelaxating sample...");
	fflush(stdout);
#endif

	// relaxate contacts
	SimLib::resetContacts(&sim);

#ifdef VERBOSE
	printf(" done.\nStarting pull strength test...");
	fflush(stdout);
#endif

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// start pull strength test
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef VERBOSE
	sim.broken_contacts = 0;
	sim.created_contacts = 0;
	sim.dissipated_contact_energy = 0;
	sim.dissipated_damping_energy = 0;
	sim.dissipated_rolling_energy = 0;
	sim.dissipated_sliding_energy = 0;
	sim.dissipated_twisting_energy = 0;
	sim.dissipated_wall_energy = 0;
#endif

	double initial_height = sim.box->height;
	error_code = SimLib::initPullStrengthTest(&sim, pull_speed);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to init pull strength test - %s\n", message);
		return EXIT_SUCCESS;
	}

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
		fprintf(log_file, "# box height (µm)   pull distance (µm)   top pressure (Pa)   bottom pressure (Pa)   n_c   contacts broken  contacts created  E_diss_tot  contact  rolling  sliding  twisting  wall  damping\n");
#else
		fprintf(log_file, "# box height (µm)   pull distance (µm)   top pressure (Pa)   bottom pressure (Pa)   n_c\n");
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

		double displacement = sim.box->height - initial_height;

		if(displacement > pull_distance)
			stop_simulation = true;

		++log_counter;
		++snapshot_counter;
		++force_counter;

		if(force_counter >= force_avg_interval)
			force_counter = 0;

		top_forces[force_counter] = sim.walls[sim.box->top_wall_id].total_force[1];
		bottom_forces[force_counter] = sim.walls[sim.box->bottom_wall_id].total_force[1];

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

			double top_pressure = 0.1 * avg_top_force / sim.box->base;
			double bottom_pressure = 0.1 * avg_bottom_force / sim.box->base;

			log_file = fopen(argv[log_filename_index], "a");

			if(log_file)
			{
#ifdef VERBOSE
				double coordination_number = 2.0 * (double)SimLib::getNumberOfContacts(&sim) / (double)sim.number_of_particles;
				double E_diss = sim.dissipated_contact_energy + sim.dissipated_rolling_energy + sim.dissipated_sliding_energy + sim.dissipated_twisting_energy + sim.dissipated_wall_energy + sim.dissipated_damping_energy;
				fprintf(log_file, "%g %g %g %g %g %i %i %g %g %g %g %g %g %g\n", 1e4 * sim.box->height, 1e4 * displacement, top_pressure, bottom_pressure, coordination_number, sim.broken_contacts, sim.created_contacts, E_diss, sim.dissipated_contact_energy, sim.dissipated_rolling_energy, sim.dissipated_sliding_energy, sim.dissipated_twisting_energy, sim.dissipated_wall_energy, sim.dissipated_damping_energy);
#else
				fprintf(log_file, "%g %g %g %g\n", 1e4 * sim.box->height, 1e4 * displacement, top_pressure, bottom_pressure);
#endif
				fclose(log_file);
			}
		}

		if(snapshot_interval > 0 && snapshot_counter >= snapshot_interval)
		{
			snapshot_counter = 0;

			char filename[300];
			sprintf(filename, "%s\\ts_%g.dat", argv[snapshot_path_index], 1e4 * displacement);

			sim.saveToFile(filename);
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
