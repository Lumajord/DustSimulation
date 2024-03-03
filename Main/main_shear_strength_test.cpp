#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double ENERGY_UNIT;
extern double gravity_modifier;
extern double particle_radius;

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

	error_code = SimLib::initCompressionBox(&sim, argv[sample_filename_index], true, true, push_speed, 0.9);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to load data - %s\n", message);
		return EXIT_SUCCESS;
	}

	sim.updateBox();
	double area = sim.box->base;

	error_code = SimLib::initCompressionBox(&sim, argv[sample_filename_index], false, true, push_speed, 0.9);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to determine base area of iniatial box - %s\n", message);
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

	while(!sim.stop_simulation)
		sim.update();

	// relaxate contacts
	SimLib::resetContacts(&sim);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// start shear strength test
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	error_code = SimLib::initShearStrengthTest(&sim, pull_speed);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to init pull strength test - %s\n", message);
		return EXIT_SUCCESS;
	}

	// determine sim time
	double sim_time = pull_distance / pull_speed;

	error_code = sim.startSimulation(sim_time, timestep, 0, snapshot_interval, NULL, argv[snapshot_path_index], false);

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

		fprintf(log_file, "# pull distance (µm)   top pressure (Pa)   bottom pressure (Pa)   top load (Pa)   bottom load (Pa)\n");

#ifdef ENABLE_WALL_GLUE
		fprintf(log_file, "# Wall glue distance / factor: %g / %g\n", wall_glue_distance, wall_glue_strength);
#endif
		fclose(log_file);
	}

	double time = 0;
	bool stop_simulation = false;
	int log_counter = 0;
	int avg_counter = 0;
	const int avg_interval = 20;
	double *top_forces = new double[avg_interval];
	double *bottom_forces = new double[avg_interval];
	double *top_normal_forces = new double[avg_interval];
	double *bottom_normal_forces = new double[avg_interval];
	memset(top_forces, 0, avg_interval * sizeof(double));
	memset(bottom_forces, 0, avg_interval * sizeof(double));
	memset(top_normal_forces, 0, avg_interval * sizeof(double));
	memset(bottom_normal_forces, 0, avg_interval * sizeof(double));

	double avg_top_force = 0;
	double avg_bottom_force = 0;
	double avg_top_normal_force = 0;
	double avg_bottom_normal_force = 0;

	while(!stop_simulation)
	{
		sim.update();

		time += timestep;
		double displacement = time * pull_speed;

		if(displacement > pull_distance)
			stop_simulation = true;

		top_forces[avg_counter] = sim.walls[sim.box->top_wall_id].total_force[0];
		bottom_forces[avg_counter] = sim.walls[sim.box->bottom_wall_id].total_force[0];
		top_normal_forces[avg_counter] = sim.walls[sim.box->top_wall_id].total_force[1];
		bottom_normal_forces[avg_counter] = sim.walls[sim.box->bottom_wall_id].total_force[1];

		++avg_counter;
		++log_counter;
		//++snapshot_counter;

		if(avg_counter >= avg_interval)
		{
			avg_counter = 0;
			avg_top_force = 0;
			avg_bottom_force = 0;
			avg_top_normal_force = 0;
			avg_bottom_normal_force = 0;

			for(int i = 0; i < avg_interval; ++i)
			{
				avg_top_force += top_forces[i];
				avg_bottom_force += bottom_forces[i];
				avg_top_normal_force += top_normal_forces[i];
				avg_bottom_normal_force += bottom_normal_forces[i];
			}

			avg_top_force /= (double)avg_interval;
			avg_bottom_force /= (double)avg_interval;
			avg_top_normal_force /= (double)avg_interval;
			avg_bottom_normal_force /= (double)avg_interval;
		}

		if(log_counter >= log_interval)
		{
			log_counter = 0;
			double top_pressure = 0.1 * avg_top_force / area;
			double bottom_pressure = 0.1 * avg_bottom_force / area;
			double top_load = 0.1 * avg_top_normal_force / area;
			double bottom_load = 0.1 * avg_bottom_normal_force / area;

			log_file = fopen(argv[log_filename_index], "a");

			if(log_file)
			{
				fprintf(log_file, "%g %g %g %g %g\n", 1e4 * displacement, top_pressure, bottom_pressure, top_load, bottom_load);
				fclose(log_file);
			}
		}

		/*if(snapshot_interval > 0 && snapshot_counter >= snapshot_interval)
		{
			snapshot_counter = 0;
			++snapshot_id;

			char filename[300];
			sprintf(filename, "%s\\ts_%i.dat", argv[snapshot_path_index], snapshot_id);

			SimLib::SimLib::printPositions(this, filename);
		}*/
	}

	if(argc == 10 || argc == 12)
		sim.saveToFile(argv[result_filename_index]);

	delete [] top_forces;
	delete [] bottom_forces;

    return EXIT_SUCCESS;
}

