#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SimulationLib.h"

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	srand(time(NULL));

	double collision_speed1;
	double collision_speed2;
	double collision_speed3;
	double impact_parameter1;
	double impact_distance1;
	double impact_distance2;
	double impact_distance3;

	int material_file_index;
	int result_file_index;
	int input_file1_index;
	int input_file2_index;
	int input_file3_index;
	int input_file4_index;

	double mass1, mass2, mass3;

	int number_of_fragements;
	int fragments_max_size;

	if(argc == 13)
	{
		collision_speed1 = atof(argv[1]);
		impact_parameter1 = atof(argv[2]);
		impact_distance1 = atof(argv[3]);
		collision_speed2 = atof(argv[4]);
		impact_distance2 = atof(argv[5]);
		input_file1_index = 6;
		input_file2_index = 7;
		input_file3_index = 8;
		result_file_index = 9;
		material_file_index = 10;
		number_of_fragements = atoi(argv[11]);
		fragments_max_size = atoi(argv[12]);
	}
	else if(argc == 16)
	{
		collision_speed1 = atof(argv[1]);
		impact_parameter1 = atof(argv[2]);
		impact_distance1 = atof(argv[3]);
		collision_speed2 = atof(argv[4]);
		impact_distance2 = atof(argv[5]);
		collision_speed3 = atof(argv[6]);
		impact_distance3 = atof(argv[7]);
		input_file1_index = 8;
		input_file2_index = 9;
		input_file3_index = 10;
		input_file4_index = 11;
		result_file_index = 12;
		material_file_index = 13;
		number_of_fragements = atoi(argv[14]);
		fragments_max_size = atoi(argv[15]);
	}
	else
	{
		printf("Wrong number of arguments! Use:\n\n");
		return EXIT_SUCCESS;
	}


	/////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	// load aggregates
	/////////////////////////////////////////////////////////////////////////////////////////////

	printf("Setting up collision between first two aggregates...\n");
	error_code = SimLib::smartCollisionCMS(&sim, argv[input_file1_index], argv[input_file2_index], collision_speed1, impact_parameter1, impact_distance1, true, true);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up collision between aggregate 1 and 2, the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	mass1 = sim.number_of_particles;

	/////////////////////////////////////////////////////////////////////////////////////////////
	// add third aggregate
	/////////////////////////////////////////////////////////////////////////////////////////////

	printf("Adding third aggregate...\n");
	Simulation sim2;
	sim2.loadMaterial(argv[material_file_index]);
	error_code = sim2.loadFromFile(argv[input_file3_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to load aggregate 3:\n%s\n", message);
		return EXIT_SUCCESS;
	}
	
	mass2 = sim2.number_of_particles;

	// translate agglomerates to make sure cms = 0 
	SimLib::centerCMS(&sim2);
	SimLib::rotateSimRandomly(&sim2);
	SimLib::filterFragments(&sim2);

	// determine spatial extension
	vec3 lower, upper;
	sim2.getEnclosingBox(&lower, &upper);

	double angle = 0.15 * M_PI;

	// shift agglomerate 2
	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
		sim2.pos_old[X_COORD(p)] += particle_radius * impact_distance2 * sin(angle);
		sim2.pos_old[Z_COORD(p)] += particle_radius * impact_distance2 * cos(angle);
	}

	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
		sim2.vel[X_COORD(p)] = - collision_speed2 * sin(angle);
		sim2.vel[Z_COORD(p)] = - collision_speed2 * cos(angle);
		sim2.vel[Y_COORD(p)] = 0;
	}

	sim.addParticlesFromSim(&sim2);

	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// add fourth aggregate
	/////////////////////////////////////////////////////////////////////////////////////////////

	if(argc == 16)
	{
		double T2 = (particle_radius * impact_distance2) / collision_speed2;
		double T3 = (particle_radius * impact_distance3) / collision_speed3;

		error_code = sim2.loadFromFile(argv[input_file4_index]);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR:\nFailed to load aggregate 3:\n%s\n", message);
			return EXIT_SUCCESS;
		}

		mass3 = sim2.number_of_particles;

		// translate agglomerates to make sure cms = 0 
		SimLib::centerCMS(&sim2);
		SimLib::rotateSimRandomly(&sim2);
		SimLib::filterFragments(&sim2);

		// calculate offset
		double dist = (T3 - T2) * (mass1 + mass2) / (mass1 + mass2 + mass3) * collision_speed2;

		double angle2 = 0.25 * M_PI;

		// shift agglomerate 2
		for(int p = 0; p < sim2.number_of_particles; ++p)
		{
			sim2.pos_old[X_COORD(p)] += dist * particle_radius * sin(angle);
			sim2.pos_old[Z_COORD(p)] += dist * particle_radius * cos(angle);

			sim2.pos_old[X_COORD(p)] -= particle_radius * impact_distance3 * sin(angle2);
			sim2.pos_old[Y_COORD(p)] += particle_radius * impact_distance3 * cos(angle2);
		}

		for(int p = 0; p < sim2.number_of_particles; ++p)
		{
			sim2.vel[X_COORD(p)] = collision_speed3 * sin(angle2);
			sim2.vel[Y_COORD(p)] = -collision_speed3 * cos(angle2);
			sim2.vel[Z_COORD(p)] = 0;
		}

		sim.addParticlesFromSim(&sim2);
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// add debris
	/////////////////////////////////////////////////////////////////////////////////////////

	if(number_of_fragements > 0)
	{
		printf("Adding debris particles...\n");
		// determine size/pos of aggregates
		vec3 agg3_pos;
		agg3_pos[0] = particle_radius * impact_distance2 * sin(angle);
		agg3_pos[1] = particle_radius * impact_distance2 * cos(angle);
		agg3_pos[2] = 0;

		double agg3_radius;
		SimLib::getSize(sim2, NULL, &agg3_radius);

		vec3 v_gas = {5.0, 1.0, 1.0};

		for(int i = 0; i < number_of_fragements; ++i)
		{
			double x = (double)(rand()%101) / 100.0;
			int num_particles = 1 + (int)floor( x*x * (double)fragments_max_size + 0.5);

			SimLib::initBAMAggregate(&sim2, NULL, num_particles, 0.5, 0, BAM_SELECT_CLOSEST);

			// determine pos
			double box_size = 1.2 * particle_radius * impact_distance2;

			vec3 pos;
			bool valid_pos = false;
			int counter = 0;

find_pos:
			++counter;
			
			if(counter > 100)
				break;

			pos[0] = (2.0 * ((double)rand() / ((double)(RAND_MAX)+1.0)) - 1.0) * box_size;
			pos[1] = (2.0 * ((double)rand() / ((double)(RAND_MAX)+1.0)) - 1.0) * box_size;
			pos[2] = (2.0 * ((double)rand() / ((double)(RAND_MAX)+1.0)) - 1.0) * box_size;

			for(int p = 0; p < sim2.number_of_particles; ++p)
			{
				vec3 new_pos;
				new_pos[0] = sim2.pos_old[X_COORD(p)] + pos[0];
				new_pos[1] = sim2.pos_old[Y_COORD(p)] + pos[1];
				new_pos[2] = sim2.pos_old[Z_COORD(p)] + pos[2];
				
				if(!sim.grid.canAddParticleAt(new_pos, sim.pos_old))
					goto find_pos;
			}

			// shift agglomerate 2
			for(int p = 0; p < sim2.number_of_particles; ++p)
			{
				sim2.pos_old[X_COORD(p)] += pos[0];
				sim2.pos_old[Y_COORD(p)] += pos[1];
				sim2.pos_old[Z_COORD(p)] += pos[2];
			}

			double phi = 2.0 * M_PI * ( (double)rand() / ((double)RAND_MAX + 1.0));
			double theta = asin( 2.0 * (double)rand() / ((double)RAND_MAX + 1.0) - 1.0 ) + 0.5 * M_PI;
			double speed = 3.0 * (double)rand() / ((double)RAND_MAX + 1.0);

			for(int p = 0; p < sim2.number_of_particles; ++p)
			{
				sim2.vel[X_COORD(p)] = v_gas[0] + speed * sin(theta) * cos(phi);
				sim2.vel[Y_COORD(p)] = v_gas[1] + speed * sin(theta) * sin(phi);
				sim2.vel[Z_COORD(p)] = v_gas[2] + speed * cos(theta);
			}

			sim.addParticlesFromSim(&sim2);
		}
	}

	sim.saveToFile(argv[result_file_index]);

	return EXIT_SUCCESS;
}