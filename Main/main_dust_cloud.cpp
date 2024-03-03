#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double particle_radius;

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	srand(time(NULL));

	int number_of_fragements;
	double box_size;
	double min_velocity;
	double max_velocity;
	int fragments_max_size;

	int agg1_filename_index;
	int agg2_filename_index;
	double impact_velocity;
	double impact_parameter;
	double impact_delay;

	int data_path_index;
	int result_file_index;
	int material_file_index;

	if(argc == 8)
	{
		number_of_fragements = atoi(argv[1]);
		fragments_max_size = atoi(argv[2]);
		box_size = 1e-4 * atof(argv[3]);
		min_velocity = atof(argv[4]);
		max_velocity = atof(argv[5]);
		result_file_index = 6;
		material_file_index = 7;
	}
	else if(argc == 13)
	{
		agg1_filename_index = 1;
		agg2_filename_index = 2;
		impact_velocity = atof(argv[3]);
		impact_parameter = atof(argv[4]);
		impact_delay = atof(argv[5]);
		data_path_index = 6;
		number_of_fragements = atoi(argv[7]);
		box_size = 1e-4 * atof(argv[8]);
		min_velocity = atof(argv[9]);
		max_velocity = atof(argv[10]);
		result_file_index = 11;
		material_file_index = 12;
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

	Simulation sim2;
	error_code = sim2.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// add debris
	/////////////////////////////////////////////////////////////////////////////////////////
	
	vec3 v_gas = {5.0, 0.0, 0.0};

	if(argc == 8)
	{
	SimLib::initBAMAggregate(&sim, NULL, 3, 0, 0, BAM_SELECT_CLOSEST);
	for(int p = 0; p < sim.number_of_particles; ++p)
	{
			sim.vel[X_COORD(p)] = 2.0 * v_gas[0];
			sim.vel[Y_COORD(p)] = v_gas[1];
			sim.vel[Z_COORD(p)] = v_gas[2];
	}

	SimLib::initBAMAggregate(&sim2, NULL, 1, 0, 0, BAM_SELECT_CLOSEST);
	// shift agglomerate 2
	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
			sim2.pos_old[X_COORD(p)] += 8.0 * particle_radius;
			sim2.pos_old[Y_COORD(p)] += 0.25 * particle_radius;
	}

	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
			sim2.vel[X_COORD(p)] = 0.6 * v_gas[0];
			sim2.vel[Y_COORD(p)] = v_gas[1];
			sim2.vel[Z_COORD(p)] = v_gas[2];
	}

	sim.addParticlesFromSim(&sim2);
	sim2.cleanUp();

	SimLib::initBAMAggregate(&sim2, NULL, 3, 0, 0, BAM_SELECT_CLOSEST);
	// shift agglomerate 2
	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
			sim2.pos_old[X_COORD(p)] += 7e-3;
			sim2.pos_old[Y_COORD(p)] += 2e-3;
			sim2.pos_old[Z_COORD(p)] += 3e-3;
	}

	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
			sim2.vel[X_COORD(p)] = 1.45 * v_gas[0];
			sim2.vel[Y_COORD(p)] = v_gas[1];
			sim2.vel[Z_COORD(p)] = v_gas[2];
	}

	sim.addParticlesFromSim(&sim2);
	sim2.cleanUp();

	SimLib::initBAMAggregate(&sim2, NULL, 4, 0, 0, BAM_SELECT_CLOSEST);
	// shift agglomerate 2
	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
			sim2.pos_old[X_COORD(p)] += 7e-3 + 18.0 * particle_radius;
			sim2.pos_old[Y_COORD(p)] += 2e-3 - 0.15 * particle_radius; 
			sim2.pos_old[Z_COORD(p)] += 3e-3;
	}

	for(int p = 0; p < sim2.number_of_particles; ++p)
	{
			sim2.vel[X_COORD(p)] = 0.75 * v_gas[0];
			sim2.vel[Y_COORD(p)] = v_gas[1];
			sim2.vel[Z_COORD(p)] = v_gas[2];
	}

	sim.addParticlesFromSim(&sim2);
	sim2.cleanUp();


	for(int i = 0; i < number_of_fragements; ++i)
	{
		double x = (double)(rand()%101) / 100.0;
		int num_particles = 1 + (int)floor( x*x * (double)fragments_max_size + 0.5);

		sim2.cleanUp();
		SimLib::initBAMAggregate(&sim2, NULL, num_particles, 0.5, 0, BAM_SELECT_CLOSEST);

		int counter = 0;

find_pos:
		++counter;

		if(counter > 100)
		{
			printf("ERROR: Could not find valid pos!\n");
			break;
		}

		vec3 pos;
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
		double speed = min_velocity + (max_velocity - min_velocity) * ((double)rand() / ((double)RAND_MAX + 1.0));

		for(int p = 0; p < sim2.number_of_particles; ++p)
		{
			sim2.vel[X_COORD(p)] = v_gas[0] + speed * sin(theta) * cos(phi);
			sim2.vel[Y_COORD(p)] = v_gas[1] + speed * sin(theta) * sin(phi);
			sim2.vel[Z_COORD(p)] = v_gas[2] + speed * cos(theta);
		}

		sim.addParticlesFromSim(&sim2);
	}

	sim.saveToFile(argv[result_file_index]);
	}
	else if(argc == 13)
	{
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// set up collision
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		printf("Setting up collision...\n");
		char filename[300];
		char filename2[300];

#ifdef WIN32
		sprintf(filename, "%s\\agglomerates\\agglomerate_%s.dat", argv[data_path_index], argv[agg1_filename_index]);
		sprintf(filename2, "%s\\agglomerates\\agglomerate_%s.dat", argv[data_path_index], argv[agg2_filename_index]);
#else
		sprintf(filename, "%s/agglomerates/agglomerate_%s.dat", argv[data_path_index], argv[agg1_filename_index]);
		sprintf(filename2, "%s/agglomerates/agglomerate_%s.dat", argv[data_path_index], argv[agg2_filename_index]);
#endif
		double impact_distance = impact_velocity * impact_delay / particle_radius;
		SimLib::smartCollisionCMS(&sim, filename, filename2, impact_velocity, impact_parameter, impact_distance, true, true);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR:\nWhile trying to setup the collison, the following error occurred:\n%s\n", message);
			return EXIT_SUCCESS;
		}

		for(int p = 0; p < sim.number_of_particles; ++p)
		{
			sim.vel[X_COORD(p)] += v_gas[0];
			sim.vel[Y_COORD(p)] += v_gas[1];
			sim.vel[Z_COORD(p)] += v_gas[2];
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// open agglomerate table
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
#ifdef WIN32
		sprintf(filename, "%s\\agglomerates.dat", argv[data_path_index]);
#else
		sprintf(filename, "%s/agglomerates.dat", argv[data_path_index]);
#endif

		FILE *file = fopen(filename, "r");

		if(!file)
		{
			printf("ERROR: Could not open %s!\n", filename);
			return EXIT_SUCCESS;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// read agglomerate table
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		int number_of_agglomerates;
		double d_dummy;
		long int i_dummy;
		fscanf(file, "%i", &number_of_agglomerates);
		fscanf(file, "%lf", &d_dummy);
		fscanf(file, "%lf %lf %lf", &d_dummy, &d_dummy, &d_dummy);
		fscanf(file, "%li", &i_dummy);

		double abundance_tot = 0;
		double *abundances = new double[number_of_agglomerates];
		int *agglomerate_sizes = new int[number_of_agglomerates];

		for(int i = 0; i < number_of_agglomerates; ++i)
		{
			fscanf(file, "%lf %lf %lf %lf %i", &d_dummy, &d_dummy, &d_dummy, &d_dummy, &(agglomerate_sizes[i]));

			double abundance = 1.0 / (double)agglomerate_sizes[i];
			
			if(i == 0)
				abundances[0] = abundance;
			else
				abundances[i] = abundances[i-1] + abundance;

			abundance_tot += abundance;
		}

		fclose(file);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// select which aggregates to add 
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		printf("Selecting frgaments...\n");

		std::vector<int> aggregte_id_list(number_of_fragements);
		int particles = 0;

		for(int i = 0; i < number_of_fragements; ++i)
		{
			double k = (double)rand() / ((double)RAND_MAX + 1.0) * abundance_tot;

			int agglomerate_id = 0;

			if(k > abundances[number_of_agglomerates/2])
				agglomerate_id = number_of_agglomerates/2;

			for(; agglomerate_id < number_of_agglomerates; ++agglomerate_id)
			{
				if(k < abundances[agglomerate_id])
					break;
			}

			aggregte_id_list[i] = agglomerate_id;
			particles += agglomerate_sizes[agglomerate_id];
		}

		int next_particle_id = sim.number_of_particles;
		sim.addParticles(particles);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// add other aggeragtes
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		printf("Adding fragments...\n");

		for(int i = 0; i < number_of_fragements; ++i)
		{
#ifdef WIN32
			sprintf(filename, "%s\\agglomerates\\agglomerate_%i.dat", argv[data_path_index], aggregte_id_list[i]+1);
#else
			sprintf(filename, "%s/agglomerates/agglomerate_%i.dat", argv[data_path_index], aggregte_id_list[i]+1);
#endif

			sim2.loadFromFile(filename);
			int counter = 0;

find_pos2:
			++counter;

			if(counter > 100)
			{
				printf("ERROR: Could not find valid pos!\n");
				break;
			}

			vec3 pos;
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
					goto find_pos2;
			}

			// add particles
			double phi = 2.0 * M_PI * ( (double)rand() / ((double)RAND_MAX + 1.0));
			double theta = asin( 2.0 * (double)rand() / ((double)RAND_MAX + 1.0) - 1.0 ) + 0.5 * M_PI;
			double speed = min_velocity + (max_velocity - min_velocity) * ((double)rand() / ((double)RAND_MAX + 1.0));

			for(int p = 0; p < sim2.number_of_particles; ++p)
			{
				vec3 shifted_pos;
				shifted_pos[0] = sim2.pos_old[X_COORD(p)] + pos[0];
				shifted_pos[1] = sim2.pos_old[Y_COORD(p)] + pos[1];
				shifted_pos[2] = sim2.pos_old[Z_COORD(p)] + pos[2];

				sim.pos_old[X_COORD(next_particle_id+p)] = shifted_pos[0];
				sim.pos_old[Y_COORD(next_particle_id+p)] = shifted_pos[1];
				sim.pos_old[Z_COORD(next_particle_id+p)] = shifted_pos[2];

				sim.grid.addParticle(shifted_pos, next_particle_id+p);

				sim.vel[X_COORD(next_particle_id+p)] = v_gas[0] + speed * sin(theta) * cos(phi);
				sim.vel[Y_COORD(next_particle_id+p)] = v_gas[1] + speed * sin(theta) * sin(phi);
				sim.vel[Z_COORD(next_particle_id+p)] = v_gas[2] + speed * cos(theta);

				//sim.contact_list[next_particle_id] = s 
	
				// start with empty contact list of current particle
				sim.contact_list[next_particle_id+p] = NULL;

				ContactListEntry* new_cl_entry = NULL;
				ContactListEntry* cl_entry = sim2.contact_list[p];

				while(cl_entry)
				{
					Contact *new_contact = new Contact();
					*new_contact = *(cl_entry->contact);

					new_contact->id1 += next_particle_id;
					new_contact->id2 += next_particle_id;

					if(sim.contact_list[next_particle_id + p] == NULL) // particle has no other contacts yet
					{
						// create contact list entry
						sim.contact_list[next_particle_id + p] = new ContactListEntry;
						sim.contact_list[next_particle_id + p]->next = NULL;
						sim.contact_list[next_particle_id + p]->id = cl_entry->id + next_particle_id;
						sim.contact_list[next_particle_id + p]->contact = new_contact;
						new_cl_entry = sim.contact_list[next_particle_id + p];
					}
					else
					{
						new_cl_entry->next = new ContactListEntry;
						new_cl_entry->next->next = NULL;
						new_cl_entry->next->id = cl_entry->id + next_particle_id;
						new_cl_entry->next->contact = new_contact;
						new_cl_entry = new_cl_entry->next;
					}

					cl_entry = cl_entry->next;
				}
			}

			next_particle_id += sim2.number_of_particles;

			// shift agglomerate 2
			/*for(int p = 0; p < sim2.number_of_particles; ++p)
			{
				sim.pos_old[X_COORD(next_particle_id)] = sim2.pos_old[X_COORD(p)] += pos[0];
				sim2.pos_old[Y_COORD(p)] += pos[1];
				sim2.pos_old[Z_COORD(p)] += pos[2];
			}

			double phi = 2.0 * M_PI * ( (double)rand() / ((double)RAND_MAX + 1.0));
			double theta = asin( 2.0 * (double)rand() / ((double)RAND_MAX + 1.0) - 1.0 ) + 0.5 * M_PI;
			double speed = min_velocity + (max_velocity - min_velocity) * ((double)rand() / ((double)RAND_MAX + 1.0));

			for(int p = 0; p < sim2.number_of_particles; ++p)
			{
				sim2.vel[X_COORD(p)] = v_gas[0] + speed * sin(theta) * cos(phi);
				sim2.vel[Y_COORD(p)] = v_gas[1] + speed * sin(theta) * sin(phi);
				sim2.vel[Z_COORD(p)] = v_gas[2] + speed * cos(theta);
			}

			sim.addParticlesFromSim(&sim2);*/

			printf("\b\b\b\b\b\b\b\b\b\b\b\b");
			printf("%2.0lf %%", 100.0 * (double)next_particle_id / (double)(sim.number_of_particles));

		}

		delete [] abundances;
		delete [] agglomerate_sizes;
	}

	sim.saveToFile(argv[result_file_index]);

	return EXIT_SUCCESS;
}