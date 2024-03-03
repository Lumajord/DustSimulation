#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include <string>

#include "Simulation.h"
#include "SimulationLib.h"

extern double density;
extern double particle_radius;

int main(int argc, char **argv)
{
	bool append_mode = false;

	if(argc == 8)
		append_mode = true;
	else if(argc != 7)
	{
		printf("Error: Wrong arguments, use -sample_filename -results_filename -material_filename -number_of_projectiles -projectiles_size -erosion_threshold\n");
		printf("or -sample_filename -results_filename -material_filename -number_of_projectiles -projectiles_size -erosion_threshold -append\n");
		return EXIT_SUCCESS;
	}

	int sample_file_index = 1;
	int result_file_index = 2;
	int material_file_index = 3;
	int number_of_projectiles = atoi(argv[4]);
	int projectile_size = atoi(argv[5]);
	int erosion_threshold = atoi(argv[6]);
	//int total_projectiles = atoi(argv[6]);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material & data
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[256];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	error_code = sim.loadFromFile(argv[sample_file_index]);

	if(error_code != EC_OK)
	{
		char message[256];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load aggregate %s, the following error occurred:\n%s\n", argv[sample_file_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// analyze file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int initial_particles = sim.number_of_particles - projectile_size * number_of_projectiles;

	std::vector<int> fragment_ids;							// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;						// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment

	SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

	int total_particles = 0;
	int max_fragment_size = 0;

	// fragments with more than 1000 monomers will not be considered fragments
	for(unsigned int agg = 0; agg < size_of_fragment.size(); ++agg)
	{
		if(size_of_fragment[agg] > max_fragment_size)
			max_fragment_size = size_of_fragment[agg];

		if(size_of_fragment[agg] > erosion_threshold)
			total_particles += size_of_fragment[agg];
	}

	double erosion_efficiency_max = (double)(initial_particles - max_fragment_size) / (double)(projectile_size * number_of_projectiles);
	double erosion_efficiency = (double)(initial_particles - total_particles) / (double)(projectile_size * number_of_projectiles);

	FILE *file;
	
	if(append_mode)
		file = fopen(argv[result_file_index], "a+");
	else
		file = fopen(argv[result_file_index], "w+");

	
	if(!file)
	{
		printf("Error: Could not open %s\n", argv[result_file_index]);
		return EXIT_SUCCESS;
	}

	if(!append_mode)
	{
		fprintf(file, "# aggregate: %s    erosion threshhold: %i\n", argv[sample_file_index], erosion_threshold);
		fprintf(file, "# velocity   erosion efficiency (biggest fragment / smart)\n");
	}

	//SimLib::filterFragments(&sim);
	//double coordination_number = 2.0 * sim.number_of_contacts / (double)sim.number_of_particles;
	//fprintf(file, "%i %g %g %g %g\n", total_projectiles, erosion_efficiency_max, erosion_efficiency, coordination_number);

	fprintf(file, "%g %g %g\n", sim.sim_info.info_storage[2], erosion_efficiency_max, erosion_efficiency);
	fclose(file);

    return EXIT_SUCCESS;
} 

