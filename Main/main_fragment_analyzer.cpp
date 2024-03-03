#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include <string>

#include "Simulation.h"
#include "SimulationLib.h"

int main(int argc, char **argv)
{
	int sample_file_index = 1;
	int result_file_index = 2;
	int material_file_index = 3;
	bool append_mode = false;

	if(argc == 5)
		append_mode = true;
	else if(argc != 4)
	{
		printf("Error: Wrong arguments, use -sample_filename -results_filename -material_filename\n");
		printf("or -sample_filename -results_filename -material_filename -append\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
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

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare logging
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// open output file
	FILE *file_results;
	
	if(append_mode)
		file_results = fopen(argv[result_file_index], "a");
	else
		file_results = fopen(argv[result_file_index], "w+");

	if(!file_results)
	{
		printf("Error: Could not open %s\n", argv[result_file_index]);
		return EXIT_SUCCESS;
	}

	if(!append_mode)
		//fprintf(file_results, "# agglomerate1 (particles, gyration_radius), agglomerate2 (particles, gyration_radius), mean_gyration_radius, sigma_gyration_radius, outer_outer_radius, sigma_outer_radius\n");
		fprintf(file_results, "# number of monomers   number of fragments   largest fragment\n");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// analyze file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	error_code = sim.loadFromFile(argv[sample_file_index]);

	if(error_code != EC_OK)
	{
		char message[256];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load data from %s, the following error occurred:\n%s\n", argv[sample_file_index], message);
		return EXIT_SUCCESS;
	}

	std::vector<int> fragment_ids;							// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;						// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment
	int largest_fragment = 0;								// id of particles of the biggest fragment that has been detected so far

	SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

	// determine biggest agglomerate
	for(unsigned int agg = 1; agg < size_of_fragment.size(); ++agg)
	{
		if(size_of_fragment[agg] > size_of_fragment[largest_fragment])
			largest_fragment = (int)agg;
	}

	fprintf(file_results, "%i %i %i\n", sim.number_of_particles, (int)size_of_fragment.size(), size_of_fragment[largest_fragment]);

	fclose(file_results);

    return EXIT_SUCCESS;
} 
