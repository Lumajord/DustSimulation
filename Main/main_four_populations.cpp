#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include <string.h>

#include "Constants.h"
#include "Simulation.h"
#include "SimulationLib.h"

extern double mass;

int main(int argc, char **argv)
{
	int result_file_index = 1;
	int material_file_index = 2;
	int coordination_number_index = 3;
	
	if(argc < 5)
	{
		printf("Wrong parameters! Use:\n-result_file -material_file -coordination_number -file1 -file2 -...\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// open output files
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *result_file = fopen(argv[result_file_index], "a+");

	if(!result_file)
	{
		printf("ERROR:\nCannot open result file %s,\n", argv[result_file_index]);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// analyze files
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> fragment_ids;				// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;			// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment

	for(int i = 4; i < argc; ++i)
	{
		error_code = sim.loadFromFile(argv[i]);

		if(error_code == EC_OK)
		{
			fragment_ids.clear();
			size_of_fragment.clear();
			particles_of_fragment.clear();

			SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

			int fragments = size_of_fragment.size();

			// determine largest and second largest fragment
			int largest_fragment = 0;
			int second_largest_fragment = 0;

			for(int f = 0; f < fragments; ++f)
			{
				int fragment_size = particles_of_fragment[f].size();

				if(fragment_size > particles_of_fragment[largest_fragment].size())
				{
					second_largest_fragment = largest_fragment;
					largest_fragment = f;
				}
				else
				{
					if(fragment_size > particles_of_fragment[second_largest_fragment].size() || second_largest_fragment == largest_fragment)
						second_largest_fragment = f;
				}
			}
			
			double a1 = (double)(particles_of_fragment[largest_fragment].size()) / (double)sim.number_of_particles;
			double a2 = (double)(particles_of_fragment[second_largest_fragment].size()) / (double)sim.number_of_particles;
			
			if(largest_fragment == second_largest_fragment)
			  a2 = 0;

			fprintf(result_file, "%s %g %g %g\n", argv[coordination_number_index], sim.sim_info.info_storage[2], a1, a2 );
		}
	}

	fprintf(result_file, "\n");
	fclose(result_file);
	return EXIT_SUCCESS;
} 
