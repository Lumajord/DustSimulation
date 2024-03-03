#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>

#include "Constants.h"
#include "Simulation.h"
#include "SimulationLib.h"

int main(int argc, char **argv)
{
	bool wall_collision;
	
	char path[300];
	char agglomerate_filename[300];
	char result_filename[300];
	char material_filename[300];

	if(argc > 3)
	{
		wall_collision = true;
		strcpy(material_filename, argv[1]);
		strcpy(path, argv[2]);
	}
	else
	{
		printf("Wrong parameters:\n");
		return EXIT_SUCCESS;
	}

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
	// analyze files
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> fragment_ids;							// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;						// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment
	int largest_fragment = 0;								// id of particles of the biggest fragment that has been detected so far

	for(int file = 3; file < argc; ++file)
	{
		// load data
		strcpy(agglomerate_filename, path);
		strcat(agglomerate_filename, '\\' );
		strcat(agglomerate_filename, argv[file]);

		error_code = sim.loadFromFile(agglomerate_filename);

		if(error_code = EC_OK)
		{
			fragment_ids.clear();
			size_of_fragment.clear();
			particles_of_fragment.clear();

			SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

			// determine result filename

			strcpy(result_filename, argv[file]);


		}
	}
	

#if defined _WIN64 || defined _WIN32
    system("pause");
#endif

    return EXIT_SUCCESS;
} 
