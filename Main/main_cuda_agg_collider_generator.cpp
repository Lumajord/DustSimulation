#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SimulationCuda.h"
#include "SimulationLib.h"

unsigned int hash(unsigned int x) {

    x = x + 123; // +123 used to avoid seed = 0 which returns hash(0)=0

    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

int main(int argc, char **argv)
{
    /////////////////////////////////////////////////////////////////////////////////////////
    // get params
    /////////////////////////////////////////////////////////////////////////////////////////

    double agg_size = 0.0;
    double agg_mass = 0.0;

    int material_file_index;
    unsigned int seed_agg = 0;


    if(argc == 4)
	{
        agg_mass = atof(argv[1]);
        material_file_index = 2;
        seed_agg = atoi(argv[3]);	}
	else
	{
        printf("Wrong number of arguments! (%d instead of 4) Use:\n\
 -agglomerate mass (µg)\
 -material\
 -seed_agg\n",
  argc);

        return EXIT_SUCCESS;
	}


    agg_size *= 1.e-4; // convert size from µm to cm


    /////////////////////////////////////////////////////////////////////////////////////
    // Output Files
    ////////////////////////////////////////////////////////////////////////////////////

    char* shpere1_file = new char[1024];
    char* shpere2_file = new char[1024];


    sprintf(shpere1_file, "../data/init/BAMsphere_");
    sprintf(shpere2_file, "../data/init/BAMsphereAgg_");


    char buf[256];


    sprintf(buf, "_size_%d_%d_%d.dat", int(agg_mass*100), int(agg_size*1.e4), seed_agg);


    size_t len = strlen(argv[material_file_index]);
    strncat(shpere1_file, argv[material_file_index]+3, len-7);
    strncat(shpere2_file, argv[material_file_index]+3, len-7);


    strcat(shpere1_file, buf);
    strcat(shpere2_file, buf);


	/////////////////////////////////////////////////////////////////////////////////////////
        // setup sim
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

	SimulationCuda sim;

    ErrorCode error_code = SimLib::generateBAMSphere(
                &sim,
                argv[material_file_index],
                agg_mass*1.e-6,
                seed_agg);


    sim.saveToFile(shpere1_file, true);

    /*
    SimulationCuda sim2;
    error_code = SimLib::generateBAMsphereAgglomerate(
                &sim2,
                argv[material_file_index],
                shpere2_file,
                agg_size,
                agg_mass*1.e-6,
                seed_agg
                );
    */

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up the collision the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

    delete [] shpere1_file;
    delete [] shpere2_file;


	return EXIT_SUCCESS;
}
