#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "SimulationCuda.h"
#include "SimulationLib.h"

extern double k_s;
extern double mass;

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

    double timestep;
    double destination_filling_factor;
    int material_file_index = -1;

    double agg_mass = 0.0;
    double agg_size = 0.0;
    int seed_agg;
    int GPU_id = -1;

    if(argc == 7)
	{
        destination_filling_factor = atof(argv[1]);
        agg_mass = atof(argv[2]);
        agg_size = atof(argv[3]);
        material_file_index = 4;
        seed_agg = atoi(argv[5]);
        GPU_id = atoi(argv[6]);
	}
	else
	{
        printf("Wrong number of arguments!  %d || 7 Use:\n\
               -destination_filling_factor\
               -agg_mass\
               -agg_size\
               -material_filename\
               -seed_agg\
               - GPU_ID\n", argc);
		return EXIT_SUCCESS;
	}


    if(GPU_id >= 0)
        cudaSetDevice(GPU_id);

    /////////////////////////////////////////////////////////////////////////////////////
    // Files
    ////////////////////////////////////////////////////////////////////////////////////


    char* shpere1_file = new char[1024];
    char* shpere2_file = new char[1024];


    sprintf(shpere1_file, "../data/init/RCRBDsphere_");
    sprintf(shpere2_file, "../data/init/RBDsphere_");

    char buf[256];


    sprintf(buf, "_size_%d_%d_%d.dat", int(agg_mass*100), int(agg_size*1.e4), seed_agg);


    size_t len = strlen(argv[material_file_index]);
    strncat(shpere1_file, argv[material_file_index]+3, len-7);
    strncat(shpere2_file, argv[material_file_index]+3, len-7);

    strcat(shpere1_file, buf);
    strcat(shpere2_file, buf);


    // more relaxing
    ErrorCode error_code;
    SimulationCuda sim;
    error_code = sim.loadMaterial(argv[material_file_index]);


    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s    %s\n", argv[material_file_index], message, argv[material_file_index]);
        return EXIT_SUCCESS;
    }




    sim.loadFromFile(shpere1_file);

    // get timestep size
#ifdef USE_VARYING_POTENTIAL_COEFFICIENTS
    double sliding_osc_period = 2.0 * M_PI * sqrt(mass / (k_s * equilibrium_contact_radius) );
#else
    double sliding_osc_period = 2.0 * M_PI * sqrt(mass / k_s);
#endif

    timestep = sliding_osc_period / 40.0;

    printf("Timestep = %e\n", timestep);

    SimLib::centerCMS(&sim);

    double sphere_radius = 0.0;
    double sphere_gradius = 0.0;

    SimLib::getSize(sim, &sphere_gradius, &sphere_radius);

    sim.m_sphere_radius = sphere_radius + particle_radius;
    sim.m_sphere_compaction_speed = 0.0;

    error_code = sim.startSimulation(1.0, timestep);
    error_code = sim.initCuda(GPU_id);
    error_code = sim.toggleGPUMode(true);





    sim.setDampingFactor(0.99);
    for(int i = 0; i < 15000; ++i)
    {
        if(i % 250 == 0)
        {

            sim.GPUprintEnergies();
            sim.copySimDataFromGPU();
            SimLib::printEnergy(&sim, NULL);
            sim.cleanUpCuda();
            SimLib::resetContacts(&sim);
            error_code = sim.initCuda(GPU_id);
            error_code = sim.toggleGPUMode(true);
            sim.setDampingFactor(0.99);
        }
        sim.update(); // relax
    }


    sim.setDampingFactor(0.95);


    printf("\n\n\n");

    for(int i = 0; i < 1000; ++i)
    {
        if(i % 250 == 0)
        {

            sim.copySimDataFromGPU();
            SimLib::printEnergy(&sim, NULL);
            printf("relaxing %.2f\ percent done coord_number = %.2f\n", 100.0*double(i) / double(1000), 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
        }
        sim.update(); // relax
    }


    if(sim.get_use_gpu())
        sim.copySimDataFromGPU();

    sim.saveToFile(shpere2_file);

    printf(" ... finished\n");

    return EXIT_SUCCESS;
}
