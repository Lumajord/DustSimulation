#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <cstring>
#include <vector>

#include "SimulationCuda.h"
#include "SimulationLib.h"

int main(int argc, char **argv)
{
        /////////////////////////////////////////////////////////////////////////////////////////
        // get params
        /////////////////////////////////////////////////////////////////////////////////////////


#ifdef SET_INPUT
        int input_file_index = -1;
        int result_file_index = -1;

        if(argc == 3)
        {
                input_file_index = 1;
                result_file_index = 2;

        }
        else
        {
                printf("Wrong number of arguments! Use:\n-input_filename -result_filename\n");
                return EXIT_SUCCESS;
        }


        /////////////////////////////////////////////////////////////////////////////////////////
        // load sim
        /////////////////////////////////////////////////////////////////////////////////////////

        printf("Loading simulation data...\n");

        SimulationCuda sim;

        ErrorCode error_code;
        printf("loading %s\n", argv[input_file_index]);

        error_code = sim.loadGPUsimFromFile(argv[input_file_index], 0);

        if(error_code != EC_OK)
        {
                char message[200];
                sim.getErrorMessage(error_code, message);
                printf("ERROR:\nDuring simulation setup the following error occurred:\n%s\n", message);
                return EXIT_SUCCESS;
        }

        printf("saving %s\n", argv[result_file_index]);

        sim.copySimDataFromGPU();
        sim.saveToFile(argv[result_file_index]);

        return EXIT_SUCCESS;

#else



    double timestep;
    double sim_time = 30.0e-6; // 30 micro seconds
    double collision_speed;
    double impact_parameter;
    double agg_size = 0.0;
    double agg_mass = 0.0;

    int material_file_index;
    int GPU_id = 0;
    unsigned int seed_agg = 0;

    int do_plot = 0;

    if(argc == 8 || argc == 9)
    {
    timestep = atof(argv[1]);
        collision_speed = atof(argv[2]);
        agg_mass = atof(argv[3]);
        impact_parameter = atof(argv[4]);
        material_file_index = 5;
        seed_agg = atoi(argv[6]);
        GPU_id = atoi(argv[7]);
    }
    else
    {
        printf("Wrong number of arguments! (%d instead of 8 || 9) Use:\n\
 -timestep\
 -collision_speed (cm/s)\
 -agglomerate mass (µg)\
 -impact_parameter\
 -material\
 -seed_agg\
 -GPU_id\
 \nor\n\
 -timestep\
 -collision_speed (cm/s)\
 -agglomerate mass (µg)\
 -impact_parameter\
 -material\
 -seed_agg\
 -GPU_id\
 -do_plot\n",
               argc);
        return EXIT_SUCCESS;
    }


    if(argc == 9)
    {
        do_plot = atoi(argv[8]);
    }

    if(GPU_id >= 0)
        cudaSetDevice(GPU_id);






    char* result_file = new char[1024];
    char* debug_file = new char[1024];
    sprintf(debug_file, "../data/GPU_2_collision_result_");
    sprintf(result_file, "../data/collision_2_result_");

    char buf[256];


    sprintf(buf, "_size_%d_%d_%d_%d_%d_%d.dat", int(timestep/1.e-12), int(agg_mass*100), int(agg_size*1.e4), int(collision_speed), int(impact_parameter*100.0), seed_agg);

    size_t len = strlen(argv[material_file_index]);

    strncat(debug_file, argv[material_file_index]+3, len-7);
    strcat(debug_file, buf);

    strncat(result_file, argv[material_file_index]+3, len-7);
    strcat(result_file, buf);

    /////////////////////////////////////////////////////////////////////////////////////////
    // load sim
    /////////////////////////////////////////////////////////////////////////////////////////

    printf("Loading simulation data...\n");

    SimulationCuda sim;

    ErrorCode error_code;
    printf("loading %s\n", debug_file);

    error_code = sim.loadGPUsimFromFile(debug_file, 0);

    if(error_code != EC_OK)
    {
            char message[200];
            sim.getErrorMessage(error_code, message);
            printf("ERROR:\nDuring simulation setup the following error occurred:\n%s   %s\n", message, debug_file);
            return EXIT_SUCCESS;
    }

    printf("saving %s\n", result_file);

    sim.copySimDataFromGPU();
    sim.saveToFile(result_file);

    delete[] result_file;
    delete[] debug_file;

    return EXIT_SUCCESS;


#endif

}
