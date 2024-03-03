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


    printf("Input:\n");
    printf("timestep = %.2e s\n", timestep);
    printf("collision_speed = %.2e cm/s\n", collision_speed);
    printf("agg_mass = %.2e µg\n", agg_mass);
    printf("impact_parameter = %.2e\n", impact_parameter);
    printf("material_file = %s\n", argv[material_file_index]);
    printf("seed_agg = %d\n", seed_agg);
    printf("GPU_id = %d\n", GPU_id);
    printf("do_plot = %d\n", do_plot);




    agg_size *= 1.e-4; // convert size from µm to cm


    /////////////////////////////////////////////////////////////////////////////////////
    // Output Files
    ////////////////////////////////////////////////////////////////////////////////////


    char* result_file = new char[1024];
    char* debug_file = new char[1024];
    char* load_file = new char[1024];


    sprintf(load_file, "../data/GPU_2_collision_result_");
    sprintf(debug_file, "../data/GPU_collision_result_");
    sprintf(result_file, "../data/collision_result_");


    char buf[256];


    sprintf(buf, "_size_%d_%d_%d.dat", int(agg_mass*100.0), int(agg_size*1.e4), seed_agg);


    size_t len = strlen(argv[material_file_index]);


    sprintf(buf, "_size_%d_%d_%d_%d_%d_%d.dat", int(timestep/1.e-12), int(agg_mass*100), int(agg_size*1.e4), int(collision_speed), int(impact_parameter*100.0), seed_agg);


    strncat(debug_file, argv[material_file_index]+3, len-7);
    strcat(debug_file, buf);

    strncat(result_file, argv[material_file_index]+3, len-7);
    strcat(result_file, buf);

    strncat(load_file, argv[material_file_index]+3, len-7);
    strcat(load_file, buf);






    char* energy_file = new char[1024];
    sprintf(energy_file, "../data/energiesGPU_");


    strncat(energy_file, argv[material_file_index]+3, len-7);
    strcat(energy_file, buf);


	/////////////////////////////////////////////////////////////////////////////////////////
	// setup sim sim
	/////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup cuda
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

        SimulationCuda sim;

        ErrorCode error_code;
        printf("loading %s\n", load_file);

        error_code = sim.loadGPUsimFromFile(load_file, 0);
	sim.end_time += sim_time;
	sim.stop_simulation=false;
        if(error_code != EC_OK)
        {
                char message[200];
                sim.getErrorMessage(error_code, message);
                printf("ERROR:\nDuring simulation setup the following error occurred:\n%s\n", message);
                return EXIT_SUCCESS;
        }



	/////////////////////////////////////////////////////////////////////////////////////////
	// start simulation
	/////////////////////////////////////////////////////////////////////////////////////////




	/////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////


    sim.print_energies_interval = 100;
    sim.energies_filename = energy_file;
    sim.print_energies_counter = 0; // to print first timestep


    sim.min_end_time = sim.current_time + 9.e-6;
    sim.check_potential_variation_counter = 0;
    sim.check_potential_variation_interval = int(5e-6 / timestep + 0.5);


    printf("Running simulation on GPU...\n\n");
	fflush(stdout);

	unsigned int snapshot_counter = 0;
	unsigned int snapshot_id_counter = 0;
    unsigned int snapshot_interval = (unsigned int)(0.25e-6 / timestep + 0.5);



	while(!sim.stop_simulation)
	{
		sim.update();


            if(snapshot_counter == snapshot_interval)
            {

                sim.printGPUsimToFile(debug_file);
                sim.copySimDataFromGPU();
                sim.saveToFile(result_file, true);

                snapshot_counter = -1;
        }

        ++snapshot_counter;

    }

	/////////////////////////////////////////////////////////////////////////////////////////
	// finalize
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("...done\n");
	fflush(stdout);

    sim.copySimDataFromGPU();
    sim.saveToFile(result_file, true);

#ifdef GPU_TRACK_NUMBER_OF_CONTACTS
    int broken_contacts, created_contacts;
    cudaMemcpy(&broken_contacts, sim.gpu_contacts_broken, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&created_contacts, sim.gpu_contacts_created, sizeof(int), cudaMemcpyDeviceToHost);
    printf("Broken contacts: %i\nCreated contacts: %i\n", broken_contacts, created_contacts);
#endif


    delete[] result_file;
    delete[] debug_file;
    delete[] load_file;

    return EXIT_SUCCESS;
}
