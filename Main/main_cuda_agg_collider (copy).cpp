#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SimulationCuda.h"
#include "SimulationLib.h"



int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

	double timestep;
	double sim_time;
	double collision_speed;
	double impact_parameter;
	double impact_distance;

	int material_file_index;
	int result_file_index;
	int input_file1_index;
	int input_file2_index;
	int GPU_id = 0;
    int seed = 0;

	unsigned int snapshot_interval = 0;
    int snapshot_path_index = 0;

    if(argc == 11 || argc == 13)
	{
		timestep = atof(argv[1]);
		sim_time = atof(argv[2]);
		collision_speed = atof(argv[3]);
		impact_parameter = atof(argv[4]);
		impact_distance = atof(argv[5]);
		input_file1_index = 6;
		input_file2_index = 7;
		result_file_index = 8;
		material_file_index = 9;
        seed = atoi(argv[10]);
	}
	else
	{
        printf("Wrong number of arguments! (%d instead of 11 || 13) Use:\n-timestep\
 -sim_time\
 -collision_speed\
 -impact_parameter\
 -impact_distance (in r_p)\
 -filename1\
 -filename2\
 -result_filename\
 -material\
\nor\n\
-timestep\
 -sim_time\
 -collision_speed\
 -impact_parameter\
 -filename1\
 -filename2\
 -result_filename\
 -material\
 -screenshot_interval\
 -screenshot_path\n",
               argc);
		return EXIT_SUCCESS;
	}

    if(argc == 13)
	{
        snapshot_interval = atoi(argv[11]);
        snapshot_path_index = 12;
	}


	/////////////////////////////////////////////////////////////////////////////////////////
	// setup sim sim
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Loading simulation data...\n");

	SimulationCuda sim;
    sim.rand_generator.seed(seed);
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);


	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}


    error_code = SimLib::collideTwoAgglomerates(&sim, argv[input_file1_index], argv[input_file2_index], collision_speed, impact_parameter, impact_distance, true, true, true, true, seed);
    //error_code = sim.loadFromFile(argv[input_file1_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up the collision the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup cuda
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Setting up CUDA...\n");
	error_code = sim.initCuda(GPU_id);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up CUDA the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	error_code = sim.toggleGPUMode(true);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// start simulation
	/////////////////////////////////////////////////////////////////////////////////////////

	//sim.setPotentialVariationStopCondition(1e-5, 2.0, 500);

    error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation!\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////

	if(snapshot_interval > 0)
	{
		char filename[300];
		char buf[30];
		strcpy(filename, argv[snapshot_path_index]);
		sprintf(buf, "positions_0.dat");
		strcat(filename, buf);

		sim.copyPositionsFromGPU();
		SimLib::printPositions(sim, filename);

		printf("Total number of frames: %i\n", (int) ( (sim_time / timestep) / (double)snapshot_interval) );
	}



    char* energy_file = new char[1024];
    sprintf(energy_file, "../energiesGPU_");
    strcat(energy_file, argv[result_file_index]+3);

    FILE *file = fopen(energy_file, "w+");
    printf("Energie file: %s\n", energy_file);
    if(file)
    {
        fprintf(file, "#time    E_tot   E_kin   E_pot\n");

        fclose(file);
    }

    sim.print_energies_interval = 500;
    sim.energies_filename = energy_file;
    sim.print_energies_counter = 1;


    printf("Running simulation on GPU...\n\n");
	fflush(stdout);

	unsigned int snapshot_counter = 0;
	unsigned int snapshot_id_counter = 0;

	while(!sim.stop_simulation)
	{
		sim.update();

        ++snapshot_counter;

		if(snapshot_interval > 0 && snapshot_counter >= snapshot_interval)
		{


			snapshot_counter = 0;
			++snapshot_id_counter;



			char filename[300];
			char buf[30];

			strcpy(filename, argv[snapshot_path_index]);
			sprintf(buf, "positions_%i.dat", snapshot_id_counter);
			strcat(filename, buf);

			sim.copyPositionsFromGPU();
			SimLib::printPositions(sim, filename);
		}

	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// finalize
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("...done\n");
	fflush(stdout);

    sim.copySimDataFromGPU();
    sim.saveToFile(argv[result_file_index], true);

#ifdef GPU_TRACK_NUMBER_OF_CONTACTS
	int broken_contacts, created_contacts;
	cudaMemcpy(&broken_contacts, sim.gpu_contacts_broken, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&created_contacts, sim.gpu_contacts_created, sizeof(int), cudaMemcpyDeviceToHost);
	printf("Broken contacts: %i\nCreated contacts: %i\n", broken_contacts, created_contacts);
#endif

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to delete Simulation\n");
        return EXIT_SUCCESS;
    }

	return EXIT_SUCCESS;
}
