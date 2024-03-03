#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <cstring>
#include <vector>

#include "SimulationCuda.h"
#include "SimulationLib.h"


extern double F_c;
extern double ENERGY_UNIT;

extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;

extern double particle_radius;
extern double density;
extern double surface_energy;
extern double nu;
extern double young_mod;
extern double crit_rolling_displacement;
extern double yield_strength;
extern double T_vis;
extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;

extern double gravity_modifier;

extern bool damping_enabled;

bool doesFileExist(const char *filename) 
{
  struct stat stFileInfo;
  int intStat;

  // attempt to get the file attributes
  intStat = stat(filename, &stFileInfo);
  
  if(intStat == 0)
	  return true;
  else 
	  return false;
}

void logPressure(SimulationCuda* sim, const char* log_filename, double force_top, double force_bottom, int log_interval)
{
	// determine pressure on top&bottom wall

    force_top /= (double)log_interval;
    force_bottom /= (double)log_interval;

	double pressure_top = 0.1 * force_top/sim->box->base; // /10 to convert from CGS to SI
	double pressure_bottom = 0.1 * force_bottom/sim->box->base; // /10 to convert from CGS to SI

	// write to log file
	FILE* log_file = fopen(log_filename, "a");

    //printf("%e  %e\n", sim->getGPUFillingFactor(), pressure_top);

	if(log_file)
	{	
        fprintf(log_file, "%.8e %.8e %.8e %.8e\n", sim->current_time, sim->getGPUFillingFactor(), pressure_top, pressure_bottom);
		fclose(log_file);
	}

    if(std::isnan(pressure_top) || std::isnan(pressure_bottom))
        sim->stop_simulation = true;

    return;
}


void logPressureNoSw(SimulationCuda* sim, const char* log_filename, double force_top, double force_bottom, int log_interval)
{

    double cross_section = 0.0;
    sim->wallInteraction.get_box_pos_from_gpu(*sim);

    sim->copyPositionsFromGPU(false);
    SimLib::getCrossSectionNoRotation(*sim, 0.1*particle_radius, cross_section);

    if(cross_section == 0.0)
        return;


    // determine pressure on top&bottom wall

    force_top /= (double)log_interval;
    force_bottom /= (double)log_interval;

    double pressure_top = 0.1 * force_top/cross_section; // /10 to convert from CGS to SI
    double pressure_bottom = 0.1 * force_bottom/cross_section; // /10 to convert from CGS to SI

    // write to log file
    FILE* log_file = fopen(log_filename, "a");

    //printf("%e  %e\n", sim->getGPUFillingFactor(), pressure_top);

    if(log_file)
    {
        fprintf(log_file, "%.8e %.8e %.8e %.8e\n", sim->current_time, (sim->number_of_particles * 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius) / (cross_section * sim->box->height), pressure_top, pressure_bottom);
        fclose(log_file);
    }

    if(std::isnan(pressure_top) || std::isnan(pressure_bottom))
        sim->stop_simulation = true;


    return;
}


int main(int argc, char **argv)
{
	// for main simulation
	double timestep;
	double wall_speed;
	double stop_filling_factor;

	double wall_compression_modifier = 1.0; 
    double wall_rolling_modifier = 0.001;
    double wall_sliding_modifier = 0.001;

	int log_interval;
	int material_filename_index;
    int do_plot = 0;
    int gpu_id = -1;
    int side_walls_ = 0;

    int seed = 0;
    int min_number_particles = 10000;

    double start_filling_factor = 0.15;
    double box_length = 0.0;


    if(argc == 10 || argc == 13)	// main mode -> measure pressure while wall is moving downwards
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		stop_filling_factor = atof(argv[3]);
        log_interval = atoi(argv[4]);
        seed = atoi(argv[5]);
        min_number_particles = atoi(argv[6]);
        material_filename_index = 7;
        start_filling_factor = atof(argv[8]);
        box_length = atof(argv[9]);
	}
	else
    {
        printf("Incorrect arguments! Use:\n\
            -timestep [s]\n\
            -wall_speed [cm/s]\n\
            -stop_filling_factor\n\
            -log_interval\n\
            -seed\n\
            -min_number_particles\n\
            -material_filename\n\
            -starting filling factor\n\
            -box length/width in [µm]\n\
            and optionally:\n\
            -do save data for plot (0 for no, other for yes)\n\
            -gpu_id (uses cpu if no gpu id is given)\n\
            -side_walls (1 for side walls, 0 for no side walls\n");
		return EXIT_SUCCESS;
	}
    if(argc == 13)
    {
        do_plot = atoi(argv[10]);
        gpu_id = atoi(argv[11]);
        side_walls_ = atoi(argv[12]);
    }



    bool side_walls = true;
    if (side_walls_ == 1)
    {
            side_walls = true;
    }
    else if (side_walls_ == 0)
    {
            side_walls = false;
    }




    double dest_fill_factor;
    if(start_filling_factor == 0.0)
        dest_fill_factor = 0.15;
    else
        dest_fill_factor = start_filling_factor;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("Performing aggregate compression on GPU using particle simulation core v%s\n", CORE_VERSION);

    SimulationCuda sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_filename_index]);

    sim.rand_generator.seed(seed);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_filename_index], message);
		return EXIT_SUCCESS;
	}

	if (gpu_id >= 0)
    {
		cudaError_t cuda_error = cudaSetDevice(gpu_id);

		if (cuda_error != cudaSuccess)
			return EC_CUDA_DEVICE_UNAVAILABLE;
	}

    /////////////////////////////////////////////////////////////////////////////////////////
    // setup agglomerate
    /////////////////////////////////////////////////////////////////////////////////////////

    double agg_size = 5.e-4;
    double agg_height = 5.e-4;


    init_agg:
    agg_size += 5.e-4;

    if(box_length > 0.0)
        agg_size = box_length*1.e-4;

    agg_height += 5.e-4;

    const double init_fill_factor = 0.175; // arbitrary chosen upper limit for the RBD filling factor

    double init_size = 20.e-4 + agg_size; // make init_size too big and cut out sample to minimize boundary effects of the agglomerate generating function
    double init_height = 40.e-4 + agg_height;



    int number_of_particles = init_size*init_size*init_height*init_fill_factor / (4.0/3.0 * 3.1415 * particle_radius*particle_radius*particle_radius);


    error_code =  SimLib::initRandomBallisticDeposition(&sim, number_of_particles, init_size, init_size, 0.0);
    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to generate the Agglomerate the following error occurred:\n%s\n", message);
        return EXIT_SUCCESS;
    }


    vec3 min_pos;
    vec3 max_pos;

    min_pos[0] = -agg_size/2;
    min_pos[1] = -agg_height/2;
    min_pos[2] = -agg_size/2;

    max_pos[0] = agg_size/2;
    max_pos[1] = agg_height/2;
    max_pos[2] = agg_size/2;


    if(!side_walls)
    {
        min_pos[0] = -agg_size/2 - 2.0*particle_radius;
        min_pos[1] = -agg_height/2;
        min_pos[2] = -agg_size/2 - 2.0*particle_radius;

        max_pos[0] = agg_size/2 + 2.0*particle_radius;
        max_pos[1] = agg_height/2;
        max_pos[2] = agg_size/2 + 2.0*particle_radius;
    }


    SimLib::sliceBox(&sim, min_pos, max_pos);


    error_code =  SimLib::reduceBoxFillingFactor(&sim, dest_fill_factor, side_walls);
    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to reduce the filling_factor the following error occurred:\n%s\n", message);
        return EXIT_SUCCESS;
    }


    if(!side_walls)
    {
        vec3 center_of_mass;
        SimLib::getCenterOfMass(&center_of_mass, sim, 0, sim.number_of_particles-1);

        SimLib::sliceCylinder(&sim, center_of_mass, 0.5*agg_size);
    }

    if(sim.number_of_particles < min_number_particles)
        goto init_agg;


    if(side_walls)
    {
        printf("Using side walls\n");
    }
    else
    {
        printf("Using no side walls\n");
    }

	double sim_time; 

	// setup sim
    error_code = SimLib::initCompressionBox(&sim, 0, side_walls, false, 0.0, stop_filling_factor, wall_compression_modifier, wall_rolling_modifier, wall_sliding_modifier, 1.0, 1.0, 1.0);


    sim_time = sim.box->height / wall_speed;


    sim.sim_info.info_storage[7] = sim.box->height * sim.getGPUFillingFactor() / stop_filling_factor;

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to setup simulation - %s\n", message);
		return EXIT_SUCCESS;
	}






	// start sim
	error_code = sim.startSimulation(sim_time, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation - %s\n", message);
		return EXIT_SUCCESS;
	}


    //////////////////////////////////////////////////////////////////////////////////
    // Bring Agglomerate to equilibrium
    //////////////////////////////////////////////////////////////////////////////////
    damping_enabled = true;
    for(int i = 0; i < 500; ++i)
    {
        sim.update();
    }

    sim.walls[sim.box->top_wall_id].velocity[1] = -wall_speed;
    damping_enabled = false;
    ///////////////////////////////////////////////////////////////////////////////////




	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare log file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    char log_filename[1024];

    if(side_walls)
        sprintf(log_filename, "../data/compression_sw_log_");
    else
        sprintf(log_filename, "../data/compression_no_sw_log_");


    size_t len = strlen(argv[material_filename_index]);
    strncat(log_filename, argv[material_filename_index]+3, len-7);


    char buf[256];

    int int_size = int(agg_height*1.e4);
    sprintf(buf, "%d_%d.dat", int_size, seed);

    strcat(log_filename, buf);

    char agg_filename[1024];

    if(side_walls)
        sprintf(agg_filename, "../data/compression_sw_agg_");
    else
        sprintf(agg_filename, "../data/compression_no_sw_agg_");


    strncat(agg_filename, argv[material_filename_index]+3, len-7);
    strcat(agg_filename, buf);

#if defined(GPU_TRACK_DISSIPATED_ENERGY) || defined(TRACK_DISSIPATED_ENERGY)

    sim.E_tot = 0.0;
    sim.dissipated_contact_energy = 0.0;
    sim.dissipated_damping_energy = 0.0;
    sim.dissipated_rolling_energy = 0.0;
    sim.dissipated_sliding_energy = 0.0;
    sim.dissipated_twisting_energy = 0.0;
    sim.dissipated_wall_energy = 0.0;
    sim.current_time = 0.0;


    char* energy_file = new char[1024];

    if(gpu_id != -1)
        sprintf(energy_file, "../data/energiesGPU_");
    else
        sprintf(energy_file, "../data/energiesCPU_");




    strncat(energy_file, argv[material_filename_index]+3, len-7);
    strcat(energy_file, buf);

    FILE *file = fopen(energy_file, "w+");
    printf("Energie file: %s\n", energy_file);
    if(file)
    {
        fprintf(file, "#time    E_tot   E_kin   E_pot\n");

        fclose(file);
    }

    sim.print_energies_interval = log_interval;
    sim.energies_filename = energy_file;
    sim.print_energies_counter = 1;

#endif




	FILE *log_file;
	log_file = fopen(log_filename, "w+");
	if (!log_file)
	{
		printf("Can't open file:	%s", log_filename);
		return EXIT_SUCCESS;
	}

    fprintf(log_file, "# number of particles: %i    number of walls: %d\n", sim.number_of_particles, sim.number_of_walls);
	fprintf(log_file, "# box base size (in m^2): %g\n", 1e-4 * sim.box->base);
	fprintf(log_file, "# box height (in m): %g \n", 0.01 * sim.box->height);
	fprintf(log_file, "# wall speed (in cm/s): %g \n", wall_speed);
	fprintf(log_file, "# initial/stop filling factor: %g %g\n", sim.getBoxFillingFactor(), stop_filling_factor);
	fprintf(log_file, "# T_vis: %g\n", T_vis);
    fprintf(log_file, "# compression/rolling/sliding modifier: %g / %g / %g\n", wall_compression_modifier, wall_rolling_modifier, wall_sliding_modifier);



    printf("# number of particles: %i   number of walls: %d\n", sim.number_of_particles, sim.number_of_walls);
    printf("# box base size (in m^2): %g\n", 1e-4 * sim.box->base);
    printf("# box height (in m): %g \n", 0.01 * sim.box->height);
    printf("# wall speed (in cm/s): %g \n", wall_speed);
    printf("# initial/stop filling factor: %g %g\n", sim.getGPUFillingFactor(), stop_filling_factor);

    if(do_plot != 0)
        printf("storing data for plotting\n");
    else
        printf("not storing data for plotting\n");

	fprintf(log_file, "#\n# time      filling factor     top pressure (in Pa)     bottom pressuren (in Pa)\n");
	fclose(log_file);


	//////////////////////////////////////////////////////////////////////////////////////
	// setup cuda again to load data from cpu
	/////////////////////////////////////////////////////////////////////////////////////
	if (gpu_id != -1)
	{
		printf("Setting up CUDA...\n");
		error_code = sim.initCuda(gpu_id);
	}
	else
	{
		printf("Setting up CPU...\n");
	}

	if (error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up CUDA the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	if (gpu_id != -1)
		error_code = sim.toggleGPUMode(true);
	else
		error_code = sim.toggleGPUMode(false);


	if (error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
		return EXIT_SUCCESS;
	}




	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("run simulation\n");
	int log_counter = 0;
    double top_force = 0.0;
    double bot_force = 0.0;

    int plot_interval = int(sim_time / (timestep * 60.0 * 40.0)) - 1;
    plot_interval *= 10;
    int timesnap = 0;
    int snap_counter = 0;


    if(do_plot != 0)
    { // print initial setup
        char buf[256];

		if (side_walls)
		{
                    sprintf(buf, "/scratch/jordan/test/sw/tmp%d.dat", snap_counter);
                    //sprintf(buf, "../data/plot/tmp%d.dat", snap_counter);
		}
		else
		{
                    sprintf(buf, "/scratch/jordan/test/no_sw/tmp%d.dat", snap_counter);
                    //sprintf(buf, "../data/plot_no_sw/tmp%d.dat", snap_counter);
		}

        if(sim.get_use_gpu())
            sim.copySimDataFromGPU();

        //sim.printGPUsimToFile(buf);
        sim.saveToFile(buf);
        //SimLib::printPositions(sim, buf, false, false);
        snap_counter++;
    }



    printf("\n Starting simulation ...\n");

    while(!sim.stop_simulation)
    {

        sim.update();


        if(do_plot != 0)
        {
            if(timesnap == plot_interval)
            {


                char buf[256];
                char buf_GPU[256];

				if (side_walls)
				{
                                    sprintf(buf, "/scratch/jordan/test/sw/tmp%d.dat", snap_counter);
                                    sprintf(buf_GPU, "/scratch/jordan/test/sw_GPU/tmp_GPU%d.dat", snap_counter);
                                    //sprintf(buf, "../data/plot/tmp%d.dat", snap_counter);
				}
				else
				{
                                    sprintf(buf, "/scratch/jordan/test/no_sw/tmp%d.dat", snap_counter);
                                    //sprintf(buf, "../data/plot_no_sw/tmp%d.dat", snap_counter);
				}


                if(sim.get_use_gpu())
                    sim.copySimDataFromGPU();

                sim.printGPUsimToFile(buf_GPU);
                sim.saveToFile(buf);
                timesnap = -1;
                snap_counter++;
                snap_counter%3;
            }

            timesnap++;
        }


        if(!sim.get_use_gpu())
        {
            top_force += sim.walls[sim.box->top_wall_id].total_force[1];
            bot_force += sim.walls[sim.box->bottom_wall_id].total_force[1];
        }


		// write data to log file
		if(log_counter == log_interval-1)
		{
            log_counter = -1;

            if(sim.get_use_gpu())
            {
                top_force += (sim.wallInteraction.getGPUTopWallForce(sim.cubPlan));
                bot_force += (sim.wallInteraction.getGPUBotWallForce(sim.cubPlan));
            }

            if(side_walls)
                logPressure(&sim, log_filename, top_force, bot_force, log_interval);
            else
                logPressureNoSw(&sim, log_filename, top_force, bot_force, log_interval);


            top_force = 0.0;
            bot_force = 0.0;

        }


        ++log_counter;
    }

    if(sim.get_use_gpu())
        sim.copySimDataFromGPU();

    sim.saveToFile(agg_filename);

    printf("\n .. simulation finished!\n");

    return EXIT_SUCCESS;
}
