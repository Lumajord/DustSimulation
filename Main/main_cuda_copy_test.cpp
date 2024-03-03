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

    sim->copyPositionsFromGPU(false);
    SimLib::getCrossSectionNoRotation(*sim, 0.2*particle_radius, cross_section);

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




    SimulationCuda sim;
    sim.loadGPUsimFromFile("./tmp4.dat", 0); // crash time: 3.0394158008255154e-03

    //sim.stop_simulation = false;

    int snap_counter = 0;

    double top_force = 0.0;
    double bot_force = 0.0;




    char log_filename[1024];

    if(sim.number_of_walls == 6)
        sprintf(log_filename, "../data/compression_sw_log_");
    else if(sim.number_of_walls == 2)
        sprintf(log_filename, "../data/compression_no_sw_log_");


    strcat(log_filename, "Alex");


    char buf[256];

    sprintf(buf, ".dat");

    strcat(log_filename, buf);

    FILE *log_file;
    log_file = fopen(log_filename, "w+");
    fclose(log_file);

    while(!sim.stop_simulation)
    {
        sim.update();


        {
            {

                printf("Timestep    %.8e    %d\n", sim.current_time, int((3.0394158008255154e-03 - sim.current_time)/sim.timestep));
                char buf[256];

                if (sim.number_of_walls == 6)
                {
                    //sprintf(buf, "/scratch/jordan/test/sw/tmp%d.dat", snap_counter);
                    sprintf(buf, "../data/plot/tmp%d.dat", snap_counter);
                }
                else if(sim.number_of_walls == 2)
                {
                    //sprintf(buf, "/scratch/jordan/test/no_sw/tmp%d.dat", snap_counter);
                    sprintf(buf, "../data/plot_no_sw/tmp%d.dat", snap_counter);
                }
                else
                    printf("Warning bad number of walls: %d\n", sim.number_of_walls);


                //if(sim.get_use_gpu())
                    //sim.copySimDataFromGPU();

                //sim.printGPUsimToFile(buf);
                snap_counter++;
            }
        }




        // write data to log file
        {

            if(sim.get_use_gpu())
            {
                top_force += (sim.wallInteraction.getGPUTopWallForce(sim.cubPlan));
                bot_force += (sim.wallInteraction.getGPUBotWallForce(sim.cubPlan));
            }

            if(sim.number_of_walls == 6)
                logPressure(&sim, log_filename, top_force, bot_force, 1);
            else
                logPressureNoSw(&sim, log_filename, top_force, bot_force, 1);


            top_force = 0.0;
            bot_force = 0.0;

        }





    }

    sim.copySimDataFromGPU();
    SimLib::printSimData(sim, "1.dat");
    sim.saveToFile("a1.dat");


    /*
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

    int seed = 0;
    double min_number_particles = 1000.0;

    int side_walls_ = 0;



    if(argc == 9 || argc == 11)	// main mode -> measure pressure while wall is moving downwards
    {
        timestep = atof(argv[1]);
        wall_speed = atof(argv[2]);
        stop_filling_factor = atof(argv[3]);
        log_interval = atoi(argv[4]);
        seed = atoi(argv[5]);
        min_number_particles = atof(argv[6]);
        material_filename_index = 7;
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
            and optionally:\n\
            -do save data for plot (0 for no, other for yes)\n\
            -gpu_id (uses cpu if no gpu id is given)\n\
            -side_walls (1 for side walls, 0 for no side walls\n");
      return EXIT_SUCCESS;
  }
  if(argc == 11)
  {
      do_plot = atoi(argv[8]);
      gpu_id = atoi(argv[9]);
      side_walls_ = atoi(argv[10]);
  }


      bool side_walls = false;
      if (side_walls_ == 1)
              side_walls = true;
      else if (side_walls_ == 0)
              side_walls = false;





    if(wall_speed <= 0.0)
    {
        printf("Warning: wall speed is zero\n");
        return 0;
    }

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




    // setup agglomerate

    double agg_size = 0.0; // convert from input [um] to simulation CGS [cm]


    init_agg:
    agg_size += 5.e-4;

    const double init_fill_factor = 0.175; // arbitrary chosen upper limit for the RBD filling factor

    const double init_size = agg_size + 10.0e-4; // make init_size bigger to guarantee that there are enough particles to fill the box
    int number_of_particles = init_size*init_size*(init_size + 20.e-4)*init_fill_factor / (4.0/3.0 * 3.1415 * particle_radius*particle_radius*particle_radius);


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
    min_pos[1] = -agg_size/2;
    min_pos[2] = -agg_size/2;

    max_pos[0] = agg_size/2;
    max_pos[1] = agg_size/2;
    max_pos[2] = agg_size/2;


    if(!side_walls)
    {
        min_pos[0] = -agg_size/2 - 2.0*particle_radius;
        min_pos[1] = -agg_size/2;
        min_pos[2] = -agg_size/2 - 2.0*particle_radius;

        max_pos[0] = agg_size/2 + 2.0*particle_radius;
        max_pos[1] = agg_size/2;
        max_pos[2] = agg_size/2 + 2.0*particle_radius;
    }


    SimLib::sliceBox(&sim, min_pos, max_pos);

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
    error_code = SimLib::initCompressionBox(&sim, 0, side_walls, false, wall_speed, stop_filling_factor, wall_compression_modifier, wall_rolling_modifier, wall_sliding_modifier, 1.0, 1.0, 1.0);


    sim_time = sim.box->height / wall_speed;

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nFailed to setup simulation - %s\n", message);
        return EXIT_SUCCESS;
    }




    /////////////////////////////////////////////////////////////////////////////////////////
    // setup cuda
    /////////////////////////////////////////////////////////////////////////////////////////

    if(gpu_id != -1)
    {
        printf("Setting up CUDA...\n");
        error_code = sim.initCuda(gpu_id);
    }
    else
    {
        printf("Setting up CPU...\n");
    }

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile setting up CUDA the following error occurred:\n%s\n", message);
        return EXIT_SUCCESS;
    }







    if(gpu_id != -1)
        error_code = sim.toggleGPUMode(true);
    else
        error_code = sim.toggleGPUMode(false);





    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
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



    sim.toggleGPUMode(true);
    error_code = sim.startSimulation(sim_time, timestep);

    printf("Num particles = %d\nWall speed = %e\n", sim.number_of_particles, sim.walls[sim.box->top_wall_id].velocity[1]);


    printf("\n\nSim1:\n");
    int i = 0;
    const int num = 1000;
    while(i<num)
    {

        sim.update();
        ++i;
    }

    printf("Temp timesteps %d\n", i);
    sim.printGPUsimToFile("tmp.dat");



    i = 0;
    while(i<num)
    {
        sim.update();
        ++i;
    }


    printf("1 timesteps %d\n", i);
    sim.copySimDataFromGPU();
    SimLib::printSimData(sim, "1.dat");
    sim.saveToFile("a1.dat");


    SimulationCuda sim2;
    sim2.loadGPUsimFromFile("tmp.dat", gpu_id);

    i = 0;
    while(i<num)
    {
        sim2.update();
        ++i;
    }

    printf("2 timesteps %d\n", i);
    sim2.copySimDataFromGPU();
    SimLib::printSimData(sim2, "2.dat");
    sim2.saveToFile("a2.dat");

    */



    return EXIT_SUCCESS;
}
