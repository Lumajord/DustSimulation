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

extern double osc_damping_factor;
extern double crit_sliding_displacement_modifier;
extern double crit_wall_sliding_displacement_modifier;

extern double gravity_modifier;


/*
Returns true if particle bounced, false otherwise
*/
bool detect_bouncing(Simulation &sim)
{


	bool velocity_up = false;

	while (!sim.stop_simulation)
	{

		if (!velocity_up && sim.vel[Y_COORD(0)] >= 0.0)
		{
			velocity_up = true;
		}
		if (velocity_up && sim.vel[Y_COORD(0)] < 0.0) // if particle has bounced of object and is still moving towards it, then it sticks
			return false;

		if (velocity_up && sim.number_of_contacts + sim.number_of_wall_contacts == 0) // if particle is moving away from target and there is no contact, then it has bounced
			return true;

		sim.update();

	}

	return false;

}




int main(int argc, char **argv)
{

    // COR = coefficient of restitution, is the velocity ratio of pre collision velocity ans post collision velocity.
    // Here we determine the the COR of spheres by shooting them at a wall with increasing speed.

    // for main simulation
    double timestep;
    int material_filename_index;

    if(argc == 3)	// main mode -> measure pressure while wall is moving downwards
    {
        timestep = atof(argv[1]);
        material_filename_index = 2;
    }
    else
    {
        printf("Incorrect arguments! Use:\n\
            -timestep [s]\n\
            -material_filename\n");
        return EXIT_SUCCESS;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // prepare simulation
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("Performing aggregate compression on GPU using particle simulation core v%s\n", CORE_VERSION);

    Simulation sim;
    ErrorCode error_code = sim.loadMaterial(argv[material_filename_index]);

    sim.rand_generator.seed(0);

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_filename_index], message);
        return EXIT_SUCCESS;
    }




#ifdef TRACK_FORCES
#undef TRACK_FORCES
#endif

    char COR_filename[1024];
    sprintf(COR_filename, "../COR_");


    size_t len = strlen(argv[material_filename_index]);
    strncat(COR_filename, argv[material_filename_index]+3, len-7);


#ifdef USE_DMT
    strcat(COR_filename, "DMT.dat");
#endif

#ifdef USE_JKR
    strcat(COR_filename, "JKR.dat");
#endif


    printf("storage filename = %s\n", COR_filename);
    FILE* COR_file = fopen(COR_filename, "w");

    if(!COR_file)
        return EXIT_SUCCESS;

    fprintf(COR_file, "contact_duration in_speed out_speed   COR\n");

    fclose(COR_file);

    double speed_increase = 100.0;
    double speed = speed_increase;
    double log_ds = 0.0;
    double max_speed = 5000.0;
	int sticking_velocity_found = 0;
    bool sticking_velocity_found_bool = false;

	printf("T_vis = %e	surface energy = %e	young_mod = %e\n", T_vis, surface_energy, young_mod);

    const double init_size = particle_radius;
    error_code =  SimLib::initRandomBallisticDeposition(&sim, 1, init_size, init_size, 0.0);

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to generate the Agglomerate the following error occurred:\n%s\n", message);
        return EXIT_SUCCESS;
    }


    sim.saveToFile("tmp.dat", false);

    sim.cleanUp();

start:

    double init_timestep = 1.e-18;

    {
        SimLib::collideAgglomerateWithWall(&sim, "tmp.dat", speed, 0.0, 1.0 + (speed*init_timestep*0.999999999)/particle_radius, false);
        printf("collide with wall\n");
    }


    error_code = sim.startSimulation(2.e-7, init_timestep);
    sim.update();
    sim.update();
    sim.current_time = 0.0;

    /////////////////////////////////////////////////////////////////////////////////////////
    // setup cpu
    /////////////////////////////////////////////////////////////////////////////////////////

    printf("Setting up CPU...\n");

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
        return EXIT_SUCCESS;
    }



    // start sim
    error_code = sim.startSimulation(5.e-7, timestep);

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nFailed to start simulation - %s\n", message);
        return EXIT_SUCCESS;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // run simulation
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("\n Starting simulation ...\n");

    printf("speed = %e   %d\n", speed, sticking_velocity_found);

    bool bounced = detect_bouncing(sim);

        // fit sticking velocity (lowest velocity at which particle does not stick anymore)
    const int accuracy = 25;
    if (!bounced && sticking_velocity_found < accuracy)
	{
		speed += speed_increase;
		goto start;
	}

    if (sticking_velocity_found < accuracy) // search for sticking velocity with higher accuracy
	{
		speed_increase /= 2.0;
                speed -= speed_increase;
		++sticking_velocity_found;
		goto start;
	}

    if(!sticking_velocity_found_bool)
    {
        double log_speed_max = log(max_speed);
        double log_speed_min = log(speed);

        log_ds = (log_speed_max - log_speed_min)/5000.0; // number of points


        sticking_velocity_found_bool = true;
    }


    double out_speed = fabs(sim.vel[Y_COORD(0)]);

    COR_file = fopen(COR_filename, "a");

    if(!COR_file)
        return EXIT_SUCCESS;

    fprintf(COR_file, "%e   %e   %e  %e\n", sim.current_time, speed, out_speed, out_speed/speed);

    fclose(COR_file);

    if(speed < max_speed)
    {
        speed = exp(log(speed) + log_ds);
        goto start;
    }


    printf("\n .. simulation finished!\n");

    return EXIT_SUCCESS;
}
