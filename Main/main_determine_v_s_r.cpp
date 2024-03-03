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
                {
                    printf("sticking\n");
                    return false;
                }
		if (velocity_up && sim.number_of_contacts + sim.number_of_wall_contacts == 0) // if particle is moving away from target and there is no contact, then it has bounced
                {
                    printf("bouncing\n");
                    return true;
                }
		sim.update();

	}

	return false;

}




int main(int argc, char **argv)
{

    // v_s_r = sticking velocity dependence of the sphere radius, is the maximum velocity at which the sphere still stick to the wall.
    // Here we determine the the v_s of spheres by shooting them at a wall with increasing speed.

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


    double particle_radius_max = 100.e-05;
    double particle_radius_min = 0.1e-05;
    double log_particle_radius_max = log(particle_radius_max); // starting radius = 10 um
    double log_particle_radius_min = log(particle_radius_min); // endingradius = 0.1 um

    double log_dr = (log_particle_radius_max - log_particle_radius_min)/15.0; // number of points


    particle_radius = exp(log_particle_radius_max);

    sim.setMaterialConstants(
                particle_radius,
                density,
                surface_energy,
                nu,
                young_mod,
                crit_rolling_displacement,
                osc_damping_factor,
                T_vis,
                rolling_modifier,
                sliding_modifier,
                twisting_modifier,
                crit_sliding_displacement_modifier,
                crit_wall_sliding_displacement_modifier
                );



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

    char V_s_filename[1024];
    sprintf(V_s_filename, "../V_s_");


    size_t len = strlen(argv[material_filename_index]);
    strncat(V_s_filename, argv[material_filename_index]+3, len-7);


#ifdef USE_DMT
    strcat(V_s_filename, "DMT.dat");
#endif

#ifdef USE_JKR
#ifdef USE_VISCOELASTIC_DAMPING
    strcat(V_s_filename, "JKRvisco.dat");
#else
    strcat(V_s_filename, "JKR.dat");
#endif
#endif

    printf("storage file: %s\n", V_s_filename);

    FILE* V_s_file = fopen(V_s_filename, "w");

    if(!V_s_file)
        return EXIT_SUCCESS;

    fprintf(V_s_file, "particle_radius  sticking_velocity\n");

    fclose(V_s_file);

	double speed_increase = 100.0;
        double speed = speed_increase;
	int sticking_velocity_found = 0;


start:

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

    {
        SimLib::collideAgglomerateWithWall(&sim, "tmp.dat", speed, 0.0, 1.000000001, false);
        printf("collide with wall\n");
    }

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

    printf("time = %e  radius = %e   speed = %e %d\n", sim.current_time/sim.end_time, particle_radius, speed, sticking_velocity_found);

    bool bounced = detect_bouncing(sim);

    // fit sticking velocity
    if (!bounced)
    {
            speed += speed_increase;
            goto start;
    }

    if (sticking_velocity_found < 17) // search for sticking velocity with higher accuracy
    {
            speed_increase /= 2.0;
            speed -= speed_increase;
            ++sticking_velocity_found;
            goto start;
    }



    V_s_file = fopen(V_s_filename, "a");

    if(!V_s_file)
        return EXIT_SUCCESS;


    fprintf(V_s_file, "%e   %e\n", particle_radius, speed);


    fclose(V_s_file);

    if(particle_radius > particle_radius_min)
    {
        particle_radius = exp(log(particle_radius) - log_dr);

        sim.setMaterialConstants(
                    particle_radius,
                    density,
                    surface_energy,
                    nu,
                    young_mod,
                    crit_rolling_displacement,
                    osc_damping_factor,
                    T_vis,
                    rolling_modifier,
                    sliding_modifier,
                    twisting_modifier,
                    crit_sliding_displacement_modifier,
                    crit_wall_sliding_displacement_modifier
                    );


                speed_increase = 100.0;
                speed = speed_increase;
		sticking_velocity_found = 0;

        goto start;
    }


    printf("\n .. simulation finished!\n");

    return EXIT_SUCCESS;
}
