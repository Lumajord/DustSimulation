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
extern double mass;
extern double crit_rolling_displacement;
extern double yield_strength;
extern double T_vis;
extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;

extern double k_s;
extern double delta_0;

extern double osc_damping_factor;
extern double crit_sliding_displacement_squared;

extern double crit_sliding_displacement_modifier;
extern double crit_wall_sliding_displacement_modifier;

extern double gravity_modifier;

extern bool damping_enabled;





/*
Returns true if particle bounced, false otherwise
*/
bool detect_bouncing(Simulation &sim)
{

    if(sim.number_of_particles == 1)
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
            if (velocity_up && sim.number_of_wall_contacts == 0) // if particle is moving away from target and there is no contact, then it has bounced
                    {
                        printf("bouncing\n");
                        return true;
                    }
            sim.update();

        }

    }


    if(sim.number_of_particles == 2)
    {
        bool velocity_up = false;


        while (!sim.stop_simulation)
        {

            if (!velocity_up && sim.vel[X_COORD(1)] >= 0.0)
            {
                velocity_up = true;
            }
            if (velocity_up && sim.vel[X_COORD(1)] < 0.0) // if particle has bounced of object and is still moving towards it, then it sticks
                    {
                        printf("sticking\n");
                        return false;
                    }
            if (velocity_up && sim.number_of_contacts == 0) // if particle is moving away from target and there is no contact, then it has bounced
                    {
                        printf("bouncing\n");
                        return true;
                    }
            sim.update();

        }

    }



    return false;

}



int main(int argc, char **argv)
{


    // for main simulation
    double timestep;
    double speed;
    int material_filename_index;
    int collision_with_wall;


    if(argc == 5)	// main mode -> measure pressure while wall is moving downwards
    {
        timestep = atof(argv[1]);
        speed = atof(argv[2]);
        material_filename_index = 3;
        collision_with_wall = atoi(argv[4]);
    }
    else
    {
        printf("Incorrect arguments! Use:\n\
            -timestep [s]\n\
            -collision speed [cm/s]\n\
            -material_filename\n\
            -collision_with_wall\n");
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




    char timestep_filename[1024];
    sprintf(timestep_filename, "../dt_");


    size_t len = strlen(argv[material_filename_index]);
    strncat(timestep_filename, argv[material_filename_index]+3, len-7);



#ifdef USE_DMT
    strcat(timestep_filename, "DMT.dat");
#endif

#ifdef USE_JKR
    strcat(timestep_filename, "JKR.dat");
#endif


    printf("storage filename = %s\n", timestep_filename);
    FILE* timestep_file = fopen(timestep_filename, "w");

    if(!timestep_file)
        return EXIT_SUCCESS;

    fprintf(timestep_file, "contact_duration in_speed out_speed   COR\n");

    fclose(timestep_file);



    double timestep_max = timestep;
    double timestep_min = 1.e-12;
    double log_max = log(timestep_max);
    double log_min = log(timestep_min);

    double log_dt = (log_max - log_min)/200.0; // number of points





    double speed_increase = speed;
    int sticking_velocity_found = 0;

    double T_vis_increase = 1.e-11;
    double surface_energy_increase = 10.0;

start:

    //printf("Time = %e   T_vis = %e    speed = %e    surface_energy = %e\n", sim.current_time/sim.end_time, T_vis, speed, surface_energy);


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



    if(collision_with_wall == 1)
    {
        SimLib::collideAgglomerateWithWall(&sim, "tmp.dat", speed, 0.0, 1.0 + (speed*1.e-18*0.999999999)/particle_radius, false);
        printf("collide with wall\n");
    }
    else
    {
        SimLib::collideTwoAgglomerates(&sim, "tmp.dat", "tmp.dat", speed, 0.0, (speed*1.e-18*0.999999999)/particle_radius, false, false, false, true, 0);
        //printf("collide two particles\n");
    }

    error_code = sim.startSimulation(2.e-10, 1.e-18);
    sim.update();
    sim.update();
    sim.current_time = 0.0;

    if(sim.getNumberOfContacts() + sim.getNumberOfWallContacts() != 1)
        printf("Waring no contact was made\n");

    /////////////////////////////////////////////////////////////////////////////////////////
    // setup cpu
    /////////////////////////////////////////////////////////////////////////////////////////

    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
        return EXIT_SUCCESS;
    }



    // start sim
    error_code = sim.startSimulation(1.e-5, timestep);

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



int mode = 5;

if(mode == 0) // fit T_vis
{
    bool bounced = detect_bouncing(sim);


    if (bounced) // bouncing
    {

        T_vis += T_vis_increase;
        //surface_energy += surface_energy_increase;

        printf("time/end_time = %e  %d\n", sim.current_time/sim.end_time, sticking_velocity_found);

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

        goto start;
    }
    else if(sticking_velocity_found < 25)
    {

        ++sticking_velocity_found;
        //surface_energy_increase /= 2.0;
        //surface_energy -= surface_energy_increase;

        T_vis -= T_vis_increase;
        T_vis_increase /= 2.0;

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


        goto start;
    }
}


if(mode == 1) // measure sticking velocity
{

    bool bounced = detect_bouncing(sim);
    printf("speed = %e\n", speed);


    if (!bounced)
    {
        speed += speed_increase;
        goto start;
    }

    if (sticking_velocity_found < 25) // search for sticking velocity with higher accuracy
    {
        speed_increase /= 2.0;
        speed -= speed_increase;

        ++sticking_velocity_found;
        goto start;
    }


}



if(mode == 2)
{

    double E_start = SimLib::getTotalEnergy(&sim);

    while(!sim.stop_simulation)
    {
        sim.update();
    }

    double E_end = SimLib::getTotalEnergy(&sim);


    timestep_file = fopen(timestep_filename, "a");

    if(!timestep_file)
        return EXIT_SUCCESS;

    fprintf(timestep_file, "%.15e  %.15e   %.15e   %.15e\n", timestep, sim.current_time, E_start,  E_end);

    fclose(timestep_file);



    //printf("timestep = %e   Energy = %e cor = %e\n", timestep, (E_start - E_end)/E_start, out_speed/speed);

    if(timestep >= timestep_min)
    {
        timestep = exp(log(timestep) - log_dt);
        goto start;
    }
}
if(mode == 3)
{


#ifdef TRACK_FORCES
    char filename[1024];
    sprintf(filename, "../surface_energy_%d_timestep_%d.dat", int(surface_energy), int(timestep/1.e-14));
    FILE* track_forces = fopen(filename, "w");

    if(!track_forces)
        return EXIT_SUCCESS;


    fprintf(track_forces, "surface energy = %e    timestep = %e\n", surface_energy, timestep);
    fprintf(track_forces, "time   compression_length    normal_force   v_rel    contact_radius  viscoelastic_force  sliding_force\n");

    fclose(track_forces);

#else // !TRACK_FORCES
    printf("Error TRACK_FORCES not enabled\n");
    //return EXIT_SUCCESS;
#endif // TRACK_FORCES



    while(!sim.stop_simulation)
        sim.update();

}
if(mode == 4)
{



    char timestep_filename[1024];
    sprintf(timestep_filename, "../energy%d_", int(timestep/1.e-15));


    size_t len = strlen(argv[material_filename_index]);
    strncat(timestep_filename, argv[material_filename_index]+3, len-7);


#ifdef USE_DMT
    strcat(timestep_filename, "DMT.dat");
#endif

#ifdef USE_JKR
    strcat(timestep_filename, "JKR.dat");
#endif

    FILE *timestep_file = fopen(timestep_filename, "w");


    if(!timestep_file)
        return EXIT_SUCCESS;


    fprintf(timestep_file, "\n\n\n\n");

    fclose(timestep_file);


    SimLib::printEnergy(&sim, timestep_filename);


    while(!sim.stop_simulation)
    {
        sim.update();
        SimLib::printEnergy(&sim, timestep_filename);
    }


}
if(mode == 5) // analyze energy conservation
{

    if(sim.getNumberOfContacts() != 1)
        printf("Waring contact broken before\n");

    damping_enabled = true;

    while(!sim.stop_simulation) // damp away all energy in system
    {

        sim.update();
        double E = SimLib::getTotalEnergy(&sim);

    }

    damping_enabled = false;
    sim.stop_simulation = false;

    double U_sliding = 0.5 * k_s * crit_sliding_displacement_squared;

    double vel = 0.5*sqrt(2.0*U_sliding/mass);

    /*
    ContactListEntry *cl_entry;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        cl_entry = sim.contact_list[p];
        if(cl_entry)
        {
        printf("%e  %e  %e\n", cl_entry->contact->n1_initial[0], cl_entry->contact->n1_initial[1], cl_entry->contact->n1_initial[2]);
        }
    }
    */

    /*
    sim.vel[X_COORD(0)] = 0.2*vel;
    sim.vel[X_COORD(1)] = -0.2*vel;
    */
    sim.vel[Y_COORD(0)] = vel;
    sim.vel[Y_COORD(1)] = -vel;

    /*
    sim.vel[Z_COORD(0)] = vel;
    sim.vel[Z_COORD(1)] = -vel;
    */

#if defined(GPU_TRACK_DISSIPATED_ENERGY) || defined(TRACK_DISSIPATED_ENERGY)

    sim.dissipated_contact_energy = 0.0;
    sim.dissipated_damping_energy = 0.0;
    sim.dissipated_energy_of_particle[0] = 0.0;
    sim.dissipated_rolling_energy = 0.0;
    sim.dissipated_sliding_energy = 0.0;
    sim.dissipated_wall_energy = 0.0;

#endif

    sim.current_time = 0.0;
    double time = 1.e-6;
    error_code = sim.startSimulation(time, timestep);

    const int num_steps =  2 + int(time/timestep + 0.5);

    double* E_tot = new double[num_steps + 2];

    int i = 2;

    const bool print = false;


    char timestep_filename[1024];
    sprintf(timestep_filename, "../energy%d_", int(timestep/1.e-15));


    size_t len = strlen(argv[material_filename_index]);
    strncat(timestep_filename, argv[material_filename_index]+3, len-7);


#ifdef USE_DMT
    strcat(timestep_filename, "DMT.dat");
#endif

#ifdef USE_JKR
    strcat(timestep_filename, "JKR.dat");
#endif

    FILE *timestep_file = fopen(timestep_filename, "w");


    if(!timestep_file)
        return EXIT_SUCCESS;


    fprintf(timestep_file, "\n\n\n\n");

    fclose(timestep_file);


    E_tot[i] = SimLib::getTotalEnergy(&sim);
    if(print)
        SimLib::printEnergy(&sim, timestep_filename);


    while(!sim.stop_simulation)
    {
        sim.update();

        ++i;
        E_tot[i] = SimLib::getTotalEnergy(&sim);
        if(print)
            SimLib::printEnergy(&sim, timestep_filename);

    }

    if(sim.getNumberOfContacts() != 1)
        printf("Waring contact broke\n");

    double max_energy = -1.e+100;
    double min_energy = 1.e+100;

    E_tot[0] = 1.e+100;
    E_tot[1] = -1.e+100;


    int E_min = 0;
    int E_max = 1;

    for(int j = 3; j < num_steps-2; ++j)
    {
        if(E_tot[j] < min_energy)
        {
            min_energy = E_tot[j];
            E_min = j;
        }

        if(E_tot[j] > max_energy)
        {
            max_energy = E_tot[j];
            E_max = j;
        }
    }


    printf("%e   %e  %e\n", min_energy, E_tot[2], max_energy);

    double err = fabs((max_energy - min_energy)/max_energy*0.5);

    printf("Timestep = %e  %e  %d   %d  %d  %e\n", timestep, err, num_steps, E_max, E_min, vel);


    delete[] E_tot;

    if(timestep >= timestep_min && err >= 5.e-3)
    //if(timestep > 7.8e-11)
    {
        timestep = exp(log(timestep) - log_dt);
        goto start;
    }



}



    printf("time/end_time = %e\n", sim.current_time/sim.end_time);

    printf("speed = %e\n", speed);

    printf("T_vis = %e  surface_energy = %e\n", T_vis, surface_energy);



    printf("\n .. simulation finished!\n");

    return EXIT_SUCCESS;
}
