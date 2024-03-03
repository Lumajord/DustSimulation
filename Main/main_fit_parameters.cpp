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

	// for main simulation
	double timestep;
	int material_filename_index;
	double young_mod_in;
	double surface_energy_in;
	double T_vis_in;

    if (argc > 6)	// main mode -> measure pressure while wall is moving downwards
	{
		timestep = atof(argv[1]);
		material_filename_index = 2;
		young_mod_in = atof(argv[3]);
		surface_energy_in = atof(argv[4]);
		T_vis_in = atof(argv[5]);
	}
	else
	{
        printf("Incorrect arguments! Use:\n\
            -timestep [s]\n\
            -material_filename\n\
			-young_modulus\n\
			-surface_energy\n\
            -T_vis\n\
            -speeds\n");

		return EXIT_SUCCESS;
	}

    double* speeds = new double[argc - 6];

    for(int i = 6; i < argc; ++i)
        speeds[i-6] = atof(argv[i]);


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_filename_index]);


	sim.setMaterialConstants(
		particle_radius,
		density,
		surface_energy_in,
		nu,
		young_mod_in,
		crit_rolling_displacement,
		osc_damping_factor,
		T_vis_in,
		rolling_modifier,
		sliding_modifier,
		twisting_modifier,
		crit_sliding_displacement_modifier,
		crit_wall_sliding_displacement_modifier
	);


	sim.rand_generator.seed(0);


	if (error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("%e", 1.e10);
		return EXIT_SUCCESS;
	}



	const double init_size = particle_radius;
	error_code = SimLib::initRandomBallisticDeposition(&sim, 1, init_size, init_size, 0.0);

	if (error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("%e", 1.e10);
		return EXIT_SUCCESS;
	}


	sim.saveToFile("tmp.dat", false);


for(int i = 0; i < argc-6; ++i)
{

    double init_timestep = 1.e-18;

    {
        SimLib::collideAgglomerateWithWall(&sim, "tmp.dat", speeds[i], 0.0, 1.0 + (speeds[i]*init_timestep*0.999999999)/particle_radius, false);
    }


    error_code = sim.startSimulation(2.e-7, init_timestep);
    sim.update();
    sim.update();
    sim.current_time = 0.0;

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup cpu
	/////////////////////////////////////////////////////////////////////////////////////////

    // start sim
	error_code = sim.startSimulation(5.e-7, timestep);

	if (error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("%e", 1.e10);
		return EXIT_SUCCESS;
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    bool bounced = detect_bouncing(sim);


    double speed_out = fabs(sim.vel[Y_COORD(0)]);
    printf("%.16e,", speed_out);

}

	return EXIT_SUCCESS;
}
