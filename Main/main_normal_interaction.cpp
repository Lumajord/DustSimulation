#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double mass;
extern double density;
extern double F_c;
extern double delta_c;
extern double ENERGY_UNIT;
extern double gravity_modifier;
extern double radius_reduced;
extern double equilibrium_radius;
extern double surface_energy;
extern double young_mod_reduced;

int main(int argc, char **argv)
{
	Simulation sim;
	sim.loadMaterial("Silicate.dat");

	FILE *file = fopen("normal_force.dat", "w+");
	fprintf(file, "# U(0) = %g,   U(-delta_c) = %g\n", sim.comp_interpolation.getJKRPotentialEnergy(0)*ENERGY_UNIT, sim.comp_interpolation.getJKRPotentialEnergy(-delta_c)*ENERGY_UNIT);
	fprintf(file, "# sticking velocity: %g\n", sqrt(2.0 * (sim.comp_interpolation.getJKRPotentialEnergy(-delta_c) - sim.comp_interpolation.getJKRPotentialEnergy(0)) / mass ) );

	double k1 = 9.0 * 0.935 * 0.25 * surface_energy / density;
	double k2 = pow( 27.0/32.0* M_PI*M_PI * surface_energy*surface_energy / (young_mod_reduced*young_mod_reduced), 1.0/3.0);
	fprintf(file, "# k1: %g, k2: %g, sqrt(k1*k2): %g\n", k1, k2, sqrt(k1*k2) ); 
	fprintf(file, "# compression length    JKR contact radius / force    DMT contact radius / force\n");

	double delta = - delta_c;

	while(delta < 50.0 * delta_c)
	{
		if(delta < 0)
			fprintf(file, "%g %g %g\n", delta / delta_c, sim.comp_interpolation.getJKRContactRadius(delta) / equilibrium_radius, sim.comp_interpolation.getJKRForce(delta) / F_c);
		else
			fprintf(file, "%g %g %g %g %g\n", delta / delta_c, sim.comp_interpolation.getJKRContactRadius(delta) / equilibrium_radius, sim.comp_interpolation.getJKRForce(delta) / F_c, sim.comp_interpolation.getDMTContactRadius(delta) / equilibrium_radius, sim.comp_interpolation.getDMTForce(delta) / F_c);

		delta += 0.05 * delta_c;
	}

	fclose(file);


	file = fopen("JKR.dat", "w+");
	fprintf(file, "# compression length    contact radius [a/a_0]     force [F / F_c]     potential [U / (F_c * delta_c)]\n");
	delta = - delta_c;

	while(delta < 10.0 * delta_c)
	{
		fprintf(file, "%g %g %g %g\n", delta / delta_c, sim.comp_interpolation.getJKRContactRadius(delta) / equilibrium_radius, sim.comp_interpolation.getJKRForce(delta) / F_c, sim.comp_interpolation.getJKRPotentialEnergy(delta) * ENERGY_UNIT );
		delta += 0.02 * delta_c;
	}
	
	fclose(file);

    return 0;
}
