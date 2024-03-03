#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include <string>

#include "Simulation.h"
#include "SimulationLib.h"

extern double density;
extern double particle_radius;

int main(int argc, char **argv)
{
	int sample_file_index = 1;
	int result_file_index = 2;
	int material_file_index = 3;
	bool append_mode = false;

	if(argc == 5)
		append_mode = true;
	else if(argc != 4)
	{
		printf("Error: Wrong arguments, use -sample_filename -results_filename -material_filename\n");
		printf("or -sample_filename -results_filename -material_filename -append\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[256];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare logging
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// open output file
	FILE *file_results;
	
	if(append_mode)
		file_results = fopen(argv[result_file_index], "a");
	else
		file_results = fopen(argv[result_file_index], "w+");

	if(!file_results)
	{
		printf("Error: Could not open %s\n", argv[result_file_index]);
		return EXIT_SUCCESS;
	}

	if(!append_mode)
		//fprintf(file_results, "# agglomerate1 (particles, gyration_radius), agglomerate2 (particles, gyration_radius), mean_gyration_radius, sigma_gyration_radius, outer_outer_radius, sigma_outer_radius, density\n");
		fprintf(file_results, "# number of monomeres    outer radius (m)    gyration radius (m)    avg dist to cms (m)    projected cross section (m^2)   density (g/cm^3)\n");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// analyze file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	error_code = sim.loadFromFile(argv[sample_file_index]);

	if(error_code != EC_OK)
	{
		char message[256];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load data from %s, the following error occurred:\n%s\n", argv[sample_file_index], message);
		return EXIT_SUCCESS;
	}

	SimLib::filterFragments(&sim);
	SimLib::centerCMS(&sim);

	double outer_radius;
	double gyration_radius;
	double avg_dist;
	double cross_section;
	int contact_histogram[13];

	if(sim.number_of_particles == 0)
	{
		outer_radius = 0;
		gyration_radius = 0;
		avg_dist = 0;
		cross_section = 0;
	}
	else if(sim.number_of_particles == 1)
	{
		outer_radius = particle_radius;
		gyration_radius = 0;
		avg_dist = 0;
		cross_section = M_PI * particle_radius * particle_radius;
	}
	else
	{
		/////////////////////////////////////////////////////////////////////////////////////////////
		// avg dist to cms, outer/gyration radius
		/////////////////////////////////////////////////////////////////////////////////////////////

		vec3 cms, delta_pos;
		SimLib::getCenterOfMass(&cms, sim, 0, sim.number_of_particles-1);

		avg_dist = 0;
		double result = 0;
		double dist, max_dist = 0;

		for(int p = 0; p < sim.number_of_particles; ++p)
		{
			delta_pos[0] = cms[0] - sim.pos_old[X_COORD(p)];
			delta_pos[1] = cms[1] - sim.pos_old[Y_COORD(p)];
			delta_pos[2] = cms[2] - sim.pos_old[Z_COORD(p)];
			dist = norm_squared(delta_pos);

			if(dist > max_dist)
				max_dist = dist;

			result += dist;
			avg_dist += sqrt(dist);
		}

		result /= (double)sim.number_of_particles;
		avg_dist /= (double)sim.number_of_particles;

		gyration_radius = sqrt(result);
		outer_radius = sqrt(max_dist);

		/////////////////////////////////////////////////////////////////////////////////////////////
		// cross section
		/////////////////////////////////////////////////////////////////////////////////////////////

		SimLib::getCrossSection(sim, 0.2 * particle_radius, 16, &cross_section, NULL);
	}

	SimLib::getContactHistogram(sim, contact_histogram);

	double r_c = sqrt(5.0/3.0) * gyration_radius;

	double agg_density = (double)sim.number_of_particles * density * (particle_radius*particle_radius*particle_radius) / (r_c*r_c*r_c);

	fprintf(file_results, "%i %g %g %g %g %g", sim.number_of_particles, 0.01 * outer_radius, 0.01 * gyration_radius, 0.01 * avg_dist, 0.0001 * cross_section, agg_density);
	for(int i = 0; i < 13; ++i)
		fprintf(file_results, " %i", contact_histogram[i]);
	fprintf(file_results, "\n");
	
	fclose(file_results);

    return EXIT_SUCCESS;
} 

