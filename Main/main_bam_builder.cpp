#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double particle_radius;

int main(int argc, char **argv)
{
	double migration_rate1;
	double migration_rate2;
	double slice_radius;
	int number_of_particles;
	int number_of_cakes;
	int result_file_index;
	int material_index;
	BAMSelectionMethod bam_selection_method;

	if(argc == 9)
	{
		number_of_cakes = atoi(argv[1]);
		number_of_particles = atoi(argv[2]);
		migration_rate1 = atof(argv[3]);
		migration_rate2 = atof(argv[4]);
		slice_radius = atof(argv[5]);
		bam_selection_method = (BAMSelectionMethod) atoi(argv[6]);
		result_file_index = 7;
		material_index = 8;
	}
	else if(argc == 7)
	{
		number_of_particles = atoi(argv[1]);
		migration_rate1 = atof(argv[2]);
		migration_rate2 = atof(argv[3]);
		slice_radius = atof(argv[4]);
		bam_selection_method = (BAMSelectionMethod) atoi(argv[5]);
		material_index = 6;
	}
	else
	{
		printf("ERROR: Wrong number of arguments!\nUse: -number_of_cakes -number_of_particles -migration_rate1 -migration_rate2 -slice_radius (in µm) -BAM_method (0,1,2) -result_file -material_filename\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	ErrorCode error_code;
	Simulation sim;

	error_code = sim.loadMaterial(argv[material_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// generate cakes
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(argc == 6)
	{
		sim.cleanUp();
		SimLib::initBAMAggregate(&sim, NULL, number_of_particles, migration_rate1, migration_rate2, bam_selection_method);
		SimLib::sliceSphere(&sim, slice_radius*1e-4, false);
		SimLib::filterFragments(&sim);

		double coord_number = 2.0 * (double)SimLib::getNumberOfContacts(&sim) / (double)sim.number_of_particles;
		double filling_factor = (double)sim.number_of_particles * particle_radius*particle_radius*particle_radius / (slice_radius*slice_radius*slice_radius*1e-12);

		char filename[300];
		sprintf(filename, "BAM_sphere_r%2.0lf_%.2lf_ff_%1.2lf.dat", slice_radius, coord_number, filling_factor);

		sim.saveToFile(filename);

		return EXIT_SUCCESS;
	}

	std::vector<double> coord_number(number_of_cakes);
	std::vector<double> filling_factor(number_of_cakes);

	for(int cake = 0; cake < number_of_cakes; ++cake)
	{
		sim.cleanUp();
		SimLib::initBAMAggregate(&sim, NULL, number_of_particles, migration_rate1, migration_rate2, bam_selection_method);
		SimLib::sliceSphere(&sim, slice_radius*1e-4, false);

		coord_number[cake] = 2.0 * (double)SimLib::getNumberOfContacts(&sim) / (double)sim.number_of_particles;
		filling_factor[cake] = (double)sim.number_of_particles * particle_radius*particle_radius*particle_radius / (slice_radius*slice_radius*slice_radius*1e-12);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate mean values
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double avg_coord_number = 0;
	double sigma_coord_number = 0;
	double avg_fill_factor = 0;
	double sigma_fill_factor = 0;

	for(int cake = 0; cake < number_of_cakes; ++cake)
	{
		avg_coord_number += coord_number[cake];
		avg_fill_factor += filling_factor[cake];

		sigma_coord_number += coord_number[cake]*coord_number[cake];
		sigma_fill_factor += filling_factor[cake]*filling_factor[cake];
	}

	avg_coord_number /= (double)number_of_cakes;
	avg_fill_factor /= (double)number_of_cakes;
	sigma_coord_number /= (double)number_of_cakes;
	sigma_fill_factor /= (double)number_of_cakes;

	sigma_coord_number = sqrt(sigma_coord_number - avg_coord_number*avg_coord_number);
	sigma_fill_factor = sqrt(sigma_fill_factor - avg_fill_factor*avg_fill_factor);

	FILE *file = fopen(argv[result_file_index], "w+");

	if(file)
	{
		fprintf(file, "# avg_coord_number  sigma_coord_number  avg_fill_factor  sigma_fill_factor\n");
		fprintf(file, "%g %g %g %g\n", avg_coord_number, sigma_coord_number, avg_fill_factor, sigma_fill_factor);
		fclose(file);
	}

    return EXIT_SUCCESS;
}
