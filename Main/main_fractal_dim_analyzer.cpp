#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include <string.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double particle_radius;

//#define CUMULATIVE_CROSS_SECTION

int main(int argc, char **argv)
{
	int sample_filename_index = 1;
	int result_filename_index = 2;
	int material_filename_index = 3;
	double inner_cutoff;
	double outer_cutoff;
	double radius_increment;

	if(argc == 7)
	{
		inner_cutoff = atof(argv[4]);
		outer_cutoff = atof(argv[5]);
		radius_increment = atof(argv[6]);
	}
	else
	{
		printf("Error: Wrong arguments, use -sample_filename -result_filename -material_filename -inner_cutoff (in agg_radius) -outer_cutoff -radius_increment (in particle radii)\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_filename_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_filename_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// analyze file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	error_code = sim.loadFromFile(argv[sample_filename_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load aggregate from %s, the following error occurred:\n%s\n", argv[sample_filename_index], message);
		return EXIT_SUCCESS;
	}

	radius_increment *= particle_radius;

	SimLib::centerCMS(&sim);

	double outer_radius;
	double gyration_radius;
	SimLib::getSize(sim, &gyration_radius, &outer_radius);

	// determine number of shells
	int shells = (int) ceil( outer_radius / radius_increment );
	std::vector<int> particles_in_shell(shells, 0);

	// sort particles in shells
	for(int p = 0; p < sim.number_of_particles; ++p)
	{
		double dist = sim.pos_old[X_COORD(p)] * sim.pos_old[X_COORD(p)] + sim.pos_old[Y_COORD(p)] * sim.pos_old[Y_COORD(p)] + sim.pos_old[Z_COORD(p)] * sim.pos_old[Z_COORD(p)];
		dist = sqrt(dist);

		int shell = (int)floor(dist / radius_increment);
		++particles_in_shell[shell];
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare log files
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	FILE* file = fopen(argv[result_filename_index], "w+");
	if(!file)
		return EC_FILE_NOT_FOUND;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// preparations for cross section
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef CUMULATIVE_CROSS_SECTION
	double cell_size = 0.2 * particle_radius;

	double x_min = -(outer_radius + 2.0 * particle_radius);
	double x_max = outer_radius + 2.0 *particle_radius;
	double z_min = -(outer_radius + 2.0 * particle_radius);
	double z_max = outer_radius + 2.0 * particle_radius;

	// allocate sufficient number of cells
	int x_cells = (int)( (x_max - x_min) / cell_size );
	int z_cells = (int)( (z_max - z_min) / cell_size );
	int *cells = new int[x_cells*z_cells];

	// range of cells that have to be checked when a particle is added (assuming that a particle covers more than 1 cell)
	int check_cells = (int)(particle_radius / cell_size)+1;

	double pos_x, pos_z;
	int x_id, z_id;

	vec3 axis;
	axis[0] = 1.0;
	axis[1] = 0;
	axis[2] = 0;

	// get positions of particles
	double *positions = new double[3*sim.number_of_particles];
	memcpy(positions, sim.pos_old, 3 * sim.number_of_particles * sizeof(double));

	std::vector<int> cell_count(shells);
	std::vector<double> projected_cross_sections(shells, 0);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// iterate through shells
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const unsigned int rotations = 16;

	for(unsigned int i = 0; i < rotations; ++i)
	{
		// rotation matrix
		double cos_a = cos(M_PI / (double)rotations);
		double sin_a = sin(M_PI / (double)rotations);
		double a11 = cos_a + axis[0] * axis[0] * (1.0 - cos_a);
		double a12 = axis[0] * axis[1] * (1.0 - cos_a) - axis[2] * sin_a;
		double a13 = axis[0] * axis[2] * (1.0 - cos_a) + axis[1] * sin_a;
		double a21 = axis[1] * axis[0] * (1.0 - cos_a) + axis[2] * sin_a;
		double a22 = cos_a + axis[1] * axis[1] * (1.0 - cos_a);
		double a23 = axis[1] * axis[2] * (1.0 - cos_a) - axis[0] * sin_a;
		double a31 = axis[2] * axis[0] * (1.0 - cos_a) - axis[1] * sin_a;
		double a32 = axis[2] * axis[1] * (1.0 - cos_a) + axis[0] * sin_a;
		double a33 = cos_a + axis[2] * axis[2] * (1.0 - cos_a);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// sort particles in shells
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// init empty cells
		memset(cells, 0, x_cells*z_cells * sizeof(int));
		
		for(int s = 0; s < shells; ++s)
			cell_count[s] = 0;

		vec3 old_pos;

		for(int p = 0; p < sim.number_of_particles; ++p)
		{
			// rotate coordinates
			old_pos[0] = positions[X_COORD(p)];
			old_pos[1] = positions[Y_COORD(p)];
			old_pos[2] = positions[Z_COORD(p)];

			positions[X_COORD(p)] = a11 * old_pos[0] + a12 * old_pos[1] + a13 * old_pos[2];
			positions[Y_COORD(p)] = a21 * old_pos[0] + a22 * old_pos[1] + a23 * old_pos[2]; 
			positions[Z_COORD(p)] = a31 * old_pos[0] + a32 * old_pos[1] + a33 * old_pos[2]; 

			// determine shell
			double dist = sim.pos_old[X_COORD(p)] * sim.pos_old[X_COORD(p)] + sim.pos_old[Y_COORD(p)] * sim.pos_old[Y_COORD(p)] + sim.pos_old[Z_COORD(p)] * sim.pos_old[Z_COORD(p)];
			dist = sqrt(dist);

			int shell = (int)floor(dist / radius_increment);

			// determine id of the cell where the center of the particle is located
			x_id = (int)( (positions[X_COORD(p)] - x_min) / cell_size );
			z_id = (int)( (positions[Z_COORD(p)] - z_min) / cell_size );

			// check surrounding cells
			pos_z = z_min + (double)(z_id - check_cells) * cell_size + 0.5 * cell_size;

			for(int z = z_id - check_cells; z < z_id + check_cells; ++z)
			{
				pos_x = x_min + (double)(x_id - check_cells) * cell_size + 0.5 * cell_size;

				for(int x = x_id - check_cells; x < x_id + check_cells; ++x)
				{
					// cell is considered to be covered by the particle if its center lies within the particle radius
					double dist_squared = (pos_x - positions[X_COORD(p)])*(pos_x - positions[X_COORD(p)]) + (pos_z - positions[Z_COORD(p)])*(pos_z - positions[Z_COORD(p)]);

					if(dist_squared < particle_radius*particle_radius)
					{
						// only count cells that have not been marked as occupied yet
						if(cells[x + x_cells * z] < shell+1)
						{
							++cell_count[shell];
							cells[x + x_cells * z] = shell+1;
						}
					}

					pos_x += cell_size;
				}

				pos_z += cell_size;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// calculate cross section for every shell
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		int total_cells = 0;

		for(int s = 0; s < shells; ++s)
		{
			total_cells += cell_count[s];
			projected_cross_sections[s] += (double)total_cells * cell_size*cell_size;
		}
	}

	fprintf(file, "# monomers: %i   outer_radius: %g   gyration_radius: %g\n", sim.number_of_particles, outer_radius, gyration_radius);
	fprintf(file, "# sphere radius    number of particles    particles/slice    projected cross section\n");
#else
	double cross_section;
	double coordination_number = 2.0 * (double)sim.number_of_contacts / (double)sim.number_of_particles;
	SimLib::getCrossSection(sim, 0.2 * particle_radius, 15, &cross_section, NULL);
	fprintf(file, "# monomers: %i   outer_radius: %g   gyration_radius: %g  projected_cross_section: %g  avg coordination number: %g\n", sim.number_of_particles, outer_radius, gyration_radius, cross_section, coordination_number);
	fprintf(file, "# sphere radius    number of particles    particles/slice\n");
#endif

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// print cumulative mass / cross section
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int s = 0; s < shells; ++s)
	{
#ifdef CUMULATIVE_CROSS_SECTION
		projected_cross_sections[s] /= (double)rotations;
#endif
		double radius = (double)s * radius_increment;

		if(s > 0)
		{
			particles_in_shell[s] += particles_in_shell[s-1];
			int diff = particles_in_shell[s] - particles_in_shell[s-1];
			double old_radius = (double)(s-1) * radius_increment;

			if(radius > outer_radius * inner_cutoff && radius < outer_radius * outer_cutoff)
#ifdef CUMULATIVE_CROSS_SECTION
				fprintf(file, "%g %i %g %g\n", radius, particles_in_shell[s], (double)diff * particle_radius*particle_radius*particle_radius / (radius*radius*radius - old_radius*old_radius*old_radius), projected_cross_sections[s]);
#else
				fprintf(file, "%g %i %g %g\n", radius, particles_in_shell[s], (double)diff * particle_radius*particle_radius*particle_radius / (radius*radius*radius - old_radius*old_radius*old_radius));
#endif
		}
		else
		{
			if(radius > outer_radius * inner_cutoff && radius < outer_radius * outer_cutoff)
#ifdef CUMULATIVE_CROSS_SECTION
				fprintf(file, "%g %i %g %g\n", radius, particles_in_shell[s], (double)particles_in_shell[s] * particle_radius*particle_radius*particle_radius / (radius*radius*radius), projected_cross_sections[s]);
#else
				fprintf(file, "%g %i %g %g\n", radius, particles_in_shell[s], (double)particles_in_shell[s] * particle_radius*particle_radius*particle_radius / (radius*radius*radius));
#endif
		}
	}

	fclose(file);

#ifdef CUMULATIVE_CROSS_SECTION
	delete [] cells;
	delete [] positions;
#endif

    return EXIT_SUCCESS;
} 

