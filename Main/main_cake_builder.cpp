#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <set>

#ifdef _WIN32
#include <windows.h>
#endif

#include "Simulation.h"
#include "SimulationLib.h"

extern double particle_radius;
extern double delta_0;

enum GRID_TYPE {UNKNOWN_GRID, SIMPLE_CUBIC_GRID, HEXAGONAL_GRID};

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

void buildGrid(GRID_TYPE grid_type, double x_size, double y_size, double z_size, double filling_factor, const char *sample_filename, const char *material_filename);

int main(int argc, char **argv)
{
#ifdef _WIN32
	Sleep(1100);
#endif

	srand(time(NULL));

	if(argc == 6)
	{
		Simulation sim;
		sim.loadMaterial(argv[5]);

		SimLib::initCylindricalRBD(&sim, atoi(argv[1]), atof(argv[2]), atof(argv[3]));

		SimLib::centerCMS(&sim);
		sim.saveToFile(argv[4]);
	}
	else if(argc == 8)
	{
		Simulation sim;
		sim.loadMaterial(argv[7]);

		SimLib::initRandomBallisticDeposition(&sim, atoi(argv[1]), atof(argv[2]), atof(argv[3]));

		vec3 lower_pos, upper_pos;
		sim.getEnclosingBox(&lower_pos, &upper_pos);

		upper_pos[1] -= (upper_pos[1] - lower_pos[1]) * atof(argv[4]);
		lower_pos[1] += (upper_pos[1] - lower_pos[1]) * atof(argv[5]);
		SimLib::sliceBox(&sim, lower_pos, upper_pos);

		SimLib::centerCMS(&sim);
		sim.saveToFile(argv[6]);
	}
	else if(argc == 10)
	{
		Simulation sim;
		sim.loadMaterial(argv[9]);

		BAMSelectionMethod bam_selection_method = (BAMSelectionMethod) atoi(argv[7]);

		SimLib::initBAMAggregate(&sim, NULL, atoi(argv[1]), atof(argv[5]), atof(argv[6]), bam_selection_method);

		double x_size = atof(argv[2]) * 1e-4;
		double y_size = atof(argv[3]) * 1e-4;
		double height = atof(argv[4]) * 1e-4;

		vec3 lower, upper;
		lower[0] = -0.5 * x_size;
		lower[1] = -0.5 * height;
		lower[2] = -0.5 * y_size;
		upper[0] = 0.5 * x_size;
		upper[1] = 0.5 * height;
		upper[2] = 0.5 * y_size;

		SimLib::sliceBox(&sim, lower, upper);
		SimLib::filterFragments(&sim);

		double filling_factor = sim.getBoxFillingFactor();
		double coord_number = 2.0 * (double)SimLib::getNumberOfContacts(&sim) / (double)sim.number_of_particles;

		char filename[300];
		sprintf(filename, "%s_ff%0.2lf_nc%1.2lf.dat", argv[8], filling_factor, coord_number);

		sim.saveToFile(filename);
	}
	/*else if(argc == 8)
	{
		if(strcmp(argv[1], "SC") == 0)
			buildGrid(SIMPLE_CUBIC_GRID, atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6], argv[7]);
		else if(strcmp(argv[1], "HEX") == 0)
			buildGrid(HEXAGONAL_GRID, atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6], argv[7]);
		else
			printf("ERROR: unkown grid type %s\n", argv[1]);
	}*/
	else
	{
		printf("Incorrect arguments! Use:\n--particles -radius -z_slice_factor -sample_filename -material_filename\n");
		printf(" or\n-particles -x_size -y_size -top_slice_factor -bottom_slice_factor -sample_filename -material_filename\n");
		printf(" or\n-particles -x_size -y_size -height -migration_rate1 -migration_rate2 -BAM_method (0,1,2) -sample_filename -material_filename\n");
		//printf(" or\n-grid_type (SC, HEX) -x_size -y_size -z_size -filling_factor -sample_filename -material_filename\n");
	}

    return 0;
} 

void buildGrid(GRID_TYPE grid_type, double x_size, double y_size, double z_size, double filling_factor, const char *sample_filename, const char *material_filename)
{
	// init sim
	Simulation sim;
	sim.loadMaterial(material_filename);

	int particles = 0;
	int estimated_particles = 0;

	// calculate estimated number of particles 
	if(grid_type == SIMPLE_CUBIC_GRID)
		estimated_particles = (int) (x_size*y_size*z_size  / (8.0 * particle_radius*particle_radius*particle_radius) );

	double *positions = new double[4 * (int)(1.25 * (double)estimated_particles)];
	double dist = 2.0 * particle_radius - delta_0;
	double x_pos = 0.5 * dist;
	double y_pos = 0.5 * dist;
	double z_pos = 0.5 * dist;
	bool stop = false;

	while(!stop)
	{
		if( filling_factor > (double)rand()/(double)RAND_MAX)
		{
			positions[4*particles] = x_pos;
			positions[4*particles+1] = z_pos;
			positions[4*particles+2] = y_pos;
			positions[4*particles+3] = 0;

			++particles;
		}

		x_pos += dist;

		if(x_pos > x_size - 0.5 * dist)
		{
			x_pos = 0.5 * dist;
			y_pos += dist;
		}

		if(y_pos > y_size - 0.5 * dist)
		{
			y_pos = 0.5 * dist;
			z_pos += dist;
		}

		if(z_pos > z_size - 0.5 * dist)
			stop = true;
	}

	// copy to sim
	sim.resizeArrays(particles, 0);

	memcpy(sim.pos_old, positions, 4 * particles * sizeof(double) );

	SimLib::centerCMS(&sim);
	sim.saveToFile(sample_filename);

	delete [] positions;
}
