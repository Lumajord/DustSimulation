#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <vector>
#include <list>
#include <math.h>
#include "Simulation.h"


int main(int argc, char **argv)
{
	if(argc != 3)
	{
		printf("ERROR: Wrong arguments!\nUse: -input_file -output_file\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	// open input & output files
	/////////////////////////////////////////////////////////////////////////////////////////////

	FILE *data_file = fopen(argv[1], "rb");

	if(!data_file)
	{
		printf("ERROR: Cannot open %s!\n", argv[1]);
		return EXIT_SUCCESS;
	}

	FILE *output_file = fopen(argv[2], "wb+");

	if(!output_file)
	{
		fclose(data_file);
		printf("ERROR: Cannot open %s!\n", argv[2]);
		return EXIT_SUCCESS;
	}


	/////////////////////////////////////////////////////////////////////////////////////////////
	// read input file
	/////////////////////////////////////////////////////////////////////////////////////////////

	char version_buffer[23];
	fread(version_buffer, sizeof(char), 22, data_file);
	version_buffer[22] = 0;

	if(strcmp(version_buffer, "DATA_FILE_VERSION_1_09") != 0)
	{
		printf("ERROR: File version %s not supported!\n", version_buffer);
		fclose(data_file);
		fclose(output_file);
		return EXIT_SUCCESS;
	}

	int num_particles, num_walls, num_contacts;
	WallBox box;
	DataFileOffsets offsets;
	SimInfo info;

	fread(&offsets, sizeof(DataFileOffsets), 1, data_file);
	fread(&num_particles, sizeof(int), 1, data_file);
	fread(&num_walls, sizeof(int), 1, data_file);
	fread(&num_contacts, sizeof(int), 1, data_file);
	fread(&info, sizeof(SimInfo), 1, data_file);

	// load box info
	if(info.sim_type >= SIM_TYPE_COMPRESSION_NO_SIDE_WALLS && info.sim_type <= SIM_TYPE_SHOCKWAVE)
		fread(&box, sizeof(WallBox), 1, data_file);
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// load particle data
	/////////////////////////////////////////////////////////////////////////////////////////////

	double *positions = new double[3*num_particles];
	double *velocities = new double[3*num_particles];
	double *angular_velocities = new double[3*num_particles];

	fread(positions, sizeof(double), 3*num_particles, data_file);
	fread(velocities, sizeof(double), 3*num_particles, data_file);
	fread(angular_velocities, sizeof(double), 3*num_particles, data_file);

	/////////////////////////////////////////////////////////////////////////////////////////////
	// read walls
	/////////////////////////////////////////////////////////////////////////////////////////////
	
	//for(int w = 0; w < number_of_walls; ++w)
	//	fread(&(walls[w]), sizeof(Wall), 1, file);

	/////////////////////////////////////////////////////////////////////////////////////////////
	// print info
	/////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(output_file, "%i %i %i\n", num_particles, num_contacts, num_walls);

	for(int p = 0; p < num_particles; ++p)
	{
		fprintf(output_file, "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", positions[X_COORD(p)], positions[Y_COORD(p)], positions[Z_COORD(p)],
			velocities[X_COORD(p)], velocities[Y_COORD(p)], velocities[Z_COORD(p)], angular_velocities[X_COORD(p)], angular_velocities[Y_COORD(p)], angular_velocities[Z_COORD(p)]);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	// clean up
	/////////////////////////////////////////////////////////////////////////////////////////////

	delete [] positions;
	delete [] velocities;
	delete [] angular_velocities;

	fclose(data_file);
	fclose(output_file);

    return EXIT_SUCCESS;
} 
