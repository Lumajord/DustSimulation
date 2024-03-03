#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include <string.h>

#include "Constants.h"
#include "Simulation.h"
#include "SimulationLib.h"

extern double mass;

int main(int argc, char **argv)
{
	int result_file_index = 1;
	int material_file_index = 2;
	double bouncing_threshold = 0.1;
	
	if(argc >= 5)
	{
		bouncing_threshold = atof(argv[3]);
	}
	else
	{
		printf("Wrong parameters! Use:\n-result_file -material_file -bouncing_threshold -file1 -file2 -...\n");
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", argv[material_file_index], message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepare result file
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *file = fopen(argv[result_file_index], "w+");

	if(!file)
	{
		printf("ERROR:\nCannot open result file %s,\n", argv[result_file_index]);
		return EXIT_SUCCESS;
	}

	fprintf(file, "# collision velocity (cm/s)     coefficient of restitution     coefficient of restitution (with rotation)\n");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// analyze files
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> fragment_ids;							// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;						// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment

	for(int i = 4; i < argc; ++i)
	{
		error_code = sim.loadFromFile(argv[i]);

		if(error_code == EC_OK)
		{
			fragment_ids.clear();
			size_of_fragment.clear();
			particles_of_fragment.clear();

			SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

			int fragments = size_of_fragment.size();
			double *fragment_velocities = new double[3*fragments];
			memset(fragment_velocities, 0, 3 * fragments * sizeof(double));

			// determine largest and second largest fragment
			int largest_fragment = 0;
			int second_largest_fragment = 0;

			for(int f = 0; f < fragments; ++f)
			{
				int fragment_size = particles_of_fragment[f].size();

				if(fragment_size > particles_of_fragment[largest_fragment].size())
				{
					second_largest_fragment = largest_fragment;
					largest_fragment = f;
				}
				else
				{
					if(fragment_size > particles_of_fragment[second_largest_fragment].size() || second_largest_fragment == largest_fragment)
						second_largest_fragment = f;
				}

				// calculate velocity of fragment
				for(std::list<int>::iterator p = particles_of_fragment[f].begin(); p != particles_of_fragment[f].end(); ++p)
				{
					fragment_velocities[3*f+0] += sim.vel[X_COORD(*p)];
					fragment_velocities[3*f+1] += sim.vel[Y_COORD(*p)];
					fragment_velocities[3*f+2] += sim.vel[Z_COORD(*p)];
				}

				fragment_velocities[3*f+0] /= (double)fragment_size;
				fragment_velocities[3*f+1] /= (double)fragment_size;
				fragment_velocities[3*f+2] /= (double)fragment_size;
			}

			if(largest_fragment != second_largest_fragment && sim.sim_info.info_storage[2] > 0)
			{
				vec3 delta_v;
				delta_v[0] = fragment_velocities[3*largest_fragment+0] - fragment_velocities[3*second_largest_fragment+0];
				delta_v[1] = fragment_velocities[3*largest_fragment+1] - fragment_velocities[3*second_largest_fragment+1];
				delta_v[2] = fragment_velocities[3*largest_fragment+2] - fragment_velocities[3*second_largest_fragment+2];
				double speed = sqrt(delta_v[0]*delta_v[0] + delta_v[1]*delta_v[1] + delta_v[2]*delta_v[2]);
				
				if(particles_of_fragment[largest_fragment].size() - particles_of_fragment[second_largest_fragment].size() < (int)(bouncing_threshold * (double)sim.number_of_particles))
				{
					double E_kin_init = (double)sim.number_of_particles * 0.25 * mass * sim.sim_info.info_storage[2]*sim.sim_info.info_storage[2];
					double E_kin_final = 0.5 * (0.5 * (double)(particles_of_fragment[largest_fragment].size() + particles_of_fragment[second_largest_fragment].size()) * mass) * speed*speed;

					double E_tot = 0;
					for(std::list<int>::iterator p = particles_of_fragment[largest_fragment].begin(); p != particles_of_fragment[largest_fragment].end(); ++p)
						E_tot += (sim.vel[X_COORD(*p)]*sim.vel[X_COORD(*p)] + sim.vel[Y_COORD(*p)]*sim.vel[Y_COORD(*p)] + sim.vel[Z_COORD(*p)]*sim.vel[Z_COORD(*p)]);

					for(std::list<int>::iterator p = particles_of_fragment[second_largest_fragment].begin(); p != particles_of_fragment[second_largest_fragment].end(); ++p)
						E_tot += (sim.vel[X_COORD(*p)]*sim.vel[X_COORD(*p)] + sim.vel[Y_COORD(*p)]*sim.vel[Y_COORD(*p)] + sim.vel[Z_COORD(*p)]*sim.vel[Z_COORD(*p)]);

					E_tot *= 0.5 * mass;
					double E_rot = E_tot - E_kin_final;
					double res_coeff = sqrt( (E_kin_final + E_rot) / E_kin_init );

					fprintf(file, "%g %g %g\n", sim.sim_info.info_storage[2], speed / sim.sim_info.info_storage[2], res_coeff);
				}
			}

			delete [] fragment_velocities;
		}
	}

    return EXIT_SUCCESS;
} 
