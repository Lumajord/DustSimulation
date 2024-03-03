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
	
        if(argc < 3)
	{
                printf("Wrong parameters! Use:\n-result_file -material_file -file1 -file2 -...\n");
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
	// open output files
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        FILE *result_file = fopen(argv[result_file_index], "a+");

        if(!result_file)
	{
                printf("ERROR:\nCannot open result file %s,\n", argv[result_file_index]);
		return EXIT_SUCCESS;
	}

        fprintf(result_file, "# filename    collision_result   collision_speed  num_particles   num_fragments   contacts_broken contacts_created    N_frag1 N_frag2\n");

        fclose(result_file);


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// analyze files
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> fragment_ids;							// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;						// number of particles of the fragment
	std::vector< std::list<int> > particles_of_fragment;	// ids of the particles of a specific fragment

        for(int i = 3; i < argc; ++i)
        {
                result_file = fopen(argv[result_file_index], "a");
		error_code = sim.loadFromFile(argv[i]);

                if(error_code != EC_OK)
                {
                    char message[200];
                    sim.getErrorMessage(error_code, message);
                    printf("File not found! %s  %s\n", argv[i], message);
                    continue;
                }

		if(error_code == EC_OK)
		{
                        printf("loading file: %s\n", argv[i]);
			fragment_ids.clear();
			size_of_fragment.clear();
			particles_of_fragment.clear();

			SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

            int fragments = size_of_fragment.size();

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
			}

                        if(largest_fragment != second_largest_fragment && particles_of_fragment.size() > 1)
			{
				CollisionResult collision_result = SimLib::getCollisionResult(particles_of_fragment[largest_fragment].size(), particles_of_fragment[second_largest_fragment].size(), sim.number_of_particles);

                                double scattering_angle = 0.0;
                                if(collision_result == COLLISION_RESULT_FRAGMENTATION || collision_result == COLLISION_RESULT_BOUNCING && true)
                                {
                                    double vx1 = 0.0;
                                    double vy1 = 0.0;

                                    double vx2 = 0.0;
                                    double vy2 = 0.0;


                                    //for(int i = 0; i < particles_of_fragment[largest_fragment].size(); ++i)
                                    for(std::list<int>::iterator it = particles_of_fragment[largest_fragment].begin(); it !=  particles_of_fragment[largest_fragment].end(); ++it)
                                    {
                                        //int p = particles_of_fragment[largest_fragment][i];
                                        int p = *it;

                                        vx1 += sim.vel[X_COORD(p)];
                                        vy1 += sim.vel[Y_COORD(p)];
                                    }

                                    //for(int i = 0; i < particles_of_fragment[largest_fragment].size(); ++i)
                                    for(std::list<int>::iterator it = particles_of_fragment[second_largest_fragment].begin(); it !=  particles_of_fragment[second_largest_fragment].end(); ++it)
                                    {
                                        //int p = particles_of_fragment[largest_fragment][i];
                                        int p = *it;

                                        vx2 += sim.vel[X_COORD(p)];
                                        vy2 += sim.vel[Y_COORD(p)];
                                    }
                                    /*
                                    for(int i = 0; i < particles_of_fragment[second_largest_fragment].size(); ++i)
                                    {
                                        int p = particles_of_fragment[second_largest_fragment][i];

                                        vx2 += sim.vel[X_COORD(p)];
                                        vy2 += sim.vel[Y_COORD(p)];
                                    }
                                    */


                                    vx1 /= particles_of_fragment[largest_fragment].size();
                                    vy1 /= particles_of_fragment[largest_fragment].size();

                                    printf("vx = %f vy = %f\n", vx1, vy1);

                                    vx2 /= particles_of_fragment[second_largest_fragment].size();
                                    vy2 /= particles_of_fragment[second_largest_fragment].size();

                                    vx1 -= vx2;
                                    vy1 -= vy2;

                                    if(vx1 < 0.0)
                                    {
                                        vx1 *= -1.0;
                                        vy1 *= -1.0;
                                    }

                                    scattering_angle = atan2 (vy1,vx1) * 180.0 / M_PI;

                                }

                                fprintf(result_file, "%s    %d  %e  %d  %d  %d  %d  %d  %d  %f\n",
                                        argv[i],
                                        collision_result,
                                        sim.sim_info.info_storage[2], // = collision speed
                                        sim.number_of_particles,
                                        fragments,
                                        sim.broken_contacts,
                                        sim.created_contacts,
                                        int(particles_of_fragment[largest_fragment].size()),
                                        int(particles_of_fragment[second_largest_fragment].size()),
                                        scattering_angle);

			}
			else
                        {
                            fprintf(result_file, "%s    %d  %e  %d  %d  %d  %d  %d  %d  %f\n",
                                    argv[i],
                                    COLLISION_RESULT_STICKING,
                                    sim.sim_info.info_storage[2], // = collision speed
                                    sim.number_of_particles,
                                    fragments,
                                    sim.broken_contacts,
                                    sim.created_contacts,
                                    int(particles_of_fragment[largest_fragment].size()),
                                    0, 0.0);
                        }
		}
            else
                {
                    printf("Could not load collision file: %s   because: %s\n", argv[i], error_code);
                }


                fclose(result_file);

	}



    printf("Exit success\n");
    return EXIT_SUCCESS;
} 
