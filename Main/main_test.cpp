#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include "cutil.h"
#include <omp.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;
extern double ENERGY_UNIT;

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	#undef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
#endif

#ifdef TRACK_DISSIPATED_ENERGY
	#undef TRACK_DISSIPATED_ENERGY
#endif


int main(int argc, char **argv)
{
	/*int runs = 10;
	unsigned int timer;
	float avg_time = 0;

	srand(time(NULL));

	Simulation sim;
	sim.loadMaterial("Silicate.mat");

	omp_set_num_threads(4);

	#pragma omp parallel
	{
		#pragma omp master
		printf("Using %i threads...\n", omp_get_num_threads());
	}

	for(int run = 0; run < runs; ++run)
	{
		sim.loadFromFile("test_collision.dat");
		sim.startSimulation(1e-4, 4e-10);

		printf("Run: %i", run+1);
	
		cutCreateTimer(&timer);
		cutStartTimer(timer);

		for(int i = 0; i < 10000; ++i)
			sim.update();

		cutStopTimer(timer);
		float time = cutGetTimerValue(timer);
		printf("  %f ms\n", time);
		avg_time += time;
		cutDeleteTimer(timer);
	}

	sim.saveToFile("result.dat");

	printf("\nSpeed test finished - avg run time: %f ms\n", avg_time / (double)runs);

	system("pause");
    return 0;*/

	//sleep(1);
	srand(time(NULL));

	int num;

	FILE *file = fopen("agglomerates.dat", "r");
	fscanf(file, "%i", &num);
	fclose(file);

	file = fopen("collision.dat", "w+");
	fprintf(file, "%i %i %lf %lf", rand()%num, rand()%num, 5.0 + 25 * (double)rand()/(double)(RAND_MAX+1), 0.5 * (double)rand()/(double)(RAND_MAX+1)); 
	fclose(file);
}
