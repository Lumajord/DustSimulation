#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>

#include "Simulation.h"
#include "SimulationLib.h"

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
  #undef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
#endif

#ifdef TRACK_DISSIPATED_ENERGY
  #undef TRACK_DISSIPATED_ENERGY
#endif

bool doesFileExist(const char *filename) 
{
  struct stat stFileInfo;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(filename, &stFileInfo);
  
  if(intStat == 0)
	  return true;
  else 
	  return false;
}

void collideChains(Simulation *sim, int runs, double sim_time, double timestep, int chain_base_size, int chain_random_size, int min_impact_speed, int max_impact_speed, char *dest_folder, char *log_filename)
{
	char result_filename[200];
	double impact_speed;
	int projectile_size;
	int chain_size;

	printf("Run		agglomerate size	impact speed	projectile size\n\n");

	for(int run = 1; run < runs; ++run)
	{
		chain_size = chain_base_size + rand()%chain_random_size;
		projectile_size = chain_base_size + rand()%chain_random_size;

		impact_speed = min_impact_speed + (max_impact_speed - min_impact_speed)/100.0 * (double)(rand()%101);

		SimLib::initChain(sim, chain_size, projectile_size, rand()%(chain_size+1), 0.45, impact_speed);
		sim->startSimulation(sim_time, timestep);

		while(sim->current_time < sim->end_time)
			sim->update();

		SimLib::removeFragments(sim);
		SimLib::centerCMS(sim);

		// save to next available file
		if(sim->number_of_particles > 2)
		{
			char c = 'a';

			while(c <= 'z')
			{
				sprintf(result_filename, "%s\\agglomerate%i%c.dat", dest_folder, sim->number_of_particles, c);

				if(!doesFileExist(result_filename))
					break;
				else
					++c;
			}

			sim->saveToFile(result_filename);

			// add entry to agglomerate list
			FILE *log_file = fopen(log_filename, "a+");
			int counter = 0;

			// if file is currently inaccessible (opened by another application) try later
			while(!log_file && counter < 100000)
			{
				++counter;
				log_file = fopen(log_filename, "a+");
			}

			if(log_file)
			{
				fprintf(log_file, "%s %s %s %lf\n", "nofile", "nofile", result_filename, impact_speed);
				fclose(log_file);
			}
		}


		printf("%i		%i		%lf		%i\n", run, sim->number_of_particles, impact_speed, projectile_size);
	}
}

bool getSampleFile(char *filename, char *foldername, int size)
{
	char number_of_samples;
	char c = 'a';

	// determine how many different samples are available
	while(c <= 'z')
	{
		sprintf(filename, "%s\\agglomerate%i%c.dat", foldername, size, c);

		if(doesFileExist(filename))
			++c;
		else
			break;
	}

	number_of_samples = c - 'a';

	if(number_of_samples > 0)
	{
		// select one sample
		sprintf(filename, "%s\\agglomerate%i%c.dat", foldername, size, 'a' + rand()%number_of_samples);

		return true;
	}
	else
		return false;
}

void performCollision(Simulation *sim, double sim_time, double timestep, char *filename1, char *filename2, double impact_speed, int run, int repetitions, char *dest_folder, char *log_filename)
{
	char result_filename[200];

	for(int r = 0; r < repetitions; ++r)
	{
		SimLib::collideTwoAgglomerates(sim, filename1, filename2, impact_speed, 0, true);
		sim->startSimulation(sim_time, timestep);

		printf("Collding %s with %s, impact speed %i, run (%i : %i)...", filename1, filename2, (int)impact_speed, run, r+1);

		// prepare progress bar
		int log_counter = 0;
		double percentage = 0;
		printf("0.00%% done");

		while(sim->current_time < sim->end_time)
		{
			sim->update();
			++log_counter;

			if(log_counter == 100)
			{
				log_counter = 0;
			
				if(percentage < 10.0)
					printf("\b\b\b\b\b\b\b\b\b\b");
				else
					printf("\b\b\b\b\b\b\b\b\b\b\b");

				percentage = 100.0 * sim->current_time/sim->end_time;
				printf("%.2lf%% done", percentage);
			}
		}

		SimLib::removeFragments(sim);
		SimLib::centerCMS(sim);

		// save to next available file
		if(sim->number_of_particles > 5)
		{
			char c = 'a';

			while(c <= 'z')
			{
				sprintf(result_filename, "%s\\agglomerate%i%c.dat", dest_folder, sim->number_of_particles, c);

				if(!doesFileExist(result_filename))
					break;
				else
					++c;
			}

			sim->saveToFile(result_filename);

			// add entry to agglomerate list
			FILE *log_file = fopen(log_filename, "a+");
			int counter = 0;

			// if file is currently inaccessible (opened by another application) try later
			while(!log_file && counter < 100000)
			{
				++counter;
				log_file = fopen(log_filename, "a+");
			}

			if(log_file)
			{
				fprintf(log_file, "%s %s %s %lf\n", filename1, filename2, result_filename, impact_speed);
				fclose(log_file);
			}
		}

		printf(" ...done\n", sim->number_of_particles);
	}
}

void collideAgglomerates(Simulation *sim, int runs, int repetitions, double sim_time, double timestep, int min_size1, int max_size1, int min_size2, int max_size2, int min_impact_speed, int max_impact_speed, char *sample_folder, char *dest_folder, char *log_filename)
{
	
	double impact_speed;
	char filename1[200];
	char filename2[200];
	
	int counter;

	bool different_sample = true;//false;

	if(!different_sample)
	{
		counter = 0;

		while(counter < 50)
		{
			if( getSampleFile(filename1, sample_folder, min_size1 + rand()%(max_size1-min_size1+1)) )
				break;
			else
				++counter;
		}
	}
	
	for(int run = 0; run < runs; ++run)
	{
		// select sample files
		if(different_sample)
		{
			counter = 0;

			while(counter < 30)
			{
				if( getSampleFile(filename1, sample_folder, min_size1 + rand()%(max_size1-min_size1+1)) )
					break;
				else
					++counter;
			}
		}

		counter = 0;

		while(counter < 30)
		{
			if( getSampleFile(filename2, sample_folder, min_size2 + rand()%(max_size2-min_size2+1)) )
				break;
			else
				++counter;
		}

		if( doesFileExist(filename1) && doesFileExist(filename2) )
		{
			impact_speed = min_impact_speed + (max_impact_speed - min_impact_speed)/100.0 * (double)(rand()%101);

			performCollision(sim, sim_time, timestep, filename1, filename2, impact_speed, run, repetitions, dest_folder, log_filename);
		}
	}
}

int main(int argc, char **argv)
{
	srand(time(NULL));

	int runs = 10; 
	int repetitions = 1; 
	int min_size1 = 10; 
	int max_size1 = 100; 
	int min_size2 = 250; 
	int max_size2 = 300;
	int min_impact_speed = 500;
	int max_impact_speed = 1500;

	double angular_irregularity = 0;
	double sim_time = 0;
	double timestep = 0;

	char dest_folder[200];
	char sample_folder[200];
	char filename1[200];
	char filename2[200];
	char material_filename[200];

	bool collide_predefined_agglomerates = false;
	bool collide_chains = false;

	if(argc == 14)
	{
		runs = atoi(argv[1]);
		repetitions = atoi(argv[2]);
		sim_time = atof(argv[3]);
		timestep = atof(argv[4]);
		min_size1 = atoi(argv[5]); 
		max_size1 = atoi(argv[6]); 
		min_size2 = atoi(argv[7]); 
		max_size2 = atoi(argv[8]);
		min_impact_speed = atoi(argv[9]);
		max_impact_speed = atoi(argv[10]);
		sprintf(sample_folder, "%s", argv[11]);
		sprintf(dest_folder, "%s", argv[12]);
		sprintf(material_filename, "%s", argv[13]);
	}
	else if(argc == 10)
	{
		collide_chains = true;

		runs = atoi(argv[1]);
		sim_time = atof(argv[2]);
		timestep = atof(argv[3]);
		min_size1 = atoi(argv[4]); 
		max_size1 = atoi(argv[5]); 
		min_impact_speed = atoi(argv[6]);
		max_impact_speed = atoi(argv[7]);
		sprintf(dest_folder, "%s", argv[8]);
		sprintf(material_filename, "%s", argv[9]);
	}
	else if(argc == 8)
	{
		collide_predefined_agglomerates = true;

		runs = atoi(argv[1]);
		sim_time = atof(argv[2]);
		timestep = atof(argv[3]);
		sprintf(filename1, "%s", argv[4]);
		sprintf(filename2, "%s", argv[5]);
		min_impact_speed = atoi(argv[6]);
		max_impact_speed = atoi(argv[6]);
		sprintf(dest_folder, "%s", argv[7]);
		sprintf(material_filename, "%s", argv[8]);
	}
	else
	{
		printf("Incorrect arguments- use:\n-runs -repetitions -sim_time -timestep -min_size1 -max_size1 -min_size2 -max_size2 -min_impact_speed -max_impact_speed -sample_folder -dest_folder -material_filename\n");
		printf("or\n-runs -sim_time -timestep -min_chain_size -max_chain_size -min_impact_speed -max_impact_speed -dest_folder -material_filename\n");
		printf("or\n-runs -sim_time -tiemstep -filename1 -filename2 -impact_speed -dest_folder -material_filename\n");
		return 0;
	}

	printf("Running agglomerate builder using particle simulation core v%s\n", CORE_VERSION);

	char result_filename[300];
	sprintf(result_filename, "%s\\agglomerate_list.txt", dest_folder);

	// load material file
	Simulation sim;
	
	if(!sim.loadMaterial(material_filename))
	{
		printf("Error: Could not load material from file %s\n", material_filename);
		return 0;
	}

	if(collide_chains)
		collideChains(&sim, runs, sim_time, timestep, min_size1, max_size1 - min_size1, min_impact_speed, max_impact_speed, dest_folder, result_filename);
	else
	{
		if(!collide_predefined_agglomerates)
			collideAgglomerates(&sim, runs, repetitions, sim_time, timestep, min_size1, max_size1, min_size2, max_size2, min_impact_speed, max_impact_speed, sample_folder, dest_folder, result_filename);
		else
			performCollision(&sim, sim_time, timestep, filename1, filename2, min_impact_speed, runs, 1, dest_folder, result_filename);
	}

    return 0;
} 


