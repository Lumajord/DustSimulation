#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

#include "Simulation.h"
#include "SimulationLib.h"

extern double particle_radius;
extern double reduced_radius;
extern double density;
extern double mass;
extern double surface_energy;
extern double crit_rolling_displacement;
extern double ENERGY_UNIT;

#ifndef TRACK_DISSIPATED_ENERGY
	#error Tracking of dissipated energy not enabled!
#endif

#ifdef ENABLE_GPU
	#include "SimulationCuda.h"
#endif

const double critical_sticking_velocity = 40.0;

void performCollision(double timestep, double sim_time, double damping_time, const char *material_filename, const char *agglomerates_path, double hit_and_stick_modifier, const char *log_filename = NULL)
{
	int id1, id2;
	double impact_velocity;
	double impact_parameter;
	double collision_time;	// not used for collision > stored in new agglomerate table for mc code
	long int seed;			// not used for collision > stored in new agglomerate table for mc code
	int number_of_collisions; // number of collisions that should be performed in one simulations (smaller aggregate will be used several times)

	printf("Performing collision:\n timestep=%g, sim_time=%g, damping_time=%g, material_filename=%s, agg_path=%s, log_filename=%s\n", timestep, sim_time, damping_time, material_filename, agglomerates_path, log_filename);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// get collision data
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *file = fopen("collision.dat", "r");

	if(!file)
	{
		printf("ERROR: Could not access file collision.dat!");
		return;
	}
	else
	{
		double temp;
		fscanf(file, "%i %i %lf %lf %li %lf %lf", &id1, &id2, &impact_velocity, &impact_parameter, &seed, &collision_time, &temp);
		fclose(file);

		number_of_collisions = (int)temp;

		/*if(impact_velocity > 1000)
		{
			FILE *error_file = fopen("error.msg", "w+");
			if(error_file)
			{
				fprintf(error_file, "ERROR:\nImpact velocity too large %g\nNumber of collisions: %i\n", impact_velocity, number_of_collisions);
				fclose(error_file);
			}
			return;
		}*/
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load material
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_GPU
	SimulationCuda sim;
#else
	Simulation sim;
#endif

	ErrorCode error_code = sim.loadMaterial(material_filename);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", material_filename, message);

		FILE *error_file = fopen("error.msg", "w+");
		if(error_file)
		{
			fprintf(error_file, "ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", material_filename, message);
			fclose(error_file);
		}

		return;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load aggregates
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	char filename1[300];
	char filename2[300];
	char buffer[50];

	strcpy(filename1, agglomerates_path);
	sprintf(buffer, "agglomerate_%i.dat", id1);
	strcat(filename1, buffer);

	strcpy(filename2, agglomerates_path);
	sprintf(buffer, "agglomerate_%i.dat", id2);
	strcat(filename2, buffer);

	printf("Colliding %s with %s...", filename1, filename2);

	// get number of particles of first agglomerate (to select remaining fragment after the collisions)
	error_code = sim.loadFromFile(filename1);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to determine the number of monomeres of the first aggregate, the following error occurred:\n%s\n", message);

		FILE *error_file = fopen("error.msg", "w+");
		if(error_file)
		{
			fprintf(error_file, "ERROR:\nWhile trying to determine the size of the first aggregate, the following error occurred:\n%s\n", message);
			fclose(error_file);
		}

		return;
	}

	int number_of_particles1 = sim.number_of_particles;

	// get size of second agglomerate (to be able to determine minimum sim time)
	error_code = sim.loadFromFile(filename2);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to determine the size of the second aggregate, the following error occurred:\n%s\n", message);

		FILE *error_file = fopen("error.msg", "w+");
		if(error_file)
		{
			fprintf(error_file, "ERROR:\nWhile trying to determine the size of the second aggregate, the following error occurred:\n%s\n", message);
			fclose(error_file);
		}

		return;
	}

	int number_of_particles2 = sim.number_of_particles;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// check hit & stick
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool hit_and_stick = false;
	bool gpu_sim = false;

	double v1 = impact_velocity / (1.0 + (double)number_of_particles1 / (double)number_of_particles2); 
	double v2 = impact_velocity / (1.0 + (double)number_of_particles2 / (double)number_of_particles1);

	double E_kin = 0.5 * mass * ((double)number_of_particles1) * v1 * v1 + 0.5 * mass * ((double)number_of_particles2) * v2 * v2;
	double E_roll = 6.0 * M_PI*M_PI * surface_energy * reduced_radius * crit_rolling_displacement;

	if( hit_and_stick_modifier * E_roll > E_kin && impact_velocity < critical_sticking_velocity)	
	{
		error_code = SimLib::hitAndStick(&sim, filename1, filename2, impact_parameter, true, true);
		
		for(int collisions = 1; collisions < number_of_collisions; ++collisions)
			error_code = SimLib::hitAndStick(&sim, NULL, filename2, impact_parameter, true, true);

		if(error_code != EC_OK)
		{
			char message[200];
			sim.getErrorMessage(error_code, message);
			printf("ERROR - Hit & stick collision failed:\n%s\n", message);

			FILE *error_file = fopen("error.msg", "w+");
			if(error_file)
			{
				fprintf(error_file, "ERROR - Hit & stick collision failed:\n%s\n", message);
				fclose(error_file);
			}

			return;
		}

		hit_and_stick = true;
		goto collision_finished;
	}
	else
	{
		if(number_of_collisions > 1)
		{
			//error_code = SimLib::collideTwoAgglomerates(&sim, filename1, filename2, impact_velocity, impact_parameter, 0, true, true, true, true);
			error_code = SimLib::impactMultiProjectiles(&sim, filename1, true, filename2, number_of_collisions, 0, impact_parameter, impact_parameter, false, 0.001 * particle_radius, 0.001 * particle_radius, true);
			//error_code = SimLib::collideMultipleAgglomerates(&sim, filename1, filename2, number_of_collisions, impact_velocity, impact_parameter, 0, true, true, true, true);
		}
		else
			error_code = SimLib::collideTwoAgglomerates(&sim, filename1, filename2, impact_velocity, impact_parameter, 0, true, true, true, true);
	}

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to setup the collision, the following error occurred:\n%s\n", message);

		FILE *error_file = fopen("error.msg", "w+");
		if(error_file)
		{
			fprintf(error_file, "ERROR:\nWhile trying to setup the collision, the following error occurred:\n%s\n", message);
			fclose(error_file);
		}

		return;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// run simulation
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	{
	double min_time = std::max(1e-5, 6.0 * particle_radius / impact_velocity);
	int check_interval = (int) (2e-6 / timestep);
	sim.setPotentialVariationStopCondition(min_time, 0.5, check_interval);

#ifdef ENABLE_GPU
	if(sim.number_of_particles > 1500)
	{
		error_code = sim.initCuda();

		if(error_code == EC_OK)
		{
			sim.toggleGPUMode(true);
			sim.startSimulation(sim_time, timestep);

			while(!sim.stop_simulation)
				sim.update();

			// damp
			sim.check_potential_variation_interval = 0;
			sim.startSimulation(damping_time, timestep);
			sim.setDampingFactor(0.995);

			while(!sim.stop_simulation)
				sim.update();

			sim.copySimDataFromGPU();
			gpu_sim = true;
			goto collision_finished;
		}
		else
			printf("ERROR: Failed to init cuda simulation\n");
	}
#endif

	sim.startSimulation(sim_time, timestep);

	while(!sim.stop_simulation)
		sim.update();

	// damp
	sim.startSimulation(damping_time, timestep);

	while(!sim.stop_simulation)
	{
		sim.update();
		sim.dampVelocities(0.995);
	}

	memset(sim.vel, 0, sim.number_of_particles * 3 * sizeof(double));
	memset(sim.vel_angular, 0, sim.number_of_particles * 3 * sizeof(double));
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// save result & update table
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

collision_finished:

	// determine fragment that contains monomere with random id
	int id = rand()%number_of_particles1;
	std::vector<int> fragment_ids;				// array storing the fragment id of every particle
	std::vector<int> size_of_fragment;			// number of particles of the fragment

	SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, NULL);
	SimLib::removeFragments(&sim, fragment_ids[id], fragment_ids, size_of_fragment);

	SimLib::centerCMS(&sim);
	sim.saveToFile(filename1, true);

	double cross_section;
	double outer_radius;
	double gyration_radius;

	SimLib::getCrossSection(sim, particle_radius / 5.0, 8, &cross_section, NULL);
	SimLib::getSize(sim, &gyration_radius, &outer_radius);

	if(cross_section < 1e-9)
	{
		sim.saveToFile("cross_section_error.dat");
		printf("ERROR: cross section failed!\n");
	}

	printf("... done. New size: %i\n", sim.number_of_particles);

	// open agglomerate table
	file = fopen("agglomerates.dat", "r");

	if(!file)
	{
		printf("ERROR: Could not open agglomerates.dat!\n");
		return;
	}
	else
	{
		// copy old file
		int number_of_agglomerates;
		double d_dummy;
		long int i_dummy;
		fscanf(file, "%i", &number_of_agglomerates);
		fscanf(file, "%lf", &d_dummy);
		fscanf(file, "%lf %lf %lf", &d_dummy, &d_dummy, &d_dummy);
		fscanf(file, "%li", &i_dummy);

		double *cross_sections = new double[number_of_agglomerates];
		double *outer_radii = new double[number_of_agglomerates];
		double *gyration_radii = new double[number_of_agglomerates];
		int *monomeres = new int[number_of_agglomerates];

		for(int i = 0; i < number_of_agglomerates; ++i)
			fscanf(file, "%lf %lf %lf %lf %i", &d_dummy, &(cross_sections[i]), &(outer_radii[i]), &(gyration_radii[i]), &(monomeres[i]));

		fclose(file);

		// log result
		if(log_filename)
		{
			FILE *log_file = fopen(log_filename, "a+");

			if(log_file)
			{
				int fragment_particles = number_of_particles1 + number_of_collisions * number_of_particles2 - sim.number_of_particles;

				if(hit_and_stick)
					fprintf(log_file, "%i %i %lf -> %i %i  hit and stick", monomeres[id1-1], monomeres[id2-1], impact_velocity, sim.number_of_particles, fragment_particles);
				else if(gpu_sim)
					fprintf(log_file, "%i %i %lf -> %i %i  ballistic gpu", monomeres[id1-1], monomeres[id2-1], impact_velocity, sim.number_of_particles, fragment_particles);
				else 
					fprintf(log_file, "%i %i %lf -> %i %i  ballistic cpu", monomeres[id1-1], monomeres[id2-1], impact_velocity, sim.number_of_particles, fragment_particles);

				if(number_of_collisions > 1)
					fprintf(log_file, "  multi collision: %i", number_of_collisions);

				fprintf(log_file, "\n");
				fclose(log_file);
			}
		}

		// replace line
		cross_sections[id1-1] = cross_section;
		outer_radii[id1-1] = outer_radius;
		gyration_radii[id1-1] = gyration_radius + particle_radius;
		monomeres[id1-1] = sim.number_of_particles;

		// write new file
		file = fopen("agglomerates.dat", "w+");
		fprintf(file, "%i\n", number_of_agglomerates);
		fprintf(file, "%g\n", collision_time);
		fprintf(file, "%g %g %g\n", density, mass, particle_radius);
		fprintf(file, "%li\n", seed);

		for(int i = 0; i < number_of_agglomerates; ++i)
		{
			fprintf(file, "%g %g %g %g %i\n", (double)monomeres[i] * mass, cross_sections[i], outer_radii[i], gyration_radii[i], monomeres[i]);

			if(cross_sections[i] < 1e-9)
			{	
				FILE *error_file = fopen("error.msg", "w+");
				if(error_file)
				{
					fprintf(error_file, "ERROR:\nCross section of agg %i too small: %g\nCollision between %i and %i -> monomers: %i cs: %g", i, cross_sections[i], id1, id2, sim.number_of_particles, cross_section);
					fclose(error_file);
				}
			}
		}

		fclose(file);

		delete [] cross_sections;
		delete [] gyration_radii;
		delete [] outer_radii;
		delete [] monomeres;
	}
}

void initEnsemble(int number_of_agglomerates, const char *material_filename, const char *agglomerates_path)
{
	// load material
	Simulation sim;
	ErrorCode error_code = sim.loadMaterial(material_filename);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s\n", material_filename, message);
		return;
	}

	// setup sim with 1 particle
	SimLib::initChain(&sim, 1, 0, 0, 0, 0);
	SimLib::centerCMS(&sim);

	// create file storing agglomerate table
	FILE *file = fopen("agglomerates.dat", "w+");

	if(!file)
	{
		printf("ERROR: Could not open agglomerates.dat!");
		return;
	}

	// generate random seed
	srand(time(NULL));
	long int seed = rand()%900000000;

	fprintf(file, "%i\n", number_of_agglomerates);
	fprintf(file, "%g\n", 0);
	fprintf(file, "%g %g %g\n", density, mass, particle_radius);
	fprintf(file, "%li\n", seed);

	char filename[300];
	char buffer[50];

	for(int id = 1; id <= number_of_agglomerates; ++id)
	{
		// add entry to agglomerate table
		fprintf(file, "%g %g %g %g %i\n", mass, M_PI * particle_radius*particle_radius, particle_radius, particle_radius, 1);

		// generate agglomerate
		strcpy(filename, agglomerates_path);
		sprintf(buffer, "agglomerate_%i.dat", id);
		strcat(filename, buffer);

		sim.saveToFile(filename);
	}

	fclose(file);
}

int main(int argc, char **argv)
{
	srand(time(NULL));

	if (argc == 3)
	{
		initEnsemble(atoi(argv[1]), argv[2], NULL);
	}
	else if(argc == 4)
	{
		initEnsemble(atoi(argv[1]), argv[2], argv[3]);
	}
	else if(argc == 6)
	{
		performCollision(atof(argv[1]), atof(argv[2]), atof(argv[3]), argv[4], argv[5], 0, NULL);
	}
	else if(argc == 7)
	{
		performCollision(atof(argv[1]), atof(argv[2]), atof(argv[3]), argv[4], argv[5], atof(argv[6]), NULL);
	}
	else if(argc == 8)
	{
		performCollision(atof(argv[1]), atof(argv[2]), atof(argv[3]), argv[4], argv[5], atof(argv[6]), argv[7]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n -timestep -sim_time -damping_time -material_filename -agglomerates_path\nor\n -timestep -sim_time -damping_time -material_filename -agglomerates_path -hit_and_stick_modifier\nor\n -timestep -sim_time -damping_time -material_filename -agglomerates_path -hit_and_stick_modifier -log_file\n");
	}

    return EXIT_SUCCESS;
}
