#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <cstring>
//#include "cutil.h"

#include "Simulation.h"
#include "SimulationLib.h"

extern double F_c;

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	#undef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
#endif

#ifdef TRACK_DISSIPATED_ENERGY
	#undef TRACK_DISSIPATED_ENERGY
#endif

#ifdef AVERAGE_WALL_FORCE
	#undef AVERAGE_WALL_FORCE
#endif

//#define USE_GRAVITY

//#define LOG_MAX_SPEED

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

int main(int argc, char **argv)
{
	srand(time(NULL));

	// will be used until initial wall<->particle contacts have been initialized (to avoid disruption of agglomerate upon wall impact for high compression velocities
	double init_timestep = 0;
	double init_wall_speed = 0;	
	double init_penetration_depth = 0;

	// for main simulation
	double timestep;
	double wall_speed;
	double stop_filling_factor;

	int log_interval;
	char sample_filename[200];
	char log_filename[200];
	char material_filename[200];

	if(argc == 8)	// main mode
	{
		timestep = atof(argv[1]);
		wall_speed = atof(argv[2]);
		stop_filling_factor = atof(argv[3]);
		sprintf(sample_filename, "%s", argv[4]);
		sprintf(log_filename, "%s", argv[5]);
		log_interval = atoi(argv[6]);
		sprintf(material_filename, "%s", argv[7]);
	}
	else if(argc == 11)	// initial wall movement at different velocity
	{
		init_timestep = atof(argv[1]);
		init_wall_speed = atof(argv[2]);
		init_penetration_depth = atof(argv[3]);
		timestep = atof(argv[4]);
		wall_speed = atof(argv[5]);
		stop_filling_factor = atof(argv[6]);
		sprintf(sample_filename, "%s", argv[7]);
		sprintf(log_filename, "%s", argv[8]);
		log_interval = atoi(argv[9]);
		sprintf(material_filename, "%s", argv[10]);
	}
	else if(argc == 7)
	{
		timestep = atof(argv[1]);
		stop_filling_factor = atof(argv[2]);
		sprintf(sample_filename, "%s", argv[3]);
		sprintf(log_filename, "%s", argv[4]);
		log_interval = atoi(argv[5]);
		sprintf(material_filename, "%s", argv[6]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n--timestep -wall_speed -stop_filling_factor -sample_filename -log_filename -log_interval -material_filename\n");
		return 0;
	}

	printf("Performing aggregate compaction using particle simulation core v%s\n", CORE_VERSION);

	if(!doesFileExist(sample_filename))
	{
		printf("Error: File %s not found!\n", sample_filename);
		return 0;
	}

	// load material file
	Simulation sim;

	if(!sim.loadMaterial(material_filename))
	{
		printf("Error: Could not load material from file %s\n", material_filename);
		return 0;
	}

	// setup box
	if(argc == 11)
	{
		SimLib::initCompressionBox(&sim, sample_filename, false, false, true, init_wall_speed, 0, stop_filling_factor);
	}
	else if(argc == 8)
	{
		SimLib::initCompressionBox(&sim, sample_filename, false, false, true, wall_speed, 0, stop_filling_factor);
		sim.startSimulation(sim.box_height / wall_speed, timestep, 0, 0, NULL, NULL);
	}
	else if(argc == 7)// continue simulation
	{
		sim.loadFromFile(sample_filename);

		wall_speed = sim.walls[sim.wall_upper_id].velocity.norm();
		sim.startSimulation(sim.box_height / wall_speed, timestep, 0, 0, NULL, NULL);
	}

	// prepare log file
	FILE *log_file;

	if(argc == 8 || argc == 11)
	{
		log_file = fopen(log_filename, "w+");
		fprintf(log_file, "# number of particles: %i\n", sim.number_of_particles);
		fprintf(log_file, "# wall speed: %lf cm/s\n", wall_speed);
		fprintf(log_file, "# filling factor: %lf\n", sim.getFillingFactor());
		fprintf(log_file, "# box base size: %g m^2\n", 1e-4 * sim.box_base);
		fprintf(log_file, "# box height: %g m\n", 0.01 * sim.box_height);
		fprintf(log_file, "# time      force top (in F_c)      top pressure      filling factor      bottom pressure \n");
		fclose(log_file);
	}

	// prepare force averaging
	double *top_forces = new double[log_interval];
	memset(top_forces, 0, sizeof(double) * log_interval);
	
	double *bottom_forces = new double[log_interval];
	memset(bottom_forces, 0, sizeof(double) * log_interval);

	// run
	int log_counter = 0;
	int save_counter = 0;
	int save_interval = 1e4;
	double percentage = 0;
	double pressure_top, pressure_bottom;
	double filling_factor;
	double force_top, force_bottom;

	double critical_filling_factor = 0.48; // if this filling factor is exceeded, the logging interval is reduced
	double initial_filling_factor = sim.getFillingFactor();

	printf("%.2lf - 0.00%% done", initial_filling_factor);

	// low velocity (not logging)
	if(argc == 11)
	{
		double a = sqrt( (wall_speed*wall_speed - init_wall_speed*init_wall_speed) / (2.0 * init_penetration_depth * particle_radius) );
		double t = ( sqrt(init_wall_speed*init_wall_speed + 2.0 * a * init_penetration_depth * particle_radius) - init_wall_speed ) / a;

		sim.startSimulation(t, init_timestep, 0, 0, NULL, NULL);

		while(!sim.stop_simulation)
		{
			sim.update();
			sim.walls[sim.wall_upper_id].velocity.y -= a * init_timestep;
		}

		sim.walls[sim.wall_upper_id].velocity = double3(0, -wall_speed, 0);
		sim.startSimulation(sim.box_height / wall_speed, timestep, 0, 0, NULL, NULL);
	}

	//unsigned int timer;
	//cutCreateTimer(&timer);
	//cutStartTimer(timer);

	while(!sim.stop_simulation)
	{
		sim.update();

		top_forces[log_counter] = sim.walls[sim.wall_upper_id].getWallForce();
		bottom_forces[log_counter] = sim.walls[sim.wall_lower_id].getWallForce();
		++log_counter;
		++save_counter;

		if(save_counter == save_interval-1)
		{
			save_counter = 0;

			// save result
			char filename[200];
			sprintf(filename, "compaction_%i.dat", sim.number_of_particles);
			sim.saveToFile(filename);
		}

		if(log_counter == log_interval-1)
		{
			log_counter = 0;

			force_top = 0;
			force_bottom = 0;

			for(int i = 0; i < log_interval; ++i)
			{
				force_top += top_forces[i];
				force_bottom += bottom_forces[i];
			}

			force_top /= (double)log_interval;
			force_bottom /= (double)log_interval;

			pressure_top = 0.1 * force_top/sim.box_base; // /10 to convert from CGS to SI
			pressure_bottom = 0.1 * force_bottom/sim.box_base; // /10 to convert from CGS to SI
			filling_factor = sim.getFillingFactor();
	
			log_file = fopen(log_filename, "a");

			if(log_file)
			{
#ifdef LOG_MAX_SPEED
				double max_speed = SimLib::getMaxSpeed(&sim);
				fprintf(log_file, "%g %g %g %g %g\n", sim.current_time,  force_top/F_c, pressure_top, filling_factor, max_speed);
#else
				fprintf(log_file, "%g %g %g %g %g\n", sim.current_time,  force_top/F_c, pressure_top, filling_factor, pressure_bottom);
#endif
				fclose(log_file);
			}

			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			percentage = 100.0 * sim.current_time/sim.end_time;
			printf("%.2lf - %3.2lf%% done", filling_factor, percentage);

			// reduce log interval
			if(filling_factor > critical_filling_factor && filling_factor < 0.55)
			{
				critical_filling_factor += 0.03;
				log_interval = (int) ( 0.7 * (double)log_interval );
			}
		}
	}

	//cutStopTimer(timer);
	//printf("%f ms\n", cutGetTimerValue(timer));
	//cutDeleteTimer(timer);

	delete [] top_forces;
	delete [] bottom_forces;

	// save result
	//char filename[200];
	//sprintf(filename, "compaction_%i.dat", sim.number_of_particles);
	//sim.saveToFile(filename);
    return 0;
}
