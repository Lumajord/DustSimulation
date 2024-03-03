#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <vector>
#include <list>
#include <math.h>


#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
  #undef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
#endif

#ifdef TRACK_DISSIPATED_ENERGY
  #undef TRACK_DISSIPATED_ENERGY
#endif

double time_bin_size;
double min_time;
double max_time;


/*double calcDeviation(std::vector<double> *mean_pressures, double p_m, double k)
{
	double result = 0;
	double f;
	double phi;

	for(int i = 0; i < mean_pressures->size(); ++i)
	{ 
		phi = min_filling_factor + (double)i * (max_filling_factor-min_filling_factor);

		f = max_filling_factor - (max_filling_factor-min_filling_factor) / ( 1.0 + exp( (log(phi) - log(p_m)) / k ) );
		result += (f - (*mean_pressures)[i]) * (f - (*mean_pressures)[i]);
	}

	return result;
}*/

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

	char result_filename[200];

	if(argc >= 6)
	{
		sprintf(result_filename, "%s", argv[1]);
		time_bin_size = atof(argv[2]);
		min_time = atof(argv[3]);
		max_time = atof(argv[4]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n-result_filename -time_bin_size -min_time -max_time -filename_1 -filename_2 ... -fileame_n\n");
		return 0;
	}

	if(min_time > max_time)
	{
		printf("Error: min_time > max_time\n");
		return 0;
	}

	// create array
	int bins = 1 +  (int) ( (max_time - min_time) / time_bin_size );
	
	std::vector< std::list<double> > top_pressure(bins);
	std::vector<double> top_mean_pressure(bins, 0);

	std::vector< std::list<double> > bottom_pressure(bins);
	std::vector<double> bottom_mean_pressure(bins, 0);

	std::vector< std::list<double> > kinetic_energy(bins);
	std::vector<double> mean_kinetic_energy(bins, 0);

	double time, t_pressure, b_pressure, kin_energy;

	for(int i = 5; i < argc; ++i)
	{
		FILE *file_curve = fopen(argv[i], "r");

		if(file_curve)
		{
			// get to 7th line (skip description header)
			for(unsigned int l = 0; l <= 6; ++l) 
			{
				char buffer[512];

				if(fgets(buffer, 512, file_curve) == NULL)
				{
					printf("Error while reading %s\n", argv[i]); 
					return -1;
				}
			}

			// read file
			while(EOF != fscanf(file_curve, "%lf %lf %lf %lf",  &time, &t_pressure, &b_pressure, &kin_energy))
			{
				// determine bin
				int bin = (int) ( (time - min_time) / time_bin_size );

				if(bin >= 0 && bin < bins)
				{
					top_pressure[bin].push_back(t_pressure);
					bottom_pressure[bin].push_back(b_pressure);
					kinetic_energy[bin].push_back(kin_energy);
				}
			}
		}
	}

	for(int bin = 0; bin < bins; ++bin)
	{
		if(top_pressure[bin].size() > 0)
		{
			for(std::list<double>::iterator p = top_pressure[bin].begin(); p != top_pressure[bin].end(); ++p)
				top_mean_pressure[bin] += *p;

			top_mean_pressure[bin] /= (double)top_pressure[bin].size();


			for(std::list<double>::iterator p = bottom_pressure[bin].begin(); p != bottom_pressure[bin].end(); ++p)
				bottom_mean_pressure[bin] += *p;

			bottom_mean_pressure[bin] /= (double)bottom_pressure[bin].size();

			for(std::list<double>::iterator p = kinetic_energy[bin].begin(); p != kinetic_energy[bin].end(); ++p)
				mean_kinetic_energy[bin] += *p;

			mean_kinetic_energy[bin] /= (double)kinetic_energy[bin].size();
		}
	}

	// printf results
	FILE *result_file = fopen(result_filename, "w+");

	if(result_file)
	{
		fprintf(result_file, "# time     top pressure     bottom pressure    kinetic energy\n");

		for(int bin = 0; bin < bins; ++bin)
			fprintf(result_file, "%g %g %g %g\n", min_time + (0.5 + (double)bin) * time_bin_size, top_mean_pressure[bin], bottom_mean_pressure[bin], mean_kinetic_energy[bin]);

		fclose(result_file);
	}
	else
		printf("Error: Cannot open %s\n", result_filename);

    return 0;
} 
