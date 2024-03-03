#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>
#include <math.h>

double bin_size;
double min_bin_value;
double max_bin_value;

int main(int argc, char **argv)
{
	int result_file_index = 1;

	int i;
	scanf("%i", &i);

	if(argc >= 6)
	{
		bin_size = atof(argv[2]);
		min_bin_value = atof(argv[3]);
		max_bin_value = atof(argv[4]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n-result_filename -bin_size -min_pull_distance (in µm) -max_pull_distance (in µm) -filename_1 -filename_2 ... -fileame_n\n");
		return EXIT_SUCCESS;
	}

	if(min_bin_value > max_bin_value)
	{
		printf("Error: min_bin_value > max_bin_value\n");
		return EXIT_SUCCESS;
	}

	int bins = 1 + (int) ( (max_bin_value - min_bin_value) / bin_size );

	////////////////////////////////////////////////////////////////////////////////////////////////
	// create arrays
	////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector< std::list<double> > pressures_top(bins);
	std::vector<double> mean_pressure_top(bins, 0);
	std::vector<double> deviation_pressure_top(bins, 0);

	std::vector< std::list<double> > pressures_bottom(bins);
	std::vector<double> mean_pressure_bottom(bins, 0);
	std::vector<double> deviation_pressure_bottom(bins, 0);

	double max_pressure_top = 0;
	double max_pressure_bottom = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// read values and sort in bins
	////////////////////////////////////////////////////////////////////////////////////////////////

	for(int i = 5; i < argc; ++i)
	{
		FILE *file_curve = fopen(argv[i], "r");

		if(!file_curve)
		{
			printf("ERROR: Could not read %s!\n", argv[i]);
		}
		else
		{
			// get to second line (skip description)
			for(unsigned int l = 0; l <= 1; ++l) 
			{
				char buffer[512];

				if(fgets(buffer, 512, file_curve) == NULL)
				{
					printf("Error while reading %s\n", argv[i]); 
					return EXIT_SUCCESS;
				}
			}

			// read file
			double box_height, pull_distance, top_pressure, bottom_pressure;

			while(EOF != fscanf(file_curve, "%lf %lf %lf %lf",  &box_height, &pull_distance, &top_pressure, &bottom_pressure))
			{
				// determine bin
				int bin = (int) ( (pull_distance - min_bin_value) / bin_size );

				if(bin >= 0 && bin < bins)
				{
					pressures_top[bin].push_back(top_pressure);
					pressures_bottom[bin].push_back(bottom_pressure);
				}
			}

			fclose(file_curve);
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate mean values and deviation
	////////////////////////////////////////////////////////////////////////////////////////////////

	for(int bin = 0; bin < bins; ++bin)
	{
		for(std::list<double>::iterator p = pressures_top[bin].begin(); p != pressures_top[bin].end(); ++p)
		{
			mean_pressure_top[bin] += *p;
			deviation_pressure_top[bin] += *p * *p;
		}

		for(std::list<double>::iterator p = pressures_bottom[bin].begin(); p != pressures_bottom[bin].end(); ++p)
		{
			mean_pressure_bottom[bin] += *p;
			deviation_pressure_bottom[bin] += *p * *p;
		}

		// calculate mean values
		if(pressures_top[bin].size() > 0)
		{
			mean_pressure_top[bin] /= (double)pressures_top[bin].size();
			deviation_pressure_top[bin] /= (double)pressures_top[bin].size();
		}

		if(pressures_bottom[bin].size() > 0)
		{
			mean_pressure_bottom[bin] /= (double)pressures_bottom[bin].size();
			deviation_pressure_bottom[bin] /= (double)pressures_bottom[bin].size();
		}

		if(mean_pressure_top[bin] < max_pressure_top)
			max_pressure_top = mean_pressure_top[bin];

		if(mean_pressure_bottom[bin] > max_pressure_bottom)
			max_pressure_bottom = mean_pressure_bottom[bin];
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// print results
	////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *result_file = fopen(argv[result_file_index], "w+");

	if(!result_file)
	{
		printf("Error: Cannot open %s\n", argv[result_file_index]);
		return 0;
	}

	fprintf(result_file, "# Max top / bottom pressure: %lg / %lg\n", fabs(max_pressure_top), fabs(max_pressure_bottom));
	fprintf(result_file, "#\n");
	fprintf(result_file, "# pull distance (µm)    mean/deviation top pressure (Pa)    mean/deviation bottom pressure (Pa)\n");

	for(int bin = 0; bin < bins; ++bin)
	{
		if(mean_pressure_top[bin] != 0)
			fprintf(result_file, "%lg %lg %lg %lg %lg\n", min_bin_value + (double)bin * bin_size, 
					mean_pressure_top[bin], sqrt(deviation_pressure_top[bin] - mean_pressure_top[bin]*mean_pressure_top[bin]), 
					mean_pressure_bottom[bin], sqrt(deviation_pressure_bottom[bin] - mean_pressure_bottom[bin]*mean_pressure_bottom[bin]) );
	}

	fclose(result_file);

    return EXIT_SUCCESS;
} 
