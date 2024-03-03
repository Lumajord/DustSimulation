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
	int num_pos_files;
	int result_file_index = 2;

	if(argc == 4)
	{
		num_pos_files = atoi(argv[1]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n-agg1_size -agg2_size -filename\n");
		return EXIT_SUCCESS;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// create arrays
	////////////////////////////////////////////////////////////////////////////////////////////////


	

	////////////////////////////////////////////////////////////////////////////////////////////////
	// read values and sort in bins
	////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *file = fopen(argv[result_file_index], "r");

	if(!file)
		printf("ERROR: Could not read %s!\n", argv[i]);
	else
	{
		// get to second line (skip description)
		char buffer[512];

		if(fgets(buffer, 512, file) == NULL)
		{
			printf("Error while reading %s\n", argv[i]); 
			return EXIT_SUCCESS;
		}

		// read file
		double box_height, pull_distance, top_pressure, bottom_pressure;

		for(unsigned int  fscanf(file, "%lf %lf %lf %lf",  &box_height, &pull_distance, &top_pressure, &bottom_pressure))
		{

		}

		fclose(file_curve);
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
