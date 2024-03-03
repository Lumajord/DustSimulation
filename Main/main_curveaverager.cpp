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


#define NEW_FORMAT

//#define CROSS_SECTION	// should be enabled when reading files from unidirectional compression
//#define MAX_TOP_PRESSURE // enabled to read max pressure

double bin_size;
double min_bin_value;
double max_bin_value;

/*double calcDeviation(std::vector<double> *mean_pressures, double p_m, double k)
{
	double result = 0;
	double f;
	double phi;

	for(int i = 0; i < mean_pressures->size(); ++i)
	{ 
		phi = min_bin_value + (double)i * (max_bin_value-min_bin_value);

		f = max_bin_value - (max_bin_value-min_bin_value) / ( 1.0 + exp( (log(phi) - log(p_m)) / k ) );
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

	if(argc >= 7)
	{
		sprintf(result_filename, "%s", argv[1]);
		bin_size = atof(argv[3]);
		min_bin_value = atof(argv[4]);
		max_bin_value = atof(argv[5]);
	}
	else
	{
		printf("Incorrect arguments! Use:\n-result_filename -bin_size -min_bin_value -max_bin_value -filename_1 -filename_2 ... -fileame_n\n");
		return 0;
	}

	if(min_bin_value > max_bin_value)
	{
		printf("Error: min_bin_value > max_bin_value\n");
		return 0;
	}

	
	int bins = 1 + (int) ( (max_bin_value - min_bin_value) / bin_size );

	////////////////////////////////////////////////////////////////////////////////////////////////
	// create arrays
	////////////////////////////////////////////////////////////////////////////////////////////////

	double time, box_filling_factor, box_pressure_top, box_pressure_bottom;

	std::vector< std::list<double> > box_pressures_top(bins);
	std::vector<double> mean_box_pressure_top(bins, 0);
	std::vector<double> deviation_box_pressure_top(bins, 0);

	std::vector< std::list<double> > box_pressures_bottom(bins);
	std::vector<double> mean_box_pressure_bottom(bins, 0);
	std::vector<double> deviation_box_pressure_bottom(bins, 0);

	std::vector< std::list<double> > box_filling_factors(bins);
	std::vector<double> mean_box_filling_factor(bins, 0);
	std::vector<double> deviation_box_filling_factor(bins, 0);

	std::vector< std::list<double> > cs_filling_factors(bins);
	std::vector<double> mean_cs_filling_factor(bins, 0);
	std::vector<double> deviation_cs_filling_factor(bins, 0);
	
#ifdef CROSS_SECTION
	double cross_section, cs_filling_factor, cs_pressure_top, cs_pressure_bottom;

	std::vector< std::list<double> > cs_pressures_top(bins);
	std::vector<double> mean_cs_pressure_top(bins, 0);
	std::vector<double> deviation_cs_pressure_top(bins, 0);

	std::vector< std::list<double> > cs_pressures_bottom(bins);
	std::vector<double> mean_cs_pressure_bottom(bins, 0);
	std::vector<double> deviation_cs_pressure_bottom(bins, 0);

	std::vector< std::list<double> > cross_sections(bins);
	std::vector<double> mean_cross_section(bins, 0);
	std::vector<double> deviation_cross_section(bins, 0);
#endif

#ifdef MAX_TOP_PRESSURE
	double max_pressure_top;

	std::vector< std::list<double> > max_pressures_top(bins);
	std::vector<double> mean_max_pressure_top(bins, 0);
	std::vector<double> deviation_max_pressure_top(bins, 0);
#endif

	////////////////////////////////////////////////////////////////////////////////////////////////
	// read values and sort in bins
	////////////////////////////////////////////////////////////////////////////////////////////////

	int line;
	int box_bin, cs_bin;
	double min_pressure = pow(10.0, min_bin_value);

	for(int i = 6; i < argc; ++i)
	{
		FILE *file_curve = fopen(argv[i], "r");

		if(!file_curve)
		{
			printf("ERROR: Could not read %s!\n", argv[i]);
		}
		else
		{
			// get to fith line (skip description)
			for(unsigned int l = 0; l <= 7; ++l) 
			{
				char buffer[512];

				if(fgets(buffer, 512, file_curve) == NULL)
				{
					printf("Error while reading %s\n", argv[i]); 
					return -1;
				}
			}

			// read file
			line = 0;

#ifdef NEW_FORMAT

#ifdef CROSS_SECTION
 #ifdef MAX_TOP_PRESSURE
			while(EOF != fscanf(file_curve, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &time, &box_filling_factor, &cs_filling_factor, &cross_section, &box_pressure_top, &cs_pressure_top, &box_pressure_bottom, &cs_pressure_bottom, &max_pressure_top))
 #else
			while(EOF != fscanf(file_curve, "%lf %lf %lf %lf %lf %lf %lf %lf\n",  &time, &box_filling_factor, &cs_filling_factor, &cross_section, &box_pressure_top, &cs_pressure_top, &box_pressure_bottom, &cs_ressure_bottom))	
 #endif
#else // CROSS_SECTION
 #ifdef MAX_TOP_PRESSURE
			while(EOF != fscanf(file_curve, "%lf %lf %lf %lf %lf\n",  &time, &box_filling_factor, &box_pressure_top, &box_pressure_bottom, &max_pressure_top))
 #else
			while(EOF != fscanf(file_curve, "%lf %lf %lf %lf\n",  &time, &box_filling_factor, &box_pressure_top, &box_pressure_bottom))	
 
 #endif
#endif // CROSS_SECTION
			{
				++line;

				if(time > 1.0)
				{
					printf("ERROR: Possible corruption of file %s in line %i!\n", argv[i], line);
					break;
				}

				if(strcmp(argv[2], "PRESSURE" ) == 0)
				{
					// determine bins
					if(box_pressure_top > min_pressure)
					{
						box_bin = (int)( (log10(box_pressure_top) - min_bin_value) / bin_size);//(int) ( (box_pressure_top - min_bin_value) / bin_size );

						if(box_bin >= 0 && box_bin < bins)
							box_filling_factors[box_bin].push_back(box_filling_factor);
					}

#ifdef CROSS_SECTION
					if(cs_pressure_top > min_pressure)
					{
						cs_bin = (int)( (log10(cs_pressure_top) - min_bin_value) / bin_size); //(int) ( (cs_pressure_top - min_bin_value) / bin_size );

						if(cs_bin >= 0 && cs_bin < bins)
							cs_filling_factors[cs_bin].push_back(cs_filling_factor);
					}
#endif
				}
				else if(strcmp(argv[2], "FILLING_FACTOR" ) == 0)
				{
					box_bin = (int) ( (box_filling_factor - min_bin_value) / bin_size );

					if(box_bin >= 0 && box_bin < bins)
					{
						box_pressures_top[box_bin].push_back(box_pressure_top);
						box_pressures_bottom[box_bin].push_back(box_pressure_bottom);

#if defined MAX_TOP_PRESSURE && !defined CROSS_SECTION
						max_pressures_top[box_bin].push_back(max_pressure_top);
#endif
					}

#ifdef CROSS_SECTION
					cs_bin = (int) ( (cs_filling_factor - min_bin_value) / bin_size );

					// sort values in bins
					if(cs_bin >= 0 && cs_bin < bins)
					{
						cs_pressures_top[cs_bin].push_back(cs_pressure_top);
						cs_pressures_bottom[cs_bin].push_back(cs_pressure_bottom);
						cross_sections[cs_bin].push_back(cross_section);
 #ifdef MAX_TOP_PRESSURE
						max_pressures_top[cs_bin].push_back(max_pressure_top);
 #endif
					}
#endif
				}
			}
#else // NEW_FORMAT
			double force;

			while(EOF != fscanf(file_curve, "%lf %lf %lf %lf %lf",  &time, &force, &cs_pressure_top, &cs_filling_factor, &cs_pressure_bottom))
			{
				// determine bin
				int bin = (int) ( (cs_filling_factor - min_bin_value) / bin_size );

				if(bin >= 0 && bin < bins)
				{
					cs_pressures_top[bin].push_back(cs_pressure_top);
					cs_pressures_bottom[bin].push_back(cs_pressure_bottom);
				}
			}
#endif // NEW_FORMAT
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate mean values and deviation
	////////////////////////////////////////////////////////////////////////////////////////////////

	for(int bin = 0; bin < bins; ++bin)
	{
#ifdef CROSS_SECTION
		for(std::list<double>::iterator p = cs_pressures_top[bin].begin(); p != cs_pressures_top[bin].end(); ++p)
		{
			mean_cs_pressure_top[bin] += *p;
			deviation_cs_pressure_top[bin] += *p * *p;
		}

		for(std::list<double>::iterator p = cs_pressures_bottom[bin].begin(); p != cs_pressures_bottom[bin].end(); ++p)
		{
			mean_cs_pressure_bottom[bin] += *p;
			deviation_cs_pressure_bottom[bin] += *p * *p;
		}

		for(std::list<double>::iterator p = cs_filling_factors[bin].begin(); p != cs_filling_factors[bin].end(); ++p)
		{
			mean_cs_filling_factor[bin] += *p;
			deviation_cs_filling_factor[bin] += *p * *p;
		}

		for(std::list<double>::iterator p = cross_sections[bin].begin(); p != cross_sections[bin].end(); ++p)
		{
			mean_cross_section[bin] += *p;
			deviation_cross_section[bin] += *p * *p;
		}

#endif

		for(std::list<double>::iterator p = box_pressures_top[bin].begin(); p != box_pressures_top[bin].end(); ++p)
		{
			mean_box_pressure_top[bin] += *p;
			deviation_box_pressure_top[bin] += *p * *p;
		}

		for(std::list<double>::iterator p = box_pressures_bottom[bin].begin(); p != box_pressures_bottom[bin].end(); ++p)
		{
			mean_box_pressure_bottom[bin] += *p;
			deviation_box_pressure_bottom[bin] += *p * *p;
		}

		for(std::list<double>::iterator p = box_filling_factors[bin].begin(); p != box_filling_factors[bin].end(); ++p)
		{
			mean_box_filling_factor[bin] += *p;
			deviation_box_filling_factor[bin] += *p * *p;
		}


#ifdef MAX_TOP_PRESSURE
		for(std::list<double>::iterator p = max_pressures_top[bin].begin(); p != max_pressures_top[bin].end(); ++p)
		{
			mean_max_pressure_top[bin] += *p;
			deviation_max_pressure_top[bin] += *p * *p;
		}
#endif
	
		// calculate mean values
		if(box_pressures_top[bin].size() > 0)
		{
			mean_box_pressure_top[bin] /= (double)box_pressures_top[bin].size();
			deviation_box_pressure_top[bin] /= (double)box_pressures_top[bin].size();
		}

		if(box_pressures_bottom[bin].size() > 0)
		{
			mean_box_pressure_bottom[bin] /= (double)box_pressures_bottom[bin].size();
			deviation_box_pressure_bottom[bin] /= (double)box_pressures_bottom[bin].size();
		}
		
		if(box_filling_factors[bin].size() > 0)
		{
			mean_box_filling_factor[bin] /= (double)box_filling_factors[bin].size();
			deviation_box_filling_factor[bin] /= (double)box_filling_factors[bin].size();
		}

#ifdef CROSS_SECTION
		if(cs_pressures_top[bin].size() > 0)
		{
			mean_cs_pressure_top[bin] /= (double)cs_pressures_top[bin].size();
			deviation_cs_pressure_top[bin] /= (double)cs_pressures_top[bin].size();
		}

		if(cs_pressures_bottom[bin].size() > 0)
		{
			mean_cs_pressure_bottom[bin] /= (double)cs_pressures_bottom[bin].size();
			deviation_cs_pressure_bottom[bin] /= (double)cs_pressures_bottom[bin].size();
		}

		if(cross_sections[bin].size() > 0)
		{
			mean_cross_section[bin] /= (double)cross_sections[bin].size();
			deviation_cross_section[bin] /= (double)max_pressures_top[bin].size();
		}

		if(cs_filling_factors[bin].size() > 0)
		{
			mean_cs_filling_factor[bin] /= (double)cs_filling_factors[bin].size();
			deviation_cs_filling_factor[bin] /= (double)cs_filling_factors[bin].size();
		}
#endif
		
#ifdef MAX_TOP_PRESSURE
		if(max_pressures_top[bin].size() > 0)
		{
			mean_max_pressure_top[bin] /= (double)max_pressures_top[bin].size();
			deviation_max_pressure_top[bin] /= (double)max_pressures_top[bin].size();
		}
#endif
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// print results
	////////////////////////////////////////////////////////////////////////////////////////////////

	FILE *result_file = fopen(result_filename, "w+");

	if(!result_file)
	{
		printf("Error: Cannot open %s\n", result_filename);
		return 0;
	}

	if(strcmp(argv[2], "PRESSURE" ) == 0)
	{
#ifdef CROSS_SECTION
		fprintf(result_file, "# top pressure     cross-section (avg/deviation)     cross-section filling factor (avg/deviation)     box filling factor (avg/deviation)\n");

		for(int bin = 0; bin < bins; ++bin)
		{
			// skip empty bins
			if(box_filling_factors[bin].size() > 0 || cs_filling_factors[bin].size() > 0)
			{
				double current_pressure = 0.5 * ( pow(10.0, min_bin_value + (double)bin * bin_size) + pow(10.0, min_bin_value + (double)(bin+1) * bin_size) );

				fprintf(result_file, "%g %g %g %g %g %g %g\n", current_pressure, mean_cross_section[bin], deviation_cross_section[bin],
															mean_cs_filling_factor[bin], sqrt(deviation_cs_filling_factor[bin] - mean_cs_filling_factor[bin]*mean_cs_filling_factor[bin]),
															mean_box_filling_factor[bin], sqrt(deviation_box_filling_factor[bin] - mean_box_filling_factor[bin]*mean_box_filling_factor[bin]));

			}
		}
#else
		fprintf(result_file, "# top pressure     filling factor (avg/deviation)\n");

		for(int bin = 0; bin < bins; ++bin)
		{
			// skip empty bins
			if(box_filling_factors[bin].size() > 0)
			{
				double current_pressure = 0.5 * ( pow(10.0, min_bin_value + (double)bin * bin_size) + pow(10.0, min_bin_value + (double)(bin+1) * bin_size) );

				fprintf(result_file, "%g %g %g\n", current_pressure, mean_box_filling_factor[bin], sqrt(deviation_box_filling_factor[bin] - mean_box_filling_factor[bin]*mean_box_filling_factor[bin]));
			}
		}
#endif
	}
	else if(strcmp(argv[2], "FILLING_FACTOR" ) == 0)
	{
#ifdef CROSS_SECTION
 #ifdef MAX_TOP_PRESSURE
		fprintf(result_file, "# filling factor    top cross-section pressure (avg/deviation)     bottom cross-section pressure (avg/deviation)      top box pressure (avg/deviation)      bottom box pressure (avg/deviation)      top max pressure (avg/deviation)\n");
 
		for(int bin = 0; bin < bins; ++bin)
		{
			if(box_pressures_top[bin].size() > 0)
				fprintf(result_file, "%g %g %g %g %g %g %g %g %g %g %g\n", min_bin_value + (double)bin * bin_size, 
														mean_cs_pressure_top[bin], sqrt(deviation_cs_pressure_top[bin] - mean_cs_pressure_top[bin]*mean_cs_pressure_top[bin]),
														mean_cs_pressure_bottom[bin],  sqrt(deviation_cs_pressure_bottom[bin] - mean_cs_pressure_bottom[bin]*mean_cs_pressure_bottom[bin]),
														mean_box_pressure_top[bin], sqrt(deviation_box_pressure_top[bin] - mean_box_pressure_top[bin]*mean_box_pressure_top[bin]),
														mean_box_pressure_bottom[bin],  sqrt(deviation_box_pressure_bottom[bin] - mean_box_pressure_bottom[bin]*mean_box_pressure_bottom[bin]),
														mean_max_pressure_top[bin], sqrt(deviation_max_pressure_top[bin] - mean_max_pressure_top[bin]*mean_max_pressure_top[bin]) );
		}
 #else
		fprintf(result_file, "# filling factor     top cross-section pressure (avg/deviation)     bottom cross-section pressure (avg/deviation)     top box pressure (avg/deviation)      bottom box pressure (avg/deviation)\n");
 
		for(int bin = 0; bin < bins; ++bin)
		{
			if(box_pressures_top[bin].size() > 0)
				fprintf(result_file, "%g %g %g %g %g %g %g %g %g\n", min_bin_value + (double)bin * bin_size, 
														mean_cs_pressure_top[bin], sqrt(deviation_cs_pressure_top[bin] - mean_cs_pressure_top[bin]*mean_cs_pressure_top[bin]),
														mean_cs_pressure_bottom[bin],  sqrt(deviation_cs_pressure_bottom[bin] - mean_cs_pressure_bottom[bin]*mean_cs_pressure_bottom[bin]),
														mean_box_pressure_top[bin], sqrt(deviation_box_pressure_top[bin] - mean_box_pressure_top[bin]*mean_box_pressure_top[bin]),
														mean_box_pressure_bottom[bin],  sqrt(deviation_box_pressure_bottom[bin] - mean_box_pressure_bottom[bin]*mean_box_pressure_bottom[bin]));
		}
#endif
#else // CROSS_SECTION
 #ifdef MAX_TOP_PRESSURE
		fprintf(result_file, "#  filling factor    top pressure (avg/deviation)     bottom pressure (avg/deviation)      top max pressure (avg/deviation)\n");

		for(int bin = 0; bin < bins; ++bin)
		{
			if(box_pressures_top[bin].size() > 0)
				fprintf(result_file, "%g %g %g %g %g %g %g\n", min_bin_value + (double)bin * bin_size, 
														mean_box_pressure_top[bin], sqrt(deviation_box_pressure_top[bin] - mean_box_pressure_top[bin]*mean_box_pressure_top[bin]),
														mean_box_pressure_bottom[bin], sqrt(deviation_box_pressure_bottom[bin] - mean_box_pressure_bottom[bin]*mean_box_pressure_bottom[bin]),
														mean_max_pressure_top[bin], sqrt(deviation_max_pressure_top[bin] - mean_max_pressure_top[bin]*mean_max_pressure_top[bin]) );
		}
 #else
		fprintf(result_file, "#  filling factor    top pressure (avg/deviation)     bottom pressure (avg/deviation)\n");

		for(int bin = 0; bin < bins; ++bin)
		{
			if(box_pressures_top[bin].size() > 0)
				fprintf(result_file, "%g %g %g %g %g\n", min_bin_value + (double)bin * bin_size, 
														mean_box_pressure_top[bin], sqrt(deviation_box_pressure_top[bin] - mean_box_pressure_top[bin]*mean_box_pressure_top[bin]),
														mean_box_pressure_bottom[bin], sqrt(deviation_box_pressure_bottom[bin] - mean_box_pressure_bottom[bin]*mean_box_pressure_bottom[bin]));
		}
 #endif
#endif // CROSS_SECTION
	}

	fclose(result_file);
    return 0;
} 
