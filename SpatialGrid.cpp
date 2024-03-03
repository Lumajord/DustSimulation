#include "SpatialGrid.h"

#include <string.h>
#include <cstdio>

extern double particle_radius;
extern double delta_0;

SpatialGrid::SpatialGrid(void)
{
	particles_in_cell = NULL;
	entries = NULL;
}

SpatialGrid::~SpatialGrid(void)
{
	if(particles_in_cell)
	{
		delete [] particles_in_cell;
		particles_in_cell = NULL;
	}
	
	if(entries)
	{
		delete [] entries;
		entries = NULL;
	}
}

void SpatialGrid::init(int cells, double cell_width, int max_particles)
{
	if(entries != NULL)
	{
		delete [] particles_in_cell;
		delete [] entries;
	}

	x_cells = y_cells = z_cells = cells;
	x_width = y_width = z_width = cell_width;
	this->max_particles = max_particles;

	particles_in_cell = new GridEntry*[x_cells * y_cells * z_cells];
	memset(particles_in_cell, 0, sizeof(GridEntry*) * (x_cells * y_cells * z_cells));

	entries = new GridEntry[max_particles];
	for(int p = 0; p < max_particles; ++p)
	{
		entries[p].next = NULL;
		entries[p].particle = p;
	}

	cell_of_particle.clear();
	cell_of_particle.resize(max_particles, 0);

	x_max = x_width * (double)(x_cells - 4);
	x_shift = 0.5 * x_max;

	y_max = y_width * (double)(y_cells - 4);
	y_shift = 0.5 * y_max;

	z_max =  z_width * (double)(z_cells - 4);
	z_shift = 0.5 * z_max;
}

void SpatialGrid::resetGrid()
{

    if(particles_in_cell == NULL)
    {
        printf("Error particles_in_cell not initialized in resetGrid()!\n");
        return;

    }

	for(int p = 0; p < max_particles; ++p)
		particles_in_cell[cell_of_particle[p]] = NULL;
}

void SpatialGrid::getCellCenter(vec3 *center_pos, double pos_x, double pos_y, double pos_z)
{
	pos_x += x_shift; 
	pos_y += y_shift;
	pos_z += z_shift;

	(*center_pos)[0] = (0.5 + floor(pos_x/x_width) ) * x_width;
	(*center_pos)[1] = (0.5 + floor(pos_y/y_width) ) * y_width;
	(*center_pos)[2] = (0.5 + floor(pos_z/z_width) ) * z_width;

	(*center_pos)[0] -= x_shift;
	(*center_pos)[1] -= y_shift;
	(*center_pos)[2] -= z_shift;
}

int SpatialGrid::getCellID(double pos_x, double pos_y, double pos_z)
{
	// determine in which cell pos is located
	int x, y, z;
	double grid_pos_x = pos_x + x_shift;
	double grid_pos_y = pos_y + y_shift;
	double grid_pos_z = pos_z + z_shift;

	if(grid_pos_x < 0)
		x = 1;
	else if(grid_pos_x > x_max)
		x = x_cells-2;
	else
		x = 1 + grid_pos_x / x_width;

	if(grid_pos_y < 0)
		y = 1;
	else if(grid_pos_y > y_max)
		y = y_cells-2;
	else
		y = 1 + grid_pos_y /y_width;

	if(grid_pos_z < 0)
		z = 1;
	else if(grid_pos_z > z_max)
		z = z_cells-2;
	else
		z = 1 + grid_pos_z / z_width;

	return (x + y * x_cells + z * x_cells * y_cells);
}

int SpatialGrid::getXID(double pos_x)
{
	double grid_pos_x = pos_x + x_shift;

	if(grid_pos_x < 0)
		return 1;
	else if(grid_pos_x > x_max)
		return x_cells-2;
	else
		return 1 + grid_pos_x / x_width;
}

int SpatialGrid::getYID(double pos_y)
{
	double grid_pos_y = pos_y + y_shift;

	if(grid_pos_y < 0)
		return 1;
	else if(grid_pos_y > y_max)
		return y_cells-2;
	else
		return 1 + grid_pos_y / y_width;
}

int SpatialGrid::getZID(double pos_z)
{
	double grid_pos_z = pos_z + z_shift;

	if(grid_pos_z < 0)
		return 1;
	else if(grid_pos_z > z_max)
		return z_cells-2;
	else
		return 1 + grid_pos_z / z_width;
}

void SpatialGrid::addParticles(double *pos, int number_of_particles)
{
	for(int p = 0; p < number_of_particles; ++p)
    {
		int cell_id = getCellID(pos[X_COORD(p)], pos[Y_COORD(p)], pos[Z_COORD(p)]);

		cell_of_particle[p] = cell_id;
        entries[p].next = NULL;
		GridEntry *next_entry = particles_in_cell[cell_id];

		if(next_entry)
		{  
			// find last particle in list
			while(next_entry->next)
				next_entry = next_entry->next;

			next_entry->next = &entries[p];
		}
		else
			particles_in_cell[cell_id] = &entries[p];
	}
}

void SpatialGrid::addParticle(vec3 &pos, int id)
{
	int cell_id = getCellID(pos[0], pos[1], pos[2]);

	// add particle to grid
	cell_of_particle[id] = cell_id;
	entries[id].next = NULL;
	GridEntry *grid_entry = particles_in_cell[cell_id];

	if(grid_entry)
	{  
		// find last particle in list
		while(grid_entry->next)
			grid_entry = grid_entry->next;

		grid_entry->next = &entries[id];
	}
	else
		particles_in_cell[cell_id] = &entries[id];
}

int SpatialGrid::getParticleAt(vec3 &pos, double *particle_positions)
{
	// determine in which cell pos is located
	int center_cell_id = getCellID(pos[0], pos[1], pos[2]);
		
	double particle_dist_squared = 4.0 * particle_radius * particle_radius;

	// check particles in neighbouring cells
	for(int z = -1; z <= 1; ++z)
	{
		for(int y = -1; y <= 1; ++y)
		{	
			int cell_id = center_cell_id - 1 - y * x_cells - z * x_cells * y_cells;

			for(int x = -1; x <= 1; ++x)
			{
				// determine particles in this cell
				GridEntry *grid_entry = particles_in_cell[cell_id];

				while(grid_entry)
				{
					// determine dist
					vec3 delta_pos;
					delta_pos[0] = particle_positions[X_COORD(grid_entry->particle)] - pos[0];
					delta_pos[1] = particle_positions[Y_COORD(grid_entry->particle)] - pos[1];
					delta_pos[2] = particle_positions[Z_COORD(grid_entry->particle)] - pos[2];

					double dist_squared = norm_squared(delta_pos);

					if(dist_squared < particle_dist_squared)
						return grid_entry->particle;
					else
						grid_entry = grid_entry->next;
				}

				++cell_id;
			}
		}
	}

	return -1;
}

bool SpatialGrid::canAddParticleAt(vec3 &pos, double *particle_positions)
{
	// determine in which cell pos is located
	int center_cell_id = getCellID(pos[0], pos[1], pos[2]);

	double particle_dist_squared = (2.0 * particle_radius - 1.1 * delta_0)*(2.0 * particle_radius - 1.1 * delta_0);
	//double particle_dist_squared = 4.0 * particle_radius * particle_radius;

	// check particles in neighbouring cells
	for(int z = -1; z <= 1; ++z)
	{
		for(int y = -1; y <= 1; ++y)
		{	
			int cell_id = center_cell_id - 1 - y * x_cells - z * x_cells * y_cells;

			for(int x = -1; x <= 1; ++x)
			{
				// determine particles in this cell
				GridEntry *grid_entry = particles_in_cell[cell_id];

				while(grid_entry)
				{
					// determine dist
					vec3 delta_pos;
					delta_pos[0] = particle_positions[X_COORD(grid_entry->particle)] - pos[0];
					delta_pos[1] = particle_positions[Y_COORD(grid_entry->particle)] - pos[1];
					delta_pos[2] = particle_positions[Z_COORD(grid_entry->particle)] - pos[2];

					double dist_squared = norm_squared(delta_pos);

					if(dist_squared < particle_dist_squared)
						return false;
					else
						grid_entry = grid_entry->next;
				}

				++cell_id;
			}
		}
	}

	return true;
}

int SpatialGrid::checkForContactAt(vec3 &pos, double *particle_positions)
{
	// determine in which cell pos is located
	int center_cell_id = getCellID(pos[0], pos[1], pos[2]);
	int cell_id;

	// check particles in neighbouring cells
	GridEntry *grid_entry;
	double particle_dist_squared = 4.0 * particle_radius*particle_radius; //(2.0 * particle_radius - delta_0)*(2.0 * particle_radius - delta_0);

	for(int z = -1; z <= 1; ++z)
	{
		for(int y = -1; y <= 1; ++y)
		{	
			cell_id = center_cell_id - 1 - y * x_cells - z * x_cells * y_cells;

			for(int x = -1; x <= 1; ++x)
			{
				// determine particles in this cell
				grid_entry = particles_in_cell[cell_id];

				while(grid_entry)
				{
					// determine dist
					vec3 delta_pos;
					delta_pos[0] = particle_positions[X_COORD(grid_entry->particle)] - pos[0];
					delta_pos[1] = particle_positions[Y_COORD(grid_entry->particle)] - pos[1];
					delta_pos[2] = particle_positions[Z_COORD(grid_entry->particle)] - pos[2];

					double dist_squared = norm_squared(delta_pos);

 					if(dist_squared < particle_dist_squared)
						return grid_entry->particle;
					else
						grid_entry = grid_entry->next;
				}

				++cell_id;
			}
		}
	}

	return -1;
}

int SpatialGrid::getNumberOfParticlesInCell(int cell_id)
{
	int result = 0;
	GridEntry *grid_entry = particles_in_cell[cell_id];

	while(grid_entry)
	{
		++result;
		grid_entry = grid_entry->next;
	}

	return result;
}

int SpatialGrid::getNumberOfParticlesInRect(double lower_x, double lower_y, double lower_z, double upper_x, double upper_y, double upper_z)
{
	int count = 0;
	int cell_id;

	int lower_x_id = getXID(lower_x);
	int lower_y_id = getYID(lower_y);
	int lower_z_id = getZID(lower_z);
	int upper_x_id = getXID(upper_x);
	int upper_y_id = getYID(upper_y);
	int upper_z_id = getZID(upper_z);

	for(int z = lower_z_id; z <= upper_z_id; ++z)
	{
		for(int y = lower_y_id; y <= upper_y_id; ++y)
		{
			cell_id = lower_x_id + y * x_cells + z * x_cells * y_cells;

			for(int x = lower_x_id; x <= upper_x_id; ++x)
			{
				count += getNumberOfParticlesInCell(cell_id);
				++cell_id;
			}	
		}
	}

	return count;
}

int SpatialGrid::getNeighbourCount(int particle_id, int search_dist)
{
	int result = 0;
	int cell_id;

	for(int z = -search_dist; z <= search_dist; ++z)
	{
		for(int y = -search_dist; y <= search_dist; ++y)
		{
			cell_id = cell_of_particle[particle_id] - search_dist - y * x_cells - z * x_cells * y_cells;

			for(int x = -search_dist; x <= search_dist; ++x)
			{
				// get particles in this cell
				if(cell_id >= 0 && cell_id < x_cells*y_cells*z_cells)	
					result += getNumberOfParticlesInCell(cell_id);

				++cell_id;
			}
		}
	}

	return result;
}

void SpatialGrid::getNeighbours(std::list<int> *neighbours, const double *particle_positions, int particle_id, double max_dist)
{
	neighbours->clear();

	GridEntry *next_entry;
	int cell_id;
	int search_dist = 3;

	for(int z = -search_dist; z <= search_dist; ++z)
	{
		for(int y = -search_dist; y <= search_dist; ++y)
		{
			cell_id = cell_of_particle[particle_id] - search_dist - y * x_cells - z * x_cells * y_cells;

			for(int x = -search_dist; x <= search_dist; ++x)
			{
				// get particles in this cell
				if(cell_id >= 0 && cell_id < x_cells*y_cells*z_cells)	
				{
					next_entry = particles_in_cell[cell_id];

					while(next_entry)
					{
						if(next_entry->particle != particle_id)
						{
							double dist_squared = (particle_positions[X_COORD(particle_id)]-particle_positions[X_COORD(next_entry->particle)]) * (particle_positions[X_COORD(particle_id)]-particle_positions[X_COORD(next_entry->particle)])
								+ (particle_positions[Y_COORD(particle_id)]-particle_positions[Y_COORD(next_entry->particle)]) * (particle_positions[Y_COORD(particle_id)]-particle_positions[Y_COORD(next_entry->particle)])
								+ (particle_positions[Z_COORD(particle_id)]-particle_positions[Z_COORD(next_entry->particle)]) * (particle_positions[Z_COORD(particle_id)]-particle_positions[Z_COORD(next_entry->particle)]);
						
							if(dist_squared < max_dist*max_dist)
								neighbours->push_back(next_entry->particle);
						}

						next_entry = next_entry->next;
					}
				}

				++cell_id;
			}
		}
	}
}

double SpatialGrid::getVolumeFillingFactor(double x_pos, double y_pos, double z_pos, double search_dist)
{
	// to get exact filling factor particles that are only partially in the considered volume need to taken into accounbt as well
	double max_particle_dist = search_dist + particle_radius;

	// determine cell
	int cell_id = getCellID(x_pos, y_pos, z_pos);

	int x_min = cell_id - max_particle_dist / x_width;
	int x_max = cell_id + max_particle_dist / x_width;
	int y_min = cell_id - max_particle_dist / y_width;
	int y_max = cell_id + max_particle_dist / y_width;
	int z_min = cell_id - max_particle_dist / z_width;
	int z_max = cell_id + max_particle_dist / z_width;

	for(int z = z_min; z <= z_max; ++z)
	{
		for(int y = y_min; y <= y_max; ++y)
		{
			cell_id = x_min + y * x_cells + z * x_cells * y_cells;

			for(int x = x_min; x <= x_max; ++x)
			{
				//count += getNumberOfParticlesInCell(cell_id);
				//++cell_id;
			}	
		}
	}

	return 0;
}



/*void SpatialGrid::printParticlesPerCell(const char *filename)
{
	std::vector<int> particle_count(x_cells * y_cells * z_cells, 0);

	for(int p = 0; p < max_particles; ++p)
		++particle_count[ cell_of_particle[p] ];
	
	FILE *file = fopen(filename, "w+");
	int total_count = 0;
	int cell_count = 0;

	if(file)
	{
		for(int cell = 0; cell < x_cells * y_cells * z_cells; ++cell)
		{
			if(particle_count[cell]>0)
			{
				fprintf(file, "%i: %i\n", cell, particle_count[cell]);
				total_count += particle_count[cell];
				++cell_count;
			}
		}

		fprintf(file, "Total count: %i  Max Particles %i  Cell Count %i\n", total_count, max_particles, cell_count);

		fclose(file);
	}
}*/
