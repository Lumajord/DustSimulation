#pragma once

#include "Constants.h"

#include <vector>
#include <list>

struct GridEntry 
{
	int particle;
	GridEntry *next;
};

class SpatialGrid
{
public:
	SpatialGrid(void);
	~SpatialGrid(void);

	void init(int cells, double cell_width, int max_particles);

	void resetGrid();

	void getCellCenter(vec3 *center_pos, double pos_x, double pos_y, double pos_z);

	int getCellID(double pos_x, double pos_y, double pos_z);

	int getXID(double pos_x);

	int getYID(double pos_y);

	int getZID(double pos_z);

	void addParticles(double *pos, int number_of_particles);

	void addParticle(vec3 &pos, int id);

	int getParticleAt(vec3 &pos, double *particle_positions);

	bool canAddParticleAt(vec3 &pos, double *particle_positions);

	int checkForContactAt(vec3 &pos, double *particle_positions);

	int getNumberOfParticlesInCell(int cell_id);

	int getNeighbourCount(int particle_id, int search_dist);

	void getNeighbours(std::list<int> *neighbours, const double *particle_positions, int particle_id, double max_dist);

	double getVolumeFillingFactor(double x_pos, double y_pos, double z_pos, double search_dist);

	int getNumberOfParticlesInRect(double lower_x, double lower_y, double lower_z, double upper_x, double upper_y, double upper_z);

	//void printParticlesPerCell(const char *filename);

	int x_cells;
	int y_cells;
	int z_cells;

	double x_width;
	double y_width;
	double z_width;

	double x_shift, x_max;
	double y_shift, y_max;
	double z_shift, z_max;

	int max_particles;

	// entries for the particles_in_cell list (to prevent unnecessary allocation/deallocation when adding partucles to the grid)
	GridEntry *entries;

	// list storing all particles within a certain cell
	GridEntry **particles_in_cell;

	// array storing the id of the cell of each particle
	std::vector<int> cell_of_particle;
};
