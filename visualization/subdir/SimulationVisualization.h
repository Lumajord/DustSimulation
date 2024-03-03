 /***************************************************************************
 *   Simulation of particle collisions/agglomeration						*
 *   Copyright (C) 2009 Alexander Seizinger									*
 *	 alexander.seizinger[...]gmx.net										*
 *																			*
 *   This program is free software; you can redistribute it and/or modify	*
 *   it under the terms of the GNU General Public License as published by	*
 *   the Free Software Foundation; either version 2 of the License, or		*
 *   (at your option) any later version.									*
 *																			*
 *   This program is distributed in the hope that it will be useful,		*
 *   but WITHOUT ANYs WARRANTY; without even the implied warranty of		*
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the			*
 *   GNU General Public License for more details.							*
 *																			*
 *   You should have received a copy of the GNU General Public License		*
 *   along with this program; if not, write to the							*
 *   Free Software Foundation, Inc.,										*
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.				*
 ***************************************************************************/

#ifndef SIMULATIONVISUALIZATION_H
#define SIMULATIONVISUALIZATION_H

#include "Simulation.h"

////////////////////////////////////////////////////////////////////////////////////////
// constants affecting visualization
////////////////////////////////////////////////////////////////////////////////////////

const double max_speed_offset = 25.0;

// colors for fragments
const int FRAGMENT_COLORS = 7;
const double FRAGMENT_COLOR_R[7] = {1.0, 0.5, 1.0, 0.5, 1.0, 1.0, 0.5};
const double FRAGMENT_COLOR_G[7] = {1.0, 1.0, 0.5, 0.5, 1.0, 0.5, 1.0};
const double FRAGMENT_COLOR_B[7] = {1.0, 0.5, 0.5, 1.0, 0.5, 1.0, 1.0};

// limits the number of frames per second
const int MAX_FPS = 30;

struct AgglomerateInfo
{
	double gyration_radius;
	double outer_radius;
	double filling_factor_sphere;
	double filling_factor_box;
	double box_height;
	double box_base;
	int fragments;
	int contact_histogram[13];
};

struct DisplayInfo
{
	int particles;
	int broken_contacts;
	int created_contacts;
	double time;
	double filling_factor;
	double wall_speed;
	double collision_speed;
	double pressure;
	double force;
	double coordination_number;
};

class SimulationVis : public Simulation
{
public:
	SimulationVis(void);
	~SimulationVis(void);

	////////////////////////////////////////////////////////////////////////////////////////
	// extended functions
	////////////////////////////////////////////////////////////////////////////////////////

	void update();

	void predictor(double timestep);

    ErrorCode initSimData(SimType type, double collision_speed = 0.0, double stop_filling_factor = 0.0, double stop_dissipation_factor = 1.0, double initial_kinetic_energy = 0.0, double shockwave_penetration_depth = 0.0);

	ErrorCode resizeArrays(int new_number_of_particles, int new_number_of_walls);

	ErrorCode addParticlesFromSim(Simulation &sim);

	ErrorCode addParticles(int particles);

	////////////////////////////////////////////////////////////////////////////////////////
	// visualization related routines
	////////////////////////////////////////////////////////////////////////////////////////

	// returns current (averaged) pressure on top wall
	double getAvgWallPressure();
	double getAvgWallForce();

	void initPressureAveraging(int averaging_steps);

	// returns the depth color modifier for a certain particle and updates max_neighbours_new for later usage
	float getDepthColor(int particle_id);

	// returns the depth color modifier of a particle with specified number of neigbours 
	float calcDepthColor(int number_of_neighbours);

	// determine various values needed to apply proper color scales
	void prepareColorScaling();

	// copies the radius of particles to specified vector
	void storeParticlePos(double *dest);

	// copies the pos and normal of all visible walls to specified vector
	void storeWalls(double *new_pos, double *normals, float *alpha_values);

	// sets particle colors according to various criteria
	void storeDissipatedContactEnergy(float *dest, bool depth);
	void storeColorParticleContacts(float *dest, bool depth);
	void storeColorWallContacts(float *dest, bool depth);
	void storeWallAngle(float *dest, bool depth);
	void storeParticleSpeed(float *dest, bool depth);
	void storeDensity(float *dest, bool depth);
	void storeDepth(float *dest, bool depth);
	void storeFragments(float *dest, bool depth);
	void storeDislocation(float *dest, bool depth);
    void store50(float *dest, bool depth);

	////////////////////////////////////////////////////////////////////////////////////////
	// color functions (static to be accessible by gl widget)
	////////////////////////////////////////////////////////////////////////////////////////

	static void getColorDensity(float *color, float density);
	static void getColorCoordinationNumber(float *color, float k);
	static void getColorDislocation(float *color, float k);
    static void getColor50(float *color, int i, int N);

	////////////////////////////////////////////////////////////////////////////////////////
	// visualization related variables
	////////////////////////////////////////////////////////////////////////////////////////

	// maximum number of neighbours/ dissipated energy per particle in the last update step (used for depth visualization)
	int max_neighbours;
	int max_neighbours_new;

	double max_diss_energy;
	double max_speed;
	vec3 avg_speed;

	// averaging of pressure/force
	double *forces_storage;
	int force_avg_size;
	int force_avg_counter;

	// ids indicate which fragment a certain particle belongs to
	int *fragment_ids;
	int fragment_id_offset;

	// array storing initial positions of the particles (used to visualize traveled distance)
	double *initial_positions;
};

#endif
