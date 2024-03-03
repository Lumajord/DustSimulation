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

#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
//#include <stdio.h>
#include <random>

#include "Constants.h"
#include "Contact.h"
#include "Wall.h"
#include "CompressionInterpolation.h"
#include "SpatialGrid.h"




struct ContactListEntry
{
	int id;					// id of the other particle of the contact
	Contact *contact;		// pointer to the contact object
	ContactListEntry *next;	// pointer to the next entry (null if none)
};

struct WallBox
{
	double height;		// in cm
	double x_size;		// in cm
	double y_size;		// in cm
	double base;		// in cm^2
	
	// lower pos of box (to determine filling factor when compressing without side walls)
	double lower_pos_x; 
	double lower_pos_y;
	double lower_pos_z;

	// ids of the lower and upper wall to determine current box height (-1 if none)
	int bottom_wall_id;
	int top_wall_id;
};

struct SimInfo
{
	SimType sim_type;				// determines type of simulation (collision, compaction in box, etc.)
	char material_identifier[12];	// identifier for the material
	double info_storage[8];			// storage for additional information
	// info_storage[0] - wall movement is stopped if this filling factor is reached
	// info_storage[1] - initial kinetic energy
	// info_storage[2] - stores initial collision speed (for visualization)
	// info_storage[3] - distance (in particle radii) the wall will move downwards to trigger shockwave (for sound speed measurement)
	// info_storage[4] - simulation is stopped if this ammount of kinetic energy (with respect to when the wall movement was stopped) has been dissipated
    // info_storage[5] - coordination number
    // info_storage[6] - coordination number of collision partner
    // info_storage[7] - wall stopping height
};

struct DataFileOffsets
{
	int pos_offset;
	int vel_offset;
	int vel_angular_offset;
	int walls_offset;
	int contacts_offset;
};

class Simulation
{
public:
	Simulation();
	~Simulation(void);


	///////////////////////////////////////////////////////////////////////
	// random number generation
	///////////////////////////////////////////////////////////////////////
private:
	std::uniform_real_distribution<double> dist_cos_theta;
	std::uniform_real_distribution<double> dist_zero_twoPi;
	std::uniform_real_distribution<double> dist_zero_one;
	std::uniform_real_distribution<double> dist_zero_one_incl;
public:
	std::mt19937_64 rand_generator;
	double get_random_zero_twoPi();
	double get_random_cos_theta();
	double get_random_zero_one();
	double get_random_zero_one_incl();

	///////////////////////////////////////////////////////////////////////
	// loading/saving
	///////////////////////////////////////////////////////////////////////

	// loads constants from a file and initializes material
	ErrorCode loadMaterial(const char *filename);

	// initializes material with given constants
    void setMaterialConstants(const double particle_radius,
            const double density,
            const double surface_energy,
            const double nu,
            const double young_mod,
            const double crit_rolling_displacement,
            const double osc_damping_factor = 0.0,
            const double T_vis_ = 1.0e-11,
            const double rolling_modifier = 1.0,
            const double sliding_modifier = 1.0,
            const double twisting_modifier = 1.0,
            const double crit_sliding_displacement_modifier = 1.0,
            const double crit_wall_sliding_displacement_modifier = 1.0
            );

	// load simulation data from file
    ErrorCode loadFromFile(const char* filename);

	// saves current pos/vel of particles to file
	ErrorCode saveToFile(const char* filename, bool save_contacts = true);

	// adds particles/contacts of another sim object
	virtual ErrorCode addParticlesFromSim(Simulation &sim);

	// prepare sim to add specified number of particles
	virtual ErrorCode addParticles(int particles);

	// inits arrays, contact list and grid for specified number of particles
	virtual ErrorCode resizeArrays(int new_number_of_particles, int new_number_of_walls);

	// allows to remove some of the particles
	void removeParticles(bool *remove_list, int removed_particles);

	// removes all walls
	void removeWalls();

	// deletes contact list and frees memory of pos, etc. arrays
	void cleanUp();

	// deletes all contacts (but leaves empty contact list)
	void deleteContacts();

    // reset all contacts as if they are newly made
    void resetContacts();

	///////////////////////////////////////////////////////////////////////
	// initialization 
	///////////////////////////////////////////////////////////////////////

	// prepare sim for a fresh run (reset various counters etc.) -> must be done before calling startSimulation(...)
    virtual ErrorCode initSimData(SimType type, double collision_speed = 0.0, double stop_filling_factor = 0.0, double stop_dissipation_factor = 1.0, double initial_kinetic_energy = 0.0, double shockwave_penetration_depth = 0.0);

	// set parameters that define how to run the simulation -> must be done before calling update()
	ErrorCode startSimulation(double duration, double timestep, int print_energies_interval = 0, int print_positions_interval = 0, const char *energy_filename = NULL, const char *print_positions_path = NULL, bool use_gravity = false);

	// set dynamic stop condition
	void setPotentialVariationStopCondition(double min_time, double potential_variation_stop_threshold, int check_potential_variation_interval);

#ifdef TRACK_CONTACTS_PER_PARTICLE
	void initInitialContacts();
#endif

#ifdef TRACK_PARTICLE_ORIENTATION
	void initParticleOrientation();
#endif

#ifdef ENABLE_FIXED_PARTICLES
	void initFixedPaticles(bool fixed);

	ErrorCode fixateParticle(int particle_id, bool fixed);
#endif

	///////////////////////////////////////////////////////////////////////
	// execution
	///////////////////////////////////////////////////////////////////////

#ifdef TRACK_PARTICLE_ORIENTATION
	void updateParticleOrientation(double timestep);
#endif

	// performs one update cycle 
	virtual void update();

	// calculates new positions
	virtual void predictor(double timestep);

	// updates the velocity/angular velocity
	void corrector(double timestep);

	// damps velocities
	void dampVelocities(double damping_factor);

	// accelerates aggregate to rotate around z-axis
	void induceRotation(double azimuthal_acceleration);

	// checks if new contacts have formed or existing contacts have broken
	void updateSticking();

	// updates contact pointer/twisting displacement
	void updateContacts(double timestep);

	// calculates forces/torques that acts upon each particle with respect to the current positions of particles
	void updateParticleInteraction(double timestep);

	// change pointers - new value are now old values and the former old values will be overwritten with new values in the next integration step
	void switchOldNewValues();

	// returns a pointer to the contact list entry that has been added or NULL if the two particles are already in contact with each other
	ContactListEntry* getContactInsertPos(int id1, int id2);

#ifdef SURFACE_ASPERITIES
	// perform damping caused by surface asperities when two particles collide
	void surfaceAsperitiesDamping(int id1, int id2);
    void surfaceAsperitiesDampingWall(int pid, int wid);
#endif

	// keeps track of energy dissipation due to potentials jumps at contact breaking
	void updateDissipatedBreakingEnergy(Contact *contact, double particle_dist);

	///////////////////////////////////////////////////////////////////////
	// box
	///////////////////////////////////////////////////////////////////////

	// initializes a box with given properties
	void initBox(double height, double x_size, double y_size, double lower_pos_x, double lower_pos_y, double lower_pos_z, double bottom_wall_id, double top_wall_id);

	// update position/height of the box
	void updateBox();

	// returns true if particle is located within current box (make sure box is initialized)
	inline bool particleInBox(int id);

	// determines the min/max coordinates of a box that contains all particles
	ErrorCode getEnclosingBox(vec3 *min_pos, vec3 *max_pos);

	///////////////////////////////////////////////////////////////////////
	// misc
	///////////////////////////////////////////////////////////////////////

	// determines the number of contacts by checking the current contact list
	int getNumberOfContacts();
    int getNumberOfWallContacts();

	// returns current cross section
	double getCrossSection();

	// returns filling factor based on cross section
	double getCrossSectionFillingFactor();

	// returns filling factor based on initial box size
	double getBoxFillingFactor();

	void calculateRotation(double *mean_omega, double *sigma_omega);

	static void getErrorMessage(ErrorCode error_code, char *message);

	// debug dump of all current data
	void debugDump(const char *filename);

	// true if properly initialized
	bool initialized;


	// energy that has dissipated since last init
	double dissipated_contact_energy;
	double dissipated_rolling_energy;
	double dissipated_sliding_energy;
	double dissipated_twisting_energy;
	double dissipated_wall_energy;
	double dissipated_damping_energy;


	// stores information about which particles are in contact with each other/with walls
	ContactListEntry **contact_list;

	// grid for faster determination of neighbouring particles
	SpatialGrid grid;

	// allows fast determination of forces in normal direction
	NormalInteraction normal_interaction;

	int number_of_particles;
	int number_of_contacts;
	int number_of_walls;
    int number_of_wall_contacts;

	// walls
	std::vector<Wall> walls;

	// position of particles
	double *pos_new;		// in cm
	double *pos_old;		// in cm

	// velocity of particles
	double *vel;			// in cm/s

	// angular velocity of particles
	double *vel_angular;

	// total force acting on a particle
    double* force_new;		// in g / (s^2 cm)      // these units descibe pressure   force should be [g * cm / s^2]
	double* force_old;		// in g / (s^2 cm)

	// total torque acting on a particle
	double* torque_new;
	double* torque_old;

#ifdef TRACK_PARTICLE_ORIENTATION
	double *orientation;
#endif

#ifdef ENABLE_FIXED_PARTICLES
	bool *fixed_particles;
#endif

	// stores the current time, the time when simulation is finished and the upate interval
	double current_time;	// in s
	double end_time;		// in s
	double min_end_time;	// in s
	double timestep;		// in s

	// true if simulation has finished (e.g. end_time has been reached)
	bool stop_simulation; 

	// stores various things describing the sim
	SimInfo sim_info;

	// box info for compression experiments
	WallBox *box;

	// keep track of the number of broken/created contacts
    unsigned int broken_contacts;
    unsigned int created_contacts;
	int created_wall_contacts;
	int broken_wall_contacts;

#ifdef TRACK_CONTACTS_PER_PARTICLE
	int *initial_number_of_contacts_of_particle;
	int *broken_contacts_of_particle;
	int *created_contacts_of_particle;
#endif

	// interval in which changes of potential energies are evaluated (allows to abort simulation if nothing happens anymore)
	int check_potential_variation_interval;
	int check_potential_variation_counter;

	double potential_variation_stop_threshold;
	double last_dissipated_energy;

	// counters to print energy/position after desired ammount of integration steps
	int print_energies_counter;
	int print_positions_counter;
	int print_positions_file_counter;

	// specifies after how many timesteps values will be written to log files (0 if none)
	int print_energies_interval;
	int print_positions_interval;

	// files where energies/positions will be stored
	//FILE *file_energies;
	//char path_file_positions[200];
	char *energies_filename;
	char *positions_path;

	// to determine cross section
	int cross_section_cells[GRID_CELLS*GRID_CELLS];

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	double *dissipated_energy_of_particle;
#endif
};

#endif
