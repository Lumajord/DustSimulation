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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Double3.h"

#ifdef _MSC_VER
#pragma warning(disable: 4005 4018 4244 4267 4996) // signed/unsigned, loss of precision, deprecation...
#endif

#define CORE_VERSION "1.21"

const char DATA_FILE_VERSION[] = "DATA_FILE_VERSION_1_10";

enum SimType {SIM_TYPE_GENERAL, SIM_TYPE_COLLISION, SIM_TYPE_WALL_COLLISION, SIM_TYPE_COMPRESSION_NO_SIDE_WALLS, SIM_TYPE_COMPRESSION_WITH_SIDE_WALLS, SIM_TYPE_COMPRESSION_RELAXATION, SIM_TYPE_SHOCKWAVE, SIM_TYPE_DYNAMIC_COMPRESSION, SIM_TYPE_COMPACTION, SIM_TYPE_OPEN_BOX, SIM_TYPE_PULL_STRENGTH_TEST, SIM_TYPE_SHEAR_STRENGTH_TEST};

enum DisplayMode {DISPLAY_NOTHING, DISPLAY_PARTICLE_SPEED, DISPLAY_DISSIPATED_ENERGY, DISPLAY_DENSITY, DISPLAY_PARTICLE_CONTACTS, DISPLAY_WALL_CONTACTS, DISPLAY_WALL_ANGLE, DISPLAY_FRAGMENTS, DISPLAY_DISLOCATION, DISPLAY_50_50};

enum ErrorCode {EC_OK, EC_INTERNAL_ERROR, EC_FILE_NOT_FOUND, EC_INVALID_FILE_VERSION, EC_INVALID_MATERIAL_FILE_VERSION, EC_NO_BOX, EC_INVALID_BOX_TYPE, EC_SIM_NOT_INITIALIZED, EC_NO_PARTICLES, EC_DIFFERING_MATERIAL, EC_INVALID_PARAMETER, EC_NO_CUDA_DEVICE, EC_CUDA_DEVICE_UNAVAILABLE, EC_INSUFFICIENT_GPU_MEMORY, EC_CUDA_NOT_INITIALIZED, EC_CUDA_INVALID_BOUNDARIES, EC_INVALID_BOUNDARIES, EC_TIMESTEP_TOO_LARGE, EC_NO_HIT_AND_STICK_POS_FOUND, EC_FIXED_PARTICLES_NOT_INITIALIZED, EC_CONTACT_LIST_CORRUPTED};

enum CollisionResult {COLLISION_RESULT_STICKING, COLLISION_RESULT_BOUNCING, COLLISION_RESULT_FRAGMENTATION};

enum BAMSelectionMethod {BAM_SELECT_RANDOM, BAM_SELECT_CLOSEST, BAM_SELECT_INTERIOR};

// addressing scheme for vectors
#define X_COORD(index) (3 * (index))
#define Y_COORD(index) (3 * (index) + 1)
#define Z_COORD(index) (3 * (index) + 2)

#define X_COORD_(index) (3 * (index))
#define Y_COORD_(index) (3 * (index) + 1)
#define Z_COORD_(index) (3 * (index) + 2)

// converts contact wall id to wall id
#define WALL_ID(index) (-index-1)

////////////////////////////////////////////////////////////////////////////////////////
// the following lines specify which things are taken into account
////////////////////////////////////////////////////////////////////////////////////////

// if enabled, the rolling, sliding, and twisting coefficients depend on the current contact radius (as in D&T1997)
//#define USE_VARYING_POTENTIAL_COEFFICIENTS

#ifdef USE_VARYING_POTENTIAL_COEFFICIENTS
	#define K_R k_r2
	#define K_S k_s2
	#define K_T k_t2
#else
	#define K_R k_r
	#define K_S k_s
	#define K_T k_t
#endif

// select type of normal interaction (one of these options must be chosen)
#define USE_JKR
//#define USE_DMT

// modifies the strength of the normal interaction
//#define USE_NORMAL_INTERACTION_MODIFIER
const double NORMAL_INTERACTION_MODIFIER = 1.0;

// interpolate strength of normal forces from precalculated table
//#define INTERPOLATE_NORMAL_FORCES
//const int NORMAL_FORCES_INTERPOLATION_POINTS = 2048;

// use fixed number of iterations to determine contact radius (for comparion with GPU)
#define USE_FIXED_CONTACT_RADIUS_ACCURACY

// if enabled, wall normal interaction is calculated using the particle<->particle interaction (-> different reduced_radius is not taken into account)
//#define USE_DEFAULT_WALL_NORMAL_INTERACTION

// enables enhanced gluing of particles that are close to walls
//#define ENABLE_WALL_GLUE

// distance (in particle radii) for additional gluing
const double wall_glue_distance = 10.0;

// modifier of normal force of particles that are close to walls
const double wall_glue_strength = 3.0;

// turn off different kinds of particle <-> particle interactions


#define ROLLING
#define SLIDING
#define TWISTING

#define INELASTIC_ROLLING
#define INELASTIC_SLIDING
#define INELASTIC_TWISTING

// turn off different kinds of particle <-> wall interactions

//#define SURFACE_ASPERITIES // additional energy dissipation mechanism. See Güttler et al. 2009


#define WALL_ROLLING
#define WALL_SLIDING
#define WALL_TWISTING

#define INELASTIC_WALL_ROLLING
#define INELASTIC_WALL_SLIDING
#define INELASTIC_WALL_TWISTING


// enable wall damping method
//#define DAMP_WALL_OSCILLATIONS

// enable damping method: simple (weak) damping of normal oscillations or viscoelastic damping proposed by Sebastian Krijt
//#define USE_CONSTANT_DAMPING
#define USE_VISCOELASTIC_DAMPING

// if enabled, the critical rolling displacement is calculated via xi_crit = k * a_eq/12 where k is obtained by taking the value that otherwise defines the cirtical rolling displacement in th material file
//#define USE_KRIJT_CRITICAL_ROLLING_DISPLACEMENT



#define TRACK_DISSIPATED_ENERGY

////////////////////// REQUIRED FOR PLOTTING
// enables tracking of the current rotation state for every particle (only affects visualization)

#define TRACK_PARTICLE_ORIENTATION

// turn off tracking of dissipated energy (saves some computing time if it is not needed)
#define TRACK_DISSIPATED_ENERGY

#define TRACK_DISSIPATED_ENERGY_PER_PARTICLE

// enables gravity (e.g. for compression simulations)
#define ENABLE_GRAVITY

////////////////////// END REQUIRED FOR PLOTTING

// track number of broken contacts per pacrticle
//#define TRACK_CONTACTS_PER_PARTICLE


// enables support for immobile particles
//#define ENABLE_FIXED_PARTICLES


// track forces for the contact number 0 and print them to file
//#define TRACK_FORCES


#define GRID_CELLS 256

// cpu IDs of side walls
#define LEFT_WALL_ID 1
#define RIGHT_WALL_ID 2
#define FRONT_WALL_ID 3
#define BACK_WALL_ID 4

#endif
