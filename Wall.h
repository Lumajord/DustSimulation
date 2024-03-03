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

#ifndef WALL_H
#define WALL_H

#include "Constants.h"

extern double particle_radius;

struct GyrationRadiusInfo
{
	float radius;
	float pos[3];
	float color[4];
};

class Wall
{
public:
	Wall(void);
	Wall(vec3 &pos, vec3 &normal);
	~Wall(void);

	// sets the spatial extension and orientation of the wall (only for visualization)
	void setSize(vec3 &x_dir, vec3 &y_dir, double x_size, double y_size);

	// calculates the edges of the wall using its current pos and size (only for visualization)
	void updateEdges();

    double getDistanceTo(const vec3 &point);

    double getDistanceTo(const double x, const double y, const double z) const;

	double getWallForce();

	double getBoxVolumeForce();

	vec3 pos;

	vec3 velocity;

	vec3 normal;

	vec3 total_force;
	vec3 total_torque;
	vec3 box_volume_force;  // force exerted on the wall by only those particles which are located within the box

	// modifiers for wall interaction
	double compression_modifier;
	double rolling_modifier;
	double sliding_modifier;
	double twisting_modifier;

	// if mass > 0, the equation of motion of the wall is solved
	double mass;

	double dissipated_energy;

	// true if force/avg force on this wall should be tracked
	bool track_force;	

	//////////////////////////////////////////////////////////////////////////////////////////
	// the following variables only affect the visualization
	//////////////////////////////////////////////////////////////////////////////////////////

	double edges[24];

	double x_size;
	double y_size;

	vec3 x_dir;
	vec3 y_dir;

	// between 0 (transparent) and 1 (opaque)
	float alpha;
};

#endif
