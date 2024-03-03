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

#include "Wall.h"
#include <string.h>

extern double particle_radius;

Wall::Wall(vec3 &pos, vec3 &normal)
{
	memcpy(this->pos, pos, sizeof(vec3));
	memcpy(this->normal, normal, sizeof(vec3));

	Wall();
}

Wall::Wall(void)
{
	total_force[0] = total_force[1] = total_force[2] = 0;
	total_torque[0] = total_torque[1] = total_torque[2] = 0;
	box_volume_force[0] = box_volume_force[1] = box_volume_force[2] = 0;
	velocity[0] = velocity[1] = velocity[2] = 0;
	mass = 0;

	compression_modifier = 1.0;
	rolling_modifier = 1.0;
	sliding_modifier = 1.0;
	twisting_modifier = 1.0;

	dissipated_energy = 0;
	alpha = 0.0f;
	x_size = 0;
	y_size = 0;
	x_dir[0] = x_dir[1] = x_dir[2] = 0;
	y_dir[0] = y_dir[1] = y_dir[2] = 0;

	track_force = false;
}

Wall::~Wall(void)
{
}

void Wall::setSize(vec3 &x_dir, vec3 &y_dir, double x_size, double y_size)
{
	this->x_size = x_size;
	this->y_size = y_size;

	memcpy(this->x_dir, x_dir, sizeof(vec3));
	memcpy(this->y_dir, y_dir, sizeof(vec3));

	updateEdges();
}

void Wall::updateEdges()
{
	// upper edges
	edges[0] = pos[0] + particle_radius * normal[0];
	edges[1] = pos[1] + particle_radius * normal[1];
	edges[2] = pos[2] + particle_radius * normal[2];

	edges[3] = edges[0] + x_size * x_dir[0];
	edges[4] = edges[1] + x_size * x_dir[1];
	edges[5] = edges[2] + x_size * x_dir[2];

	edges[6] = edges[3] + y_size * y_dir[0];
	edges[7] = edges[4] + y_size * y_dir[1];
	edges[8] = edges[5] + y_size * y_dir[2];

	edges[9] = edges[0] + y_size * y_dir[0];
	edges[10] = edges[1] + y_size * y_dir[1];
	edges[11] = edges[2] + y_size * y_dir[2];

	// lower edges
	edges[12] = pos[0] - particle_radius * normal[0];
	edges[13] = pos[1] - particle_radius * normal[1];
	edges[14] = pos[2] - particle_radius * normal[2];

	edges[15] = edges[12] + x_size * x_dir[0];
	edges[16] = edges[13] + x_size * x_dir[1];
	edges[17] = edges[14] + x_size * x_dir[2];

	edges[18] = edges[15] + y_size * y_dir[0];
	edges[19] = edges[16] + y_size * y_dir[1];
	edges[20] = edges[17] + y_size * y_dir[2];

	edges[21] = edges[12] + y_size * y_dir[0];
	edges[22] = edges[13] + y_size * y_dir[1];
	edges[23] = edges[14] + y_size * y_dir[2];
}

double Wall::getDistanceTo(const vec3 &point)
{
	double dist = (point[0] - pos[0]) * normal[0] + (point[1] - pos[1]) * normal[1] + (point[2] - pos[2]) * normal[2];
	return fabs(dist);
}

double Wall::getDistanceTo(const double x, const double y, const double z) const
{
	return fabs( (x-pos[0]) * normal[0] + (y-pos[1]) * normal[1] + (z-pos[2]) * normal[2] );
}

double Wall::getWallForce()
{
	return fabs( dot_product(normal, total_force) );
}

double Wall::getBoxVolumeForce()
{
	return fabs( dot_product(normal, box_volume_force) );
}
