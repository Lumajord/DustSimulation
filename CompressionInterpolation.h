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

#ifndef NormalInteraction_H
#define NormalInteraction_H

#include "Constants.h"
#include <vector>

class NormalInteraction
{
public:
	// inits table for specified range
#ifdef INTERPOLATE_NORMAL_FORCES
	void init(int interpolation_points);
#else
	void init();
#endif

	// returns approx force for specified compression length
	double getJKRApproxForce(double compression_length);

	// calculates contact radius / force / potential energy for JKR particle<->particle interaction
    double getJKRContactRadius(double compression_length);
    double getJKRForce(double contact_radius);
	double getJKRPotentialEnergy(double compression_length);

	// calculates contact radius / force / potential energy for JKR particle<->particle interaction
	double getJKRWallContactRadius(double compression_length);
    double getJKRWallForce(double contact_radius);
	double getJKRWallPotentialEnergy(double compression_length);

	// calculates contact radius / force for DMT particle<->particle interaction
	double getDMTForce(double compression_length);
	double getDMTContactRadius(double compression_length);
    double getDMTWallForce(double compression_length);
    double getDMTWallContactRadius(double compression_length);

	// returns viscous damping force according to Brilliatov et al. 2007
	double getViscousDissipationForce(double current_contact_radius, double old_contact_radius, double timestep);

//private:
public:
	// number of interpolation points
	int interpolation_points;

	// array storing the forces
	std::vector<double> values;

	double min_compression_length;
	double max_compression_length;

	double interval_width;

	// some constants needed to calculate the JKR forces
	double k1, k2, k3, k4;
	double c1_contact_radius, c2_contact_radius;
	double c1_contact_radius_wall, c2_contact_radius_wall;
	double k1_wall, k2_wall;

	// some constants for Brilliantov
	double br_k1, br_k2;
};

#endif
