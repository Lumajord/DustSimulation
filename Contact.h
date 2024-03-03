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

#ifndef CONTACT_H
#define CONTACT_H

#include "Double3.h"

class RotParam
{
public:
	RotParam(){e0 = 1.0; e1 = 0.0; e2 = 0.0; e3 = 0.0;};

	double e0, e1, e2, e3;

	void Normalize()
	{
		double norm = e0*e0 + e1*e1 + e2*e2 + e3*e3;
		norm = 1.0 / sqrt(norm);
		e0 *= norm;
		e1 *= norm;
		e2 *= norm;
		e3 *= norm;
	};
};

class Contact
{
public:

	Contact();

	Contact(int id1, int id2, vec3 &n2_initial);

	Contact(int particle_id, int wall_id, vec3 &wall_normal, vec3 &initial_pos);

	void updateRotation(vec3 &omega1, vec3 &omega1_dot, double timestep);

	void updateRotation(vec3 &omega1, vec3 &omega2, vec3 &omega1_dot, vec3 &omega2_dot, double timestep);

	void getCurrentN1(vec3 *n1);
	void getCurrentN2(vec3 *n2);

	void updateN1Initial(vec3 &n1);
	void updateN2Initial(vec3 &n2);

	// ids of the particles forming the contact
	int id1;
	int id2;

	// rotation parameters of the contact pointers
	RotParam rot1;
	RotParam rot2;
	
	// initial contact pointers
	vec3 n1_initial;
	vec3 n2_initial;

	vec3 old_contact_normal;

	// contact radius in the last update step
	//double contact_radius;

	// twisting displacement
	double twisting_displacement;

	// 
	double compression_length;
};

#endif
