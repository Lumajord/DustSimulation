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

#include "Contact.h"

Contact::Contact(int id1, int id2, vec3 &n2_initial)
{
	this->id1 = id1;
	this->id2 = id2;

	n1_initial[0] = -n2_initial[0];
	n1_initial[1] = -n2_initial[1];
	n1_initial[2] = -n2_initial[2];

	this->n2_initial[0] = n2_initial[0];
	this->n2_initial[1] = n2_initial[1];
	this->n2_initial[2] = n2_initial[2];

	old_contact_normal[0] = n2_initial[0];
	old_contact_normal[1] = n2_initial[1];
	old_contact_normal[2] = n2_initial[2];

    compression_length = 0;
	twisting_displacement = 0;
}

Contact::Contact(int particle_id, int wall_id, vec3 &wall_normal, vec3 &initial_pos)
{
	id1 = particle_id;
	id2 = -wall_id-1;

	n1_initial[0] = -wall_normal[0];
	n1_initial[1] = -wall_normal[1];
	n1_initial[2] = -wall_normal[2];

	// n2 stores the location, where the contact has been established
	n2_initial[0] = initial_pos[0];
	n2_initial[1] = initial_pos[1];
    n2_initial[2] = initial_pos[2];


    old_contact_normal[0] = 0.0; // unused for wall contacts
    old_contact_normal[1] = 0.0;
    old_contact_normal[2] = 0.0;


    compression_length = 0;
	twisting_displacement = 0;
}

Contact::Contact()
{
	compression_length = 0;
	twisting_displacement = 0;
}

void Contact::updateRotation(vec3 &omega1, vec3 &omega2, vec3 &omega1_dot, vec3 &omega2_dot, double timestep)
{
	double e_dot[4];
	double e_ddot[4];

	// first particle
	e_dot[0] = - 0.5 * (rot1.e1 * omega1[0] + rot1.e2 * omega1[1] + rot1.e3 * omega1[2]);
	e_dot[1] = 0.5 * (rot1.e0 * omega1[0] - rot1.e2 * omega1[2] + rot1.e3 * omega1[1]);
	e_dot[2] = 0.5 * (rot1.e0 * omega1[1] - rot1.e3 * omega1[0] + rot1.e1 * omega1[2]);
	e_dot[3] = 0.5 * (rot1.e0 * omega1[2] - rot1.e1 * omega1[1] + rot1.e2 * omega1[0]);

    double temp = 0.5 * e_dot[0];

	e_ddot[0] = - 0.25 * (rot1.e0 * norm_squared(omega1) + 2.0 * (rot1.e1 * omega1_dot[0] + rot1.e2 * omega1_dot[1] + rot1.e3 * omega1_dot[2]));
	e_ddot[1] = temp * omega1[0] + 0.5 * (rot1.e0 * omega1_dot[0] - rot1.e2 * omega1_dot[2] + rot1.e3 * omega1_dot[1]);
	e_ddot[2] = temp * omega1[1] + 0.5 * (rot1.e0 * omega1_dot[1] - rot1.e3 * omega1_dot[0] + rot1.e1 * omega1_dot[2]);
	e_ddot[3] = temp * omega1[2] + 0.5 * (rot1.e0 * omega1_dot[2] - rot1.e1 * omega1_dot[1] + rot1.e2 * omega1_dot[0]);

	rot1.e0 += timestep * e_dot[0] + 0.5 * timestep * timestep * e_ddot[0];
	rot1.e1 += timestep * e_dot[1] + 0.5 * timestep * timestep * e_ddot[1];
	rot1.e2 += timestep * e_dot[2] + 0.5 * timestep * timestep * e_ddot[2];
	rot1.e3 += timestep * e_dot[3] + 0.5 * timestep * timestep * e_ddot[3];

	// second particle
	e_dot[0] = - 0.5 * (rot2.e1 * omega2[0] + rot2.e2 * omega2[1] + rot2.e3 * omega2[2]);
	e_dot[1] = 0.5 * (rot2.e0 * omega2[0] - rot2.e2 * omega2[2] + rot2.e3 * omega2[1]);
	e_dot[2] = 0.5 * (rot2.e0 * omega2[1] - rot2.e3 * omega2[0] + rot2.e1 * omega2[2]);
	e_dot[3] = 0.5 * (rot2.e0 * omega2[2] - rot2.e1 * omega2[1] + rot2.e2 * omega2[0]);

    temp = 0.5 * e_dot[0];

	e_ddot[0] = - 0.25 * (rot2.e0 * norm_squared(omega2) + 2.0 * (rot2.e1 * omega2_dot[0] + rot2.e2 * omega2_dot[1] + rot2.e3 * omega2_dot[2]));
	e_ddot[1] = temp * omega2[0] + 0.5 * (rot2.e0 * omega2_dot[0] - rot2.e2 * omega2_dot[2] + rot2.e3 * omega2_dot[1]);
	e_ddot[2] = temp * omega2[1] + 0.5 * (rot2.e0 * omega2_dot[1] - rot2.e3 * omega2_dot[0] + rot2.e1 * omega2_dot[2]);
	e_ddot[3] = temp * omega2[2] + 0.5 * (rot2.e0 * omega2_dot[2] - rot2.e1 * omega2_dot[1] + rot2.e2 * omega2_dot[0]);

	rot2.e0 += timestep * e_dot[0] + 0.5 * timestep * timestep * e_ddot[0];
	rot2.e1 += timestep * e_dot[1] + 0.5 * timestep * timestep * e_ddot[1];
	rot2.e2 += timestep * e_dot[2] + 0.5 * timestep * timestep * e_ddot[2];
	rot2.e3 += timestep * e_dot[3] + 0.5 * timestep * timestep * e_ddot[3];

	// make sure that e0^2 + .. + e3^2 = 1
	rot1.Normalize();
	rot2.Normalize();
}

#include <stdio.h>

void Contact::updateRotation(vec3 &omega1, vec3 &omega1_dot, double timestep)
{
	double e_dot[4];
	double e_ddot[4];


	// first particle
    e_dot[0] =-0.5 * (rot1.e1 * omega1[0] + rot1.e2 * omega1[1] + rot1.e3 * omega1[2]);
	e_dot[1] = 0.5 * (rot1.e0 * omega1[0] - rot1.e2 * omega1[2] + rot1.e3 * omega1[1]);
	e_dot[2] = 0.5 * (rot1.e0 * omega1[1] - rot1.e3 * omega1[0] + rot1.e1 * omega1[2]);
	e_dot[3] = 0.5 * (rot1.e0 * omega1[2] - rot1.e1 * omega1[1] + rot1.e2 * omega1[0]);

	double temp = - 0.25 * (rot1.e1 * omega1[0] + rot1.e2 * omega1[1] + rot1.e3 * omega1[2]);

	e_ddot[0] = - 0.25 * (rot1.e0 * norm_squared(omega1) + 2.0 * (rot1.e1 * omega1_dot[0] + rot1.e2 * omega1_dot[1] + rot1.e3 * omega1_dot[2]));
	e_ddot[1] = temp * omega1[0] + 0.5 * (rot1.e0 * omega1_dot[0] - rot1.e2 * omega1_dot[2] + rot1.e3 * omega1_dot[1]);
	e_ddot[2] = temp * omega1[1] + 0.5 * (rot1.e0 * omega1_dot[1] - rot1.e3 * omega1_dot[0] + rot1.e1 * omega1_dot[2]);
	e_ddot[3] = temp * omega1[2] + 0.5 * (rot1.e0 * omega1_dot[2] - rot1.e1 * omega1_dot[1] + rot1.e2 * omega1_dot[0]);


	rot1.e0 += timestep * e_dot[0] + 0.5 * timestep * timestep * e_ddot[0];
	rot1.e1 += timestep * e_dot[1] + 0.5 * timestep * timestep * e_ddot[1];
	rot1.e2 += timestep * e_dot[2] + 0.5 * timestep * timestep * e_ddot[2];
	rot1.e3 += timestep * e_dot[3] + 0.5 * timestep * timestep * e_ddot[3];


	// make sure that e0^2 + .. + e3^2 = 1
	rot1.Normalize();
}

void Contact::getCurrentN1(vec3 *n1)
{
	(*n1)[0] = 2.0 * ( (rot1.e0 * rot1.e0 +  rot1.e1 * rot1.e1 - 0.5) * n1_initial[0] + (rot1.e1 * rot1.e2 - rot1.e3 * rot1.e0) * n1_initial[1] + (rot1.e1 * rot1.e3 + rot1.e2 * rot1.e0) * n1_initial[2]);
	(*n1)[1] = 2.0 * ( (rot1.e1 * rot1.e2 + rot1.e3 * rot1.e0) * n1_initial[0] + (rot1.e0 * rot1.e0 + rot1.e2 * rot1.e2 - 0.5) * n1_initial[1] + (rot1.e2 * rot1.e3 - rot1.e1 * rot1.e0) * n1_initial[2]);
	(*n1)[2] = 2.0 * ( (rot1.e1 * rot1.e3 - rot1.e2 * rot1.e0) * n1_initial[0] + (rot1.e2 * rot1.e3 + rot1.e1 * rot1.e0) * n1_initial[1] + (rot1.e0 * rot1.e0 +  rot1.e3 * rot1.e3 - 0.5) * n1_initial[2]);
}


void Contact::getCurrentN2(vec3 *n2)
{
	(*n2)[0] = 2.0 * ( (rot2.e0 * rot2.e0 +  rot2.e1 * rot2.e1 - 0.5) * n2_initial[0] + (rot2.e1 * rot2.e2 - rot2.e3 * rot2.e0) * n2_initial[1] + (rot2.e1 * rot2.e3 + rot2.e2 * rot2.e0) * n2_initial[2]);
	(*n2)[1] = 2.0 * ( (rot2.e1 * rot2.e2 + rot2.e3 * rot2.e0) * n2_initial[0] + (rot2.e0 * rot2.e0 + rot2.e2 * rot2.e2 - 0.5) * n2_initial[1] + (rot2.e2 * rot2.e3 - rot2.e1 * rot2.e0) * n2_initial[2]);
	(*n2)[2] = 2.0 * ( (rot2.e1 * rot2.e3 - rot2.e2 * rot2.e0) * n2_initial[0] + (rot2.e2 * rot2.e3 + rot2.e1 * rot2.e0) * n2_initial[1] + (rot2.e0 * rot2.e0 +  rot2.e3 * rot2.e3 - 0.5) * n2_initial[2]);
}

void Contact::updateN1Initial(vec3 &n1)
{
	n1_initial[0] = 2.0 * ( (0.5 - rot1.e2 * rot1.e2 - rot1.e3 * rot1.e3) * n1[0] + (rot1.e1 * rot1.e2 + rot1.e3 * rot1.e0) * n1[1] + (rot1.e1 * rot1.e3 - rot1.e2 * rot1.e0) * n1[2]);
	n1_initial[1] = 2.0 * ( (rot1.e1 * rot1.e2 - rot1.e3 * rot1.e0) * n1[0] + (0.5 - rot1.e1 * rot1.e1 - rot1.e3 * rot1.e3) * n1[1] + (rot1.e2 * rot1.e3 + rot1.e1 * rot1.e0) * n1[2]);
	n1_initial[2] = 2.0 * ( (rot1.e1 * rot1.e3 + rot1.e2 * rot1.e0) * n1[0] + (rot1.e2 * rot1.e3 - rot1.e1 * rot1.e0) * n1[1] + (0.5 - rot1.e1 * rot1.e1 - rot1.e2 * rot1.e2) * n1[2]);
}

void Contact::updateN2Initial(vec3 &n2)
{
	n2_initial[0] = 2.0 * ( (0.5 - rot2.e2 * rot2.e2 - rot2.e3 * rot2.e3) * n2[0] + (rot2.e1 * rot2.e2 + rot2.e3 * rot2.e0) * n2[1] + (rot2.e1 * rot2.e3 - rot2.e2 * rot2.e0) * n2[2]);
	n2_initial[1] = 2.0 * ( (rot2.e1 * rot2.e2 - rot2.e3 * rot2.e0) * n2[0] + (0.5 - rot2.e1 * rot2.e1 - rot2.e3 * rot2.e3) * n2[1] + (rot2.e2 * rot2.e3 + rot2.e1 * rot2.e0) * n2[2]);
	n2_initial[2] = 2.0 * ( (rot2.e1 * rot2.e3 + rot2.e2 * rot2.e0) * n2[0] + (rot2.e2 * rot2.e3 - rot2.e1 * rot2.e0) * n2[1] + (0.5 - rot2.e1 * rot2.e1 - rot2.e2 * rot2.e2) * n2[2]);
}
