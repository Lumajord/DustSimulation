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

#include "CompressionInterpolation.h"

#include <cstdio>

extern double particle_radius;
extern double reduced_radius;
extern double equilibrium_contact_radius;
extern double surface_energy;
extern double young_mod_reduced;
extern double delta_c;
extern double nu;
extern double young_mod;
extern double viscous_constant;
extern double F_c;
extern double wall_reduced_radius;
extern double wall_equilibrium_contact_radius;
extern double wall_delta_c;
extern double wall_F_c;

#ifdef INTERPOLATE_NORMAL_FORCES
void NormalInteraction::init(int interpolation_points)
#else
void NormalInteraction::init()
#endif
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate constants for JKR interaction
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	k1 = 12.0 * M_PI * surface_energy * reduced_radius / (equilibrium_contact_radius * equilibrium_contact_radius * equilibrium_contact_radius);
	k2 = 12.0 * M_PI * surface_energy * reduced_radius / pow(equilibrium_contact_radius, 1.5);
	k3 = 4.0 / 3.0 * young_mod_reduced / reduced_radius;
	k4 = 4.0 * M_PI * reduced_radius * surface_energy;

	c1_contact_radius = sqrt(equilibrium_contact_radius);
	c2_contact_radius = 0.25 * 2.0 / 3.0 * pow(equilibrium_contact_radius, 1.5);

#ifdef USE_DEFAULT_WALL_NORMAL_INTERACTION
	c1_contact_radius_wall = c1_contact_radius;
	c2_contact_radius_wall = c2_contact_radius;
	k1_wall = k1;
	k2_wall = k2;
#else
	c1_contact_radius_wall = sqrt(wall_equilibrium_contact_radius);
	c2_contact_radius_wall = 0.25 * 2.0 / 3.0 * pow(wall_equilibrium_contact_radius, 1.5);
	k1_wall = 12.0 * M_PI * surface_energy * wall_reduced_radius / (wall_equilibrium_contact_radius * wall_equilibrium_contact_radius * wall_equilibrium_contact_radius);
	k2_wall = 12.0 * M_PI * surface_energy * wall_reduced_radius / pow(wall_equilibrium_contact_radius, 1.5);
#endif

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate constants for Brilliantov interaction
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double D = 1.5 * (1 - nu*nu) / young_mod;

	br_k1 = 3.0 / (D * reduced_radius);
	br_k2 = 1.5 * sqrt(12.0 * M_PI * surface_energy / D);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// init interpolation
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INTERPOLATE_NORMAL_FORCES
	double contact_radius;
	this->interpolation_points = interpolation_points;

	min_compression_length = - delta_c;
	max_compression_length = 0.25 * particle_radius;

	interval_width = (max_compression_length - min_compression_length) / (double)interpolation_points;

	values.resize(interpolation_points+1);

	// calculate values in the desired compression_length interval
	for(int i = 0; i < interpolation_points+1; ++i)
	{
		contact_radius = getJKRContactRadius(min_compression_length + (double)i * interval_width);
		values[i] = ( k1 * contact_radius*contact_radius*contact_radius - k2 * pow(contact_radius, 1.5) );
	}
#endif
}

double NormalInteraction::getJKRApproxForce(double compression_length)
{
#ifdef INTERPOLATE_NORMAL_FORCES
	if(compression_length < max_compression_length)
	{
		double pos;
		double offset = (compression_length - min_compression_length) / interval_width;

		// TODO: better interpolation
		offset = modf(offset, &pos);

		if(pos < 0)
			return 0;
		else
			return values[(int)pos] + offset * (values[(int)pos + 1] - values[(int)pos]);
	}
	else
	{
		double contact_radius = getJKRContactRadius(compression_length);
		return getJKRForce(contact_radius);
	}
#else
	double contact_radius = getJKRContactRadius(compression_length);
	return getJKRForce(contact_radius);
#endif
}

double NormalInteraction::getDMTContactRadius(double compression_length)
{
	if(compression_length < 0)
		return 0;
	else
		return sqrt(compression_length * reduced_radius);
}

double NormalInteraction::getDMTForce(double compression_length)
{
	if(compression_length > 0)
	{
		double a = sqrt(compression_length * reduced_radius);

#ifdef USE_NORMAL_INTERACTION_MODIFIER
		return NORMAL_INTERACTION_MODIFIER * (k3 * a*a*a - k4);
#else
		return k3 * a*a*a - k4;
#endif
	}
	else
		return 0;
}


double NormalInteraction::getDMTWallContactRadius(double compression_length)
{
    if(compression_length < 0)
        return 0;
    else
        return sqrt(compression_length * particle_radius); // for wall contact: reduced radius = particle radius
}


double NormalInteraction::getDMTWallForce(double compression_length)
{
    if(compression_length > 0)
    {
        double a = sqrt(compression_length * particle_radius); // for wall contact: reduced radius = particle radius

#ifdef USE_NORMAL_INTERACTION_MODIFIER
        return NORMAL_INTERACTION_MODIFIER * (k3 * a*a*a - k4);
#else
        return 0.5 * k3 * a*a*a - 2.0 * k4; // adjust for different reduced radius for wall contact: k3 contains reduced radius in denominator and k4 in the numerator
#endif
    }
    else
        return 0;
}


double NormalInteraction::getJKRContactRadius(double compression_length)
{
	// contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
	// use equilibrium contact radius as starting value
	double k = compression_length * reduced_radius / 3.0;
	double x_pow3;
	double x_new;
	double x_old = c1_contact_radius;

    // use Newton-Raphson method to find root
#ifdef USE_FIXED_CONTACT_RADIUS_ACCURACY
    for(int i = 0; i < 20; ++i)
    {
		x_pow3 = x_old * x_old * x_old;
		x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius);

        if(std::abs(x_new - x_old) / x_new < 1.e-14)
            break;

		x_old = x_new;
    }
#else
	int i = 100;
	do
	{
		x_pow3 = x_old * x_old * x_old;
		x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius);

		if(std::abs(x_new - x_old) / particle_radius < 0.0001)
			break;

		x_old = x_new;
	} while (--i > 0);
#endif

	return x_new * x_new;
}

double NormalInteraction::getJKRForce(double contact_radius)
{
	double x = sqrt(contact_radius*contact_radius*contact_radius);

#ifdef USE_NORMAL_INTERACTION_MODIFIER
	return NORMAL_INTERACTION_MODIFIER * x * (k1 * x - k2);
#else
	return x * (k1 * x - k2);
#endif
}

double NormalInteraction::getJKRPotentialEnergy(double compression_length)
{
	double a = getJKRContactRadius(compression_length);
	a /= equilibrium_contact_radius;

#ifdef USE_NORMAL_INTERACTION_MODIFIER
	return NORMAL_INTERACTION_MODIFIER * 4.0 * pow(6.0, 1.0 / 3.0) * F_c * delta_c * (4.0 / 5.0 * a * a * a * a * a - 4.0 / 3.0 * pow(a, 3.5) + 1.0 / 3.0 * a * a);
#else
    return 4.0 * pow(6.0, 1.0 / 3.0) * F_c * delta_c * (4.0 / 5.0 * a*a*a*a*a - 4.0 / 3.0 * a*a*a*sqrt(a) + 1.0 / 3.0 *a*a);
#endif
}

double NormalInteraction::getViscousDissipationForce(double current_contact_radius, double old_contact_radius, double timestep)
{
	double a_dot = (current_contact_radius - old_contact_radius) / timestep;
	return viscous_constant * a_dot * (br_k1 * current_contact_radius*current_contact_radius - br_k2 * sqrt(current_contact_radius));
}


double NormalInteraction::getJKRWallContactRadius(double compression_length)
{
	// contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
	// use equilibrium contact radius as starting value
	double x_new;
	double x_pow3;
	double x_old = c1_contact_radius_wall;
	double k = compression_length * wall_reduced_radius / 3.0;

#ifdef USE_FIXED_CONTACT_RADIUS_ACCURACY
    for(int i = 0; i < 20; ++i)
    {
        x_pow3 = x_old * x_old * x_old;
        x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius_wall);

        if(std::abs(x_new - x_old) / x_new < 1.e-14)
            break;


        x_old = x_new;
    }
#else
	int i = 100;
	do
	{
		x_pow3 = x_old * x_old * x_old;
		x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius_wall);

		if(std::abs(x_new - x_old) / particle_radius < 0.0001)
			break;

		x_old = x_new;
	} while(--i > 0);
#endif

	return x_new * x_new;
}

/*
double NormalInteraction::getJKRWallContactRadius(double compression_length)
{
    // contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
    // use equilibrium contact radius as starting value
    double x_new;
    double x_pow3;
    double x_old = c1_contact_radius_wall;
    double k = compression_length * wall_reduced_radius / 3.0;

    for(int i = 0; i < 20; ++i)
    {
        x_pow3 = x_old * x_old * x_old;
        x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius_wall);

        if(std::abs(x_new - x_old) / x_new < 1.e-14)
            break;

        x_old = x_new;
    }

    return x_new * x_new;
}
*/

double NormalInteraction::getJKRWallForce(double contact_radius)
{

#ifdef USE_NORMAL_INTERACTION_MODIFIER
	return NORMAL_INTERACTION_MODIFIER * (k1_wall * contact_radius*contact_radius*contact_radius - k2_wall * pow(contact_radius, 1.5));
#else
    double x = sqrt(contact_radius*contact_radius*contact_radius);
    return x * (k1_wall * x - k2_wall);
#endif	
}


double NormalInteraction::getJKRWallPotentialEnergy(double compression_length)
{
	double a = getJKRWallContactRadius(compression_length) / wall_equilibrium_contact_radius;

#ifdef USE_NORMAL_INTERACTION_MODIFIER
	return NORMAL_INTERACTION_MODIFIER * 4.0 * pow(6.0, 1.0 / 3.0) * wall_F_c * wall_delta_c * (4.0 / 5.0 * a * a * a * a * a - 4.0 / 3.0 * pow(a, 3.5) + 1.0 / 3.0 * a * a);
#else	
	return 4.0 * pow(6.0, 1.0 / 3.0) * wall_F_c * wall_delta_c * (4.0 / 5.0 * a * a * a * a * a - 4.0 / 3.0 * pow(a, 3.5) + 1.0 / 3.0 * a * a);
#endif
}
