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

#ifndef DOUBLE3_H
#define DOUBLE3_H

#if defined _WIN32 || _WIN64
	#define _USE_MATH_DEFINES
	#include <math.h>
#else
	#include <cmath>
#endif

typedef double vec3[3];

inline double dot_product(vec3 &a, vec3 &b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

inline double norm_squared(vec3 &vec) { return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);}

inline double norm(vec3 &vec) {return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]); }//{return sqrtFast(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]); }

inline void normalize(vec3 *vec) { double n = 1.0 / sqrt( (*vec)[0] * (*vec)[0] + (*vec)[1] * (*vec)[1] + (*vec)[2] * (*vec)[2]); (*vec)[0] *= n; (*vec)[1] *= n; (*vec)[2] *= n; }

#endif
