/***************************************************************************
 *   Vector class                                                          *
 *   Copyright (C) 2010 by Tobias Müller                                   *
 *   Tobias_Mueller@twam.info                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _VECTOR_H_
#define _VECTOR_H_
#include "math.h"
#include <cstdarg>

/**
	Defines a vector with dimension dim
*/
template <typename T, int dim>
class Vector
{
	public:
		// contructor/destructor
		Vector();
		Vector(const T(&)[dim]);
		Vector(unsigned int, ...);
		~Vector();

		// entry operators
		//inline T& operator[](unsigned int i) { return entries[i]; }
		//inline T operator[](unsigned int i) const { return entries[i]; }
		inline T& operator()(unsigned int i) { return entries[i]; }
		inline T operator()(unsigned int i) const { return entries[i]; }
		inline operator const T*(void) const { return entries; }
		inline operator T*(void) { return entries; }

		// unary operators
		Vector operator-();

		// scalar operators
		template <typename T_, int dim_> friend Vector<T_, dim_> operator*(T_, const Vector<T_, dim_>&);
		template <typename T_, int dim_> friend Vector<T_, dim_> operator*(const Vector<T_, dim_>&, T_);
		template <typename T_, int dim_> friend Vector<T_, dim_> operator/(const Vector<T_, dim_>&, T_);
		Vector& operator*=(T);
		Vector& operator/=(T);

		// vector operators
		Vector operator+(const Vector&) const;
		Vector operator-(const Vector&) const;
		T operator*(const Vector&) const;
		Vector<T, 3> operator^(const Vector<T ,3>&) const;
		Vector& operator+=(const Vector&);
		Vector& operator-=(const Vector&);
		Vector<T ,3>& operator^=(const Vector<T, 3>&);

		// other stuff
		inline T norm() const { return sqrt(norm2()); };
		inline T norm2() const { T result = 0; for (unsigned int i = 0; i < dim; ++i) { result += entries[i]*entries[i]; } return result; };
		void normalize();

	private:
		T entries[dim];
};

/**
	default constructor
*/
template <typename T, int dim>
Vector<T, dim>::Vector()
{
	for (unsigned int i = 0; i < dim; ++i) {
		entries[i] = 0;
	}
}

/**
	initialize with array
*/
template <typename T, int dim>
Vector<T, dim>::Vector(const T(&data)[dim])
{
	for (unsigned int i = 0; i < dim; ++i) {
		entries[i] = data[i];
	}	
}

/**
	initialize with variable argument list. first argument must be number of following arguments
*/
template <typename T, int dim>
Vector<T, dim>::Vector(unsigned int N, ...)
{
	va_list arg_list;
	va_start(arg_list, N);
	for (unsigned int i = 0; i < (N < dim ? N : dim); ++i) {
		entries[i] = va_arg(arg_list, T);
	}
	va_end(arg_list);
}

/**
	destructor
*/
template <typename T, int dim>
Vector<T, dim>::~Vector()
{
	
}

/**
	unary minus
*/
template <typename T, int dim>
Vector<T, dim> Vector<T, dim>::operator-()
{
	Vector<T, dim> s;

	for (unsigned int i = 0; i < dim; ++i) {
		s[i] = -entries[i];
	}

	return s;
}

/**
	scalar multiplication
*/
template <typename T, int dim>
Vector<T, dim> operator*(T c, const Vector<T, dim> &r)
{
	Vector<T, dim> s;

	for (unsigned int i = 0; i < dim; ++i) {
		s[i] = c*r[i];
	}

	return s;
}

/**
	scalar multiplication
*/
template <typename T, int dim>
Vector<T, dim> operator*(const Vector<T, dim> &r, T c)
{
	Vector<T, dim> s;

	for (unsigned int i = 0; i < dim; ++i) {
		s[i] = r[i]*c;
	}

	return s;
}

/**
	scalar division
*/
template <typename T, int dim>
Vector<T, dim> operator/(const Vector<T, dim> &r, T c)
{
	Vector<T, dim> s;

	for (unsigned int i = 0; i < dim; ++i) {
		s[i] = r[i]/c;
	}

	return s;
}

/**
	scalar multiplication
*/
template <typename T, int dim>
Vector<T, dim>& Vector<T, dim>::operator*=(T c)
{
	for (unsigned int i = 0; i < dim; ++i) {
		entries[i] *= c;
	}

	return *this;
}

/**
	scalar division
*/
template <typename T, int dim>
Vector<T, dim>& Vector<T, dim>::operator/=(T c)
{
	for (unsigned int i = 0; i < dim; ++i) {
		entries[i] /= c;
	}

	return *this;
}

/**
	addition
*/
template <typename T, int dim>
Vector<T, dim> Vector<T, dim>::operator+(const Vector<T, dim> &r) const
{
	Vector<T, dim> s;

	for (unsigned int i = 0; i < dim; ++i) {
		s[i] = entries[i] + r[i];
	}

	return s;
}

/**
	substraction
*/
template <typename T, int dim>
Vector<T, dim> Vector<T, dim>::operator-(const Vector<T, dim> &r) const
{
	Vector<T, dim> s;

	for (unsigned int i = 0; i < dim; ++i) {
		s[i] = entries[i] - r[i];
	}

	return s;
}

/**
	(inner) dot product
*/
template <typename T, int dim>
T Vector<T, dim>::operator*(const Vector<T, dim> &r) const
{
	T result;
	
	for (unsigned int i = 0; i < dim; ++i) {
		result += entries[i]*r[i];
	}
	
	return result;
}

/**
	addition
*/
template <typename T, int dim>
Vector<T, dim>& Vector<T, dim>::operator+=(const Vector<T, dim> &r)
{
	for (unsigned int i = 0; i < dim; ++i) {
		entries[i] += r[i];
	}

	return *this;
}

/**
	subtraction
*/
template <typename T, int dim>
Vector<T, dim>& Vector<T, dim>::operator-=(const Vector<T, dim> &r)
{
	for (unsigned int i = 0; i < dim; ++i) {
		entries[i] -= r[i];
	}

	return *this;
}

/**
	normalize vector
*/
template <typename T, int dim>
void Vector<T, dim>::normalize()
{
	double n = norm();

	// we can't normalize if |r|=0
	if (n == 0)
		return;

	for (unsigned int i = 0; i < dim; ++i) {
		entries[i] /= n;
	}
}

/**
	(exterior) cross product
*/
template <typename T, int dim>
Vector<T, 3> Vector<T, dim>::operator^(const Vector<T, 3> &r) const
{
	Vector<T, 3> s;

	s[0] = entries[1] * r[2] - entries[2] * r[1];
	s[1] = entries[2] * r[0] - entries[0] * r[2];
	s[2] = entries[0] * r[1] - entries[1] * r[0];

	return s;
}

/**
        (exterior) cross product
*/
template <typename T, int dim>
Vector<T, 3>& Vector<T, dim>::operator^=(const Vector<T, 3> &r)
{
	double x=entries[0], y=entries[1], z=entries[2];

	entries[0] = y * r[2] - z * r[1];
	entries[1] = z * r[0] - x * r[2];
	entries[2] = x * r[1] - y * r[0];

	return *this;
}

#endif