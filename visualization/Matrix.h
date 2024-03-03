/***************************************************************************
 *   Matrix class                                                          *
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

#ifndef _MATRIX_H_
#define _MATRIX_H_
#include "math.h"
#include "Vector.h"

/**
	Defines a matrix with m rows and n columns.
*/
template <typename T, int m, int n>
class Matrix
{
	public:
		Matrix();
		~Matrix();

		// entry operators
		inline T& operator[](unsigned int i) { return entries[i]; }
		inline T operator[](unsigned int i) const { return entries[i]; }
		inline T& operator()(unsigned int i, unsigned int j) { return entries[i*n+j]; }
		inline T operator()(unsigned int i, unsigned int j) const { return entries[i*n+j]; }
		inline operator const T*(void) const { return entries; }

		// unary operators
		Matrix operator-();

		// scalar operators
		Matrix& operator*=(T);
		Matrix& operator/=(T);

		Matrix operator+(const Matrix&) const;
		Matrix operator-(const Matrix&) const;
		Matrix& operator+=(const Matrix&);
		Matrix& operator-=(const Matrix&);

		// column / row access
		Vector<T, n> getRow(unsigned int i) const;
		Vector<T, m> getColumn(unsigned int j) const;
		void setRow(unsigned int i, Vector<T, n> &v) const;
		void setColumn(unsigned int j, Vector<T, m> &v) const;

	private:
		T entries[m*n];
};

/**
	constructor
*/
template <typename T, int m, int n>
Matrix<T, m, n>::Matrix(){
	for (unsigned int k = 0; k < m*n; ++k) {
		entries[k] = 0;
	}
}

/**
	destructor
*/
template <typename T, int m, int n>
Matrix<T, m, n>::~Matrix(){
	
}

/**
	unary minus
*/
template <typename T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator-()
{
	Matrix<T, m, n> s;

	for (unsigned int k = 0; k < m*n; ++k) {
		s[k] = -entries[k];
	}

	return s;
}

/**
	scalar multiplication
*/
template <typename T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator*=(T c)
{
	for (unsigned int k = 0; k < m*n; ++k) {
		entries[k] *= c;
	}

	return *this;	
}

/**
	scalar division
*/
template <typename T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator/=(T c)
{
	for (unsigned int k = 0; k < m*n; ++k) {
		entries[k] /= c;
	}

	return *this;	
}

/**
	addition
*/
template <typename T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator+(const Matrix<T, m, n> &r) const
{
	Matrix<T, m, n> s;

	for (unsigned int k = 0; k < m*n; ++k) {
		s[k] = entries[k] + r[k];
	}

	return s;
}

/**
	substraction
*/
template <typename T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator-(const Matrix<T, m, n> &r) const
{
	Matrix<T, m, n> s;

	for (unsigned int k = 0; k < m*n; ++k) {
		s[k] = entries[k] - r[k];
	}

	return s;
}

/**
	addition
*/
template <typename T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator+=(const Matrix<T, m, n> &r)
{
	for (unsigned int k = 0; k < m*n; ++k) {
		entries[k] += r[k];
	}

	return *this;
}

/**
	subtraction
*/
template <typename T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator-=(const Matrix<T, m, n> &r)
{
	for (unsigned int k = 0; k < m*n; ++k) {
		entries[k] -= r[k];
	}

	return *this;
}


/**
	get a row
*/
template <typename T, int m, int n>
Vector<T, n> Matrix<T, m, n>::getRow(unsigned int i) const 
{
	Vector<T, n> v;

	for (unsigned int j = 0; j < n; ++j) {
		v[j] = entries[i*n+j];
	}

	return v;
}

/**
	get a column
*/
template <typename T, int m, int n>
Vector<T, m> Matrix<T, m, n>::getColumn(unsigned int j) const 
{
	Vector<T, m> v;

	for (unsigned int i = 0; i < m; ++i) {
		v[i] = entries[i*n+j];
	}

	return v;
}

/**
	set a row
*/
template <typename T, int m, int n>
void Matrix<T, m, n>::setRow(unsigned int i, Vector<T, n> &v) const
{
	for (unsigned int j = 0; j < n; ++j) {
		entries[i*n+j] = v[j];
	}
}

/**
	set a column
*/
template <typename T, int m, int n>
void Matrix<T, m, n>::setColumn(unsigned int j, Vector<T, m> &v) const
{
	for (unsigned int i = 0; i < m; ++i) {
		entries[i*n+j] = v[i];
	}
}

/**
	matrix-vector multiplication
*/
template <typename T, int m, int n>
Vector<T, m> operator*(const Matrix<T, m, n> &_m, const Vector<T, n> &_v)
{
	Vector<T, m> v;

	for (unsigned int i = 0; i < m; ++i) {
		v[i] = 0.0;
		for (unsigned int j = 0; j < n; ++j) {
			v[i] += _v(j) * _m(i,j);
		}
	}

	return v;
}

/**
	Matrix-matrix multiplication
*/
template <typename T, int m, int n, int o>
Matrix<T, m, n> operator*(const Matrix<T, m, o> &_m1, const Matrix<T, o, n> &_m2)
{
	Matrix<T, m, n> r;

	for (unsigned int i = 0; i < m; ++i) {
		for (unsigned int j = 0; j < n; ++j) {
			r(i,j) = 0;
			for (unsigned int k = 0; k < o; ++k) {
				r(i,j) += _m1(i,k) * _m2(k,j);
			}
		}
	}

	return r;
}

/**
	create rotation matrix
*/
template <typename T>
Matrix<T, 3, 3> rotationMatrix(const Vector<T, 3> &axis, double alpha)
{
	Matrix<T, 3, 3> r;

	double c = cos(alpha);
	double s = sin(alpha);

	r(0,0) = c + axis(0)*axis(0)*(1-c);
	r(0,1) = axis(0)*axis(1)*(1-c)-axis(2)*s;
	r(0,2) = axis(0)*axis(2)*(1-c)+axis(1)*s;

	r(1,0) = axis(1)*axis(0)*(1-c)+axis(2)*s;
	r(1,1) = c + axis(1)*axis(1)*(1-c);
	r(1,2) = axis(1)*axis(2)*(1-c)-axis(0)*s;

	r(2,0) = axis(2)*axis(0)*(1-c)-axis(1)*s;
	r(2,1) = axis(2)*axis(1)*(1-c)+axis(0)*s;
	r(2,2) = c + axis(2)*axis(2)*(1-c);

	return r;
}



#endif