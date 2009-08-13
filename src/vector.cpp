/******************************************************************************
 *  _ _ _      ____                                                           *
 * | (_) |__  / ___|_ __ __ _ ___ ___ _ __ ___   __ _ _ __  _ __              *
 * | | | '_ \| |  _| '__/ _` / __/ __| '_ ` _ \ / _` | '_ \| '_ \             *
 * | | | |_) | |_| | | | (_| \__ \__ \ | | | | | (_| | | | | | | |            *
 * |_|_|_.__/ \____|_|  \__,_|___/___/_| |_| |_|\__,_|_| |_|_| |_|            *
 *                                                                            *
 *                                                                            *
 * The files in this directory and elsewhere which refer to this LICENCE      *
 * file are part of Grassmann, the library for the high-level management of   *
 * linear algebra structs.                                                    *
 *                                                                            *
 * Copyright (C) 2009 by BlackLight, <blacklight@autistici.org>               *
 * Web: http://0x00.ath.cx                                                    *
 *                                                                            *
 * Grassmann is free software; you can redistribute it and/or modify it under *
 * the terms of the GNU General Public License as published by the Free       *
 * Software Foundation; either version 3 or (at your option) any later        *
 * version.                                                                   *
 *                                                                            *
 * Grassmann is distributed in the hope that it will be useful, but WITHOUT   *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License      *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with Grassmann; if not, write to the Free Software Foundation, Inc.,       *
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.                     *
 *                                                                            *
 * As a special exception, if other files instantiate templates or use        *
 * macros or inline functions from these files, or you compile these          *
 * files and link them with other works to produce a work based on these      *
 * files, these files do not by themselves cause the resulting work to be     *
 * covered by the GNU General Public License. However the source code for     *
 * these files must still be made available in accordance with section (3)    *
 * of the GNU General Public License.                                         *
 *                                                                            *
 * This exception does not invalidate any other reasons why a work based on   *
 * this file might be covered by the GNU General Public License.              *
 *******************************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cmath>

#include "grassmann.hpp"
#include "grassmann_exception.hpp"

using std::vector;
using std::string;
using std::stringstream;
using std::exception;

namespace grassmann  {
double Vector::norm2()  {
	double res = 0;

	for (size_t i=0; i < v.size(); i++)
		res += v[i]*v[i];

	return sqrt(res);
}

Vector operator+ (Vector& a, Vector& b) throw(UnequalVectorSizeException)  {
	Vector c;
	
	if (a.size() != b.size())
		throw UnequalVectorSizeException();

	for (size_t i=0; i < a.size(); i++)
		c.push_back(a[i] + b[i]);

	return c;
}

Vector operator- (Vector& a, Vector& b) throw(UnequalVectorSizeException)  {
	Vector c;
	
	if (a.size() != b.size())
		throw UnequalVectorSizeException();

	for (size_t i=0; i < a.size(); i++)
		c.push_back(a[i] - b[i]);

	return c;
}

double operator* (Vector& a, Vector& b) throw(UnequalVectorSizeException)  {
	double res = 0;

	if (a.size() != b.size())
		throw UnequalVectorSizeException();

	for (size_t i=0; i < a.size(); i++)
		res += a[i]*b[i];

	return res;
}

Vector operator* (double l, Vector v)  {
	Vector vv(v.size());

	for (size_t i=0; i < v.size(); i++)
		vv[i] = v[i]*l;

	return vv;
}

Vector operator* (Vector v, double l)  {
	return l*v;
}

double& Vector::operator[] (size_t i) throw(IndexOutOfBoundsException)  {
	if (i >= v.size())
		throw IndexOutOfBoundsException();

	return v[i];
}

Vector& Vector::operator<< (double element)  {
	v.push_back(element);
	return *this;
}

bool operator== (Vector a, Vector b)  {
	if (a.size() != b.size())
		return false;

	for (size_t i=0; i < a.size(); i++)
		if (a[i] != b[i])
			return false;

	return true;
}

bool operator!= (Vector a, Vector b)  {
	return !(a == b);
}

void Vector::operator= (string s)  {
	*this = Vector(s,',');
}

size_t Vector::size()  {
	return v.size();
}

double Vector::max()  {
	if (v.size() == 0)
		return 0;

	double M = v[0];

	for (size_t i=0; i<v.size(); i++)  {
		if (v[i]>M)
			M=v[i];
	}

	return M;
}

Vector::Vector()  { v.clear(); }

Vector::Vector (vector<double> v)  {
	this->v = v;
}

Vector::Vector (double *vect, size_t size)  {
	v.clear();

	for (size_t i=0; i < size; i++)
		v.push_back(vect[i]);
}

Vector::Vector (size_t size, ...)  {
	v.clear();
	va_list elements;
	va_start(elements,size);

	for (size_t i=0; i < size; i++)
		v.push_back(va_arg(elements, double));

	va_end(elements);
}

Vector::Vector (string s, char d)  {
	char delim[2];
	delim[0] = d;
	delim[1] = 0;
	v.clear();

	char *tok;
	char *str = new char[s.length()];
	size_t len = 0;

	for (size_t i=0; i < s.length(); i++)  {
		if (s[i] != ' ' && s[i] != '\t' && s[i] != '\r' && s[i] != '\n')
			str[len++] = s[i];
	}

	tok = strtok(str, delim);

	while (tok)  {
		v.push_back(atof(tok));
		tok = strtok(NULL, delim);
	}

	delete [] str;
}

void Vector::push_back (double el)  {
	v.push_back(el);
}

bool Vector::isNull()  {
	for (size_t i=0; i < size(); i++)
		if (v[i]!=0)
			return false;
	return true;
}

bool Vector::multiple(Vector v1)  {
	if (size() != v1.size())
		return false;

	if (size() == 0)
		return true;

	double pivot = v[0]/v1[0];

	for (size_t i=0; i < size(); i++)  {
		if (v1[i]==0)  {
			if (v[i]!=0)
				return false;
		} else {
			if ((v[i]/v1[i]) != pivot)
				return false;
		}
	}
	
	return true;
}

void Vector::clear()  {
	v.clear();
}

string Vector::toString()  {
	stringstream ss (stringstream::in | stringstream::out);

	ss << "(";

	for (size_t i=0; i < size(); i++)  {
		if (i > 0 && i < size())
			ss << " ";

		ss << " " << v[i] << " ";
	}

	ss << ")";
	return ss.str();
}

Vector nullVector (size_t n)  {
	vector<double> v(n, 0.0);
	return Vector(v);
}

double polynomialInterpolation (double x, vector<Vector> points)  {
	size_t n = points.size();
	double px = 0.0;
	vector<double> lagrangian(n, 1.0);

	for (size_t i=0; i < n; i++)
		for (size_t j=0; j < n; j++)
			if (i != j)
				lagrangian[i] *= ( (x - points[j][0]) / (points[i][0] - points[j][0]) );

	for (size_t i=0; i < n; i++)
		px += (points[i][1] * lagrangian[i]);

	return px;
}

namespace interpolateNS  {
	bool vectorCmp (grassmann::Vector a, grassmann::Vector b)  {
		return a[0] < b[0];
	}
}

double linearInterpolation (double x, vector<Vector> points) throw(ValueOutOfRangeException)  {
	int first_index = -1, last_index = -1;
	sort (points.begin(), points.end(), interpolateNS::vectorCmp);

	if ((x < points[0][0]) || (x > points[points.size()-1][0]))
		throw ValueOutOfRangeException();

	for (size_t i=0; i < points.size()-1; i++)  {
		if (x >= points[i][0] && x <= points[i+1][0])  {
			first_index = i;
			last_index  = i+1;
			break;
		}
	}

	double x0 = points[first_index][0],
		  x1 = points[last_index][0],
		  y0 = points[first_index][1],
		  y1 = points[last_index][1];

	return ( (y1-y0)/(x1-x0) )*(x-x0) + y0;
}

// end of grassmann namespace, game over
}

