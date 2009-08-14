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

#include <algorithm>

#include "grassmann.hpp"
#include "grassmann_exception.hpp"

using std::vector;

namespace grassmann  {
	double polynomialInterpolation (double x, vector<Vector> points) throw(CoincidentPointsException)  {
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

	double linearInterpolation (double x, vector<Vector> points)
			throw(ValueOutOfRangeException, CoincidentPointsException)  {
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
}

