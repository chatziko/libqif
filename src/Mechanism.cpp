#include "Mechanism.h"
/*
This file belongs to the LIBQIF library.
A Quantitative Information Flow C++ Toolkit Library.
Copyright (C) 2013  Universidad Nacional de Río Cuarto(National University of Río Cuarto).
Author: Martinelli Fernán - fmartinelli89@gmail.com - Universidad Nacional de Río Cuarto (Argentina)
LIBQIF Version: 1.0
Date: 12th Nov 2013
========================================================================
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

=========================================================================
*/

#include "aux.h"
#include "Metric.h"


template<typename eT>
bool Mechanism<eT>::is_private(eT epsilon) {
	auto mtv = metrics::mult_total_variation<eT, Prob<eT>>();

	for(uint i = 0; i < C.n_rows; i++) {
		for(uint j = i+1; j < C.n_rows; j++) {
			eT mp = mtv(C.row(i), C.row(j));

			if(!less_than_or_eq(mp, epsilon * d(i, j)))
				return false;
		}
	}
	return true;
}


template<typename eT>
eT Mechanism<eT>::smallest_epsilon() {
	auto mtv = metrics::mult_total_variation<eT, Prob<eT>>();

	eT res(0);
	for(uint i = 0; i < C.n_rows; i++) {
		for(uint j = i+1; j < C.n_rows; j++) {
			eT ratio = mtv(C.row(i), C.row(j)) / d(i, j);
			if(less_than(res, ratio))
				res = ratio;
		}
	}
	return res;
}



template class Mechanism<double>;
template class Mechanism<float>;

//template class Metric<double, uint>;
//template class Metric<float, uint>;
//template class Euclidean<double, uint>;
//template class Euclidean<float, uint>;

