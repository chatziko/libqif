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
#include "MinEntropy.h"
#include "aux.h"

template<typename eT>
eT MinEntropy<eT>::vulnerability(const Prob<eT>& pi) {
	eT max_x = eT(0);
	for(uint x = 0; x < pi.n_cols; x++) {
		eT el = pi.at(x);
		if(el > max_x)
			max_x = el;
	}
	return max_x;

//	return accu(max(this->C, 1));
}

//sum y max x pi(x) C[x,y]
template<typename eT>
eT MinEntropy<eT>::cond_vulnerability(const Prob<eT>& pi) {
	this->check_prior(pi);

	eT sum_y = eT(0);

	for(uint y = 0; y < this->C.n_cols; y++) {
		eT max_x = eT(0);
		for(uint x = 0; x < this->C.n_rows; x++) {
			eT el = pi.at(x) * this->C.at(x, y);
			if(el > max_x) {
				max_x = el;
			}
		}
		sum_y += max_x;
	}
	return sum_y;

//	eT s = eT(0);
//	for(uint y = 0; y < this->C.n_cols; y++)
//		s += max(pi % trans(this->C.col(y)));
}

template<typename eT>
eT MinEntropy<eT>::max_mult_leakage() {
	eT sum_y = eT(0);

	for(uint y = 0; y < this->C.n_cols; y++) {
		eT max_x = 0;
		for(uint x = 0; x < this->C.n_rows; x++) {
			eT el = this->C.at(x, y);
			if(el > max_x)
				max_x = el;
		}
		sum_y += max_x;
	}
	return sum_y;
}

template class MinEntropy<double>;
template class MinEntropy<float>;
template class MinEntropy<urat>;

