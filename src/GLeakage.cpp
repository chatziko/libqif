#include "GLeakage.h"
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

// max_w sum_x pi[x] G[w, x]
//
template<typename eT>
eT GLeakage<eT>::vulnerability(const Prob<eT>& pi) {
	this->check_prior(pi, true);

	return arma::max(this->G * trans(pi));
}

// sum_y max_w sum_x pi[x] C[x, y] G[w, x]
//
template<typename eT>
eT GLeakage<eT>::cond_vulnerability(const Prob<eT>& pi) {
	this->check_prior(pi);

	eT s = eT(0);
	for(uint y = 0; y < this->C.n_cols; y++)
		s += arma::max(this->G * (trans(pi) % this->C.col(y)));
	return s;
}

template class GLeakage<double>;
template class GLeakage<float>;
template class GLeakage<rat>;

