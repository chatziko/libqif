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

#include "Shannon.h"
#include <armadillo>
#include <cmath>

// H(X) = - sum_x pi[x] log2(pi[x])
//
template<typename eT>
double Shannon<eT>::entropy(const Prob<eT>& pi) {
	double sum_x = 0;
	for(uint x = 0; x < pi.n_cols; x++) {
		double el = pi.at(x);
		sum_x -= el > 0 ? el * log2(el) : 0;
	}

	return sum_x;
}

// computes H(X|Y)
//
// we use the formula
//    H(X|Y) = H(Y|X) + H(X) - H(Y)
// since H(Y|X) is easier to compute:
//    H(Y|X) = sum_x pi[x] H(C[x,-])   (entropy of row x)
//
template<typename eT>
double Shannon<eT>::cond_entropy(const Prob<eT>& pi) {
	this->check_prior(pi);

	double Hyx = 0;
	for(uint x = 0; x < this->C.n_rows; x++)
		Hyx += pi.at(x) * entropy(this->C.row(x));

	return Hyx + entropy(pi) - entropy(pi * this->C);
}

//Blahut-Arimoto Algorithm
//
template<typename eT>
double Shannon<eT>::capacity() {
	uint m = this->C.n_rows;
	uint n = this->C.n_cols;

	Prob<eT> F(m), Px(m), Py(m);
	uniform(Px);

	while(1) {
		// Py = output dist
		Py = Px * this->C;

		// update F
		for(uint i = 0; i < m; i++) {
			double s = 0;
			for(uint j = 0; j < n; j++) {
				double el = this->C.at(i, j);
				s += el > 0 ? el * log(el / Py.at(j)) : 0;		// NOTE: this is e-base log, not 2-base!
			}
			F.at(i) = exp(s);
		}

		// check stop condition
		double d = dot(F, Px);
		double IL = log2(d);
		double IU = log2(max(F));

		if(IU - IL < this->precision)
			return IL;

		// update Px
		Px %= F / d;		// % is element-wise mult
	}
}

template class Shannon<double>;
template class Shannon<float>;
