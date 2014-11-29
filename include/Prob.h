#ifndef _QIF_Prob_h_
#define _QIF_Prob_h_
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

#include "types.h"
#include "Channel.h"

// Note: using EnableIf alias for SFINAE to make things more readable.
// Avoiding SFINAE in function parameters also allows inference of T
// i.e. we can call f(param) instead of f<T>(param)
//

template<typename T, typename = EnableIf<is_Prob<T>>>
inline
T& uniform(T& pi) {
	typedef typename T::elem_type eT;

	pi.fill( eT(1) / eT(pi.n_cols) );
	return pi;
}

template<typename T, typename = EnableIf<is_Prob<T>>>
inline T uniform(uint n) {
	T pi(n);
	return uniform(pi);
}


template<typename T, typename = EnableIf<is_Prob<T>>>
inline
T& dirac(T& pi, uint i = 0) {
	typedef typename T::elem_type eT;

	pi.zeros();
	pi.at(i) = eT(1);
	return pi;
}

template<typename T, typename = EnableIf<is_Prob<T>>>
inline
T dirac(uint n, uint i = 0) {
	T pi(n);
	return dirac(pi, i);
}


template<typename T, typename = EnableIf<is_Prob<T>>>
inline
T& randu(T& pi) {
	typedef typename T::elem_type eT;

	Channel<eT> c(1, pi.n_cols);
	c.randu();
	pi = c.row(0);
	return pi;
}

template<typename T, typename = EnableIf<is_Prob<T>>>
inline
T randu(uint n) {
	T pi(n);
	return randu(pi);
}


template<typename T, typename = EnableIf<is_Prob<T>>>
inline
bool is_proper(const T& pi) {
	typedef typename T::elem_type eT;

	eT sum = 0;
	for(uint j = 0; j < pi.n_cols; j++) {
		// elements should be non-negative
		const eT& elem = pi.at(j);
		if(less_than(elem, eT(0)))
			return false;

		sum += elem;
	}

	// sum should be 1
	if(!equal(sum, eT(1)))
		return false;

	return true;
}

template<typename T, typename = EnableIf<is_Prob<T>>>
inline
void check_proper(const T& pi) {
	if(!is_proper<T>(pi))
		throw 1;
}


#endif
