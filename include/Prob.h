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

#include <exception>
#include "types.h"
#include "aux.h"

// Note: using EnableIf alias for SFINAE to make things more readable.
// We use the last method from http://loungecpp.wikidot.com/tips-and-tricks:enable-if-for-c-11
// that uses parameter packs, i.e. EnableIf<cond>...
// This allows to have functions with exactly the same signature, differing only in the EnableIf conditions
//

template<typename T, EnableIf<is_Prob<T>>...>
inline
T& uniform(T& pi) {
	typedef typename T::elem_type eT;

	pi.fill( eT(1) / eT(pi.n_cols) );
	return pi;
}

template<typename T, EnableIf<is_Prob<T>>...>
inline T uniform(uint n) {
	T pi(n);
	return uniform(pi);
}


template<typename T, EnableIf<is_Prob<T>>...>
inline
T& dirac(T& pi, uint i = 0) {
	typedef typename T::elem_type eT;

	pi.zeros();
	pi.at(i) = eT(1);
	return pi;
}

template<typename T, EnableIf<is_Prob<T>>...>
inline
T dirac(uint n, uint i = 0) {
	T pi(n);
	return dirac(pi, i);
}


// Generate a distribution of n elements uniformly (among all the elements of the n-1 simplex)
//
// Algorithm: Generate n-1 numbers uniformly in [0,1], add 0 and 1 to the list,
// sort, and take the difference between consecutive elements.
// This algorithm has advantages over the simple "normalize n uniform elements" tehcnique:
//
// 1. The normalizing technique is not uniform! See the url below.
//
// 2. The algorithm involves no divisions, and the resulting sum is much closer to exactly 1.0
//    Using this algorithm with float the Kantorovich tests over random dists always pass, while with
//    the normalizing algorithm there are instabilities.
//
// 3. We avoid to sum all elements, which creates a big denominator under rats
//
// Note: if the n logn complexity is a problem, there's a linear algorithm involving logs in the following url:
// http://stats.stackexchange.com/questions/14059/generate-uniformly-distributed-weights-that-sum-to-unity
//
template<typename T, EnableIf<is_Prob<T>>...>
inline
T& randu(T& pi) {
	typedef typename T::elem_type eT;

	pi.randu();
	pi(pi.n_cols-1) = eT(1);		// add 1 to the list. We don't really need to add 0
	pi = arma::sort(pi);

	for(uint i = pi.n_cols-1; i > 0; i--)
		pi(i) -= pi(i-1);

	// naif normalize algorithm
	//pi.randu();
	//pi /= arma::accu(pi);

	return pi;
}

template<typename T, EnableIf<is_Prob<T>>...>
inline
T randu(uint n) {
	T pi(n);
	return randu(pi);
}


template<typename T, EnableIf<is_Prob<T>>...>
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

template<typename T, EnableIf<is_Prob<T>>...>
inline
void check_proper(const T& pi) {
	if(!is_proper<T>(pi))
		throw std::runtime_error("not a proper dist");
}


#endif
