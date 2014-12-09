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


template<typename T, EnableIf<is_Prob<T>>...>
inline
T normalise_prob(const T& pi) {
	return pi / sum(pi);
}


template<typename T, EnableIf<is_Prob<T>>...>
inline
T& randu(T& pi) {
	pi.randu();
	pi = normalise_prob(pi);
	return pi;
}

// randu for rationals
// Problem: how to select uniformly a rational
// We could generate a double and convert to rational, however
// having too many rationals with different denominators is problematic,
// summing them creates an overflow and is_proper returns false.
// So we simply take a common denominator (4096) and generate nominator
// uniformly in [0,den]
//
template<>
inline
rprob& randu<rprob>(rprob& pi) {
	const int den = 4096;
	Mat<int> m(1, 1);

	for(auto& e : pi) {
		m = arma::randi<Mat<int>>(1, 1, arma::distr_param(0, den));		// use whatever random number generator armadillo is using
		e = rat(m.at(0,0));
		e /= den;
	}
	pi = normalise_prob(pi);

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
		throw 1;
}



template<typename T, EnableIf<is_Prob<T>>...>
inline
typename T::elem_type total_variation(const T& x, const T& y) {
	if(x.n_cols != y.n_cols) throw "size mismatch";

	return arma::max(arma::abs(x - y)) / 2;
}


// max_i | ln x[i] - ln y[i] |
//
template<typename T, EnableIf<is_Prob<T>>...>
inline
typename T::elem_type multiplicative_distance(const T& x, const T& y) {
	typedef typename T::elem_type eT;

	if(x.n_cols != y.n_cols) throw "size mismatch";

	eT res = eT(0);
	for(uint i = 0; i < x.n_cols; i++) {
		eT a = x(i),
		   b = y(i);
		bool a_is_zero = equal(a, eT(0)),
			 b_is_zero = equal(b, eT(0));

		if(a_is_zero && b_is_zero)			// both are the same, no diff
			continue;
		else if(a_is_zero || b_is_zero)		// only one is zero, diff is infty
			return infinity<eT>();
		else {
			eT diff = std::abs(std::log(a) - std::log(b));
			if(less_than(res, diff))
				res = diff;
		}
	}
	return res;
}


template<typename T, EnableIf<is_Prob<T>>...>
inline
typename T::elem_type bounded_entropy_distance(const T& x, const T& y) {
	typedef typename T::elem_type eT;

	if(x.n_cols != y.n_cols) throw "size mismatch";

	eT res = eT(0);
	for(uint i = 0; i < x.n_cols; i++) {
		eT m = std::max(x.at(i), y.at(i));
		if(equal(m, eT(0)))
			continue;

		eT d = abs(x.at(i) - y.at(i)) / m;
		if(d > res)
			res = d;
	}
	return res;
}


#endif
