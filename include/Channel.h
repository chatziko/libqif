#ifndef _QIF_Channel_h_
#define _QIF_Channel_h_
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
#include "Prob.h"

// Note: using EnableIf alias for SFINAE to make things more readable.
// We use the last method from http://loungecpp.wikidot.com/tips-and-tricks:enable-if-for-c-11
// that uses parameter packs, i.e. EnableIf<cond>...
// This allows to have functions with exactly the same signature, differing only in the EnableIf/DisableIf conditions
//

template<typename T, EnableIf<is_Channel<T>>...>
inline
T& identity(T& C) {
	if(!C.is_square()) throw "not square";
	C.eye();
	return C;
}

template<typename T, EnableIf<is_Channel<T>>...>
inline
T identity(uint n) {
	T C(n, n);
	return identity(C);
}


template<typename T, EnableIf<is_Channel<T>>...>
inline
T& no_interference(T& C) {
	typedef typename T::elem_type eT;

	C.zeros();
	C.col(0).fill(eT(1));
	return C;
}

template<typename T, EnableIf<is_Channel<T>>...>
inline
T no_interference(uint n) {
	T C(n, 1);
	return no_interference(C);
}


template<typename T, EnableIf<is_Channel<T>>...>
inline
T& randu(T& C) {
	typedef typename T::elem_type eT;

	for(uint i = 0; i < C.n_rows; i++)
		C.row(i) = randu<Prob<eT>>(C.n_cols);

	return C;
}

template<typename T, EnableIf<is_Channel<T>>...>
inline
T randu(uint n) {
	T C(n, n);
	return randu(C);
}

template<typename T, EnableIf<is_Channel<T>>...>
inline
T randu(uint n, uint m) {
	T C(n, m);
	return randu(C);
}


template<typename T, EnableIf<is_Channel<T>>...>
inline
bool is_proper(const T& C) {
	typedef typename T::elem_type eT;

	for(uint i = 0; i < C.n_rows; i++)
		if(!is_proper<Prob<eT>>(C.row(i)))
			return false;

	return true;
}


template<typename T, EnableIf<is_Channel<T>>...>
inline
void check_proper(const T& C) {
	if(!is_proper<T>(C))
		throw 1;
}


#endif
