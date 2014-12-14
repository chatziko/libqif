#ifndef _QIF_Chan_h_
#define _QIF_Chan_h_
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
#include "LinearProgram.h"

// Note: using EnableIf alias for SFINAE to make things more readable.
// We use the last method from http://loungecpp.wikidot.com/tips-and-tricks:enable-if-for-c-11
// that uses parameter packs, i.e. EnableIf<cond>...
// This allows to have functions with exactly the same signature, differing only in the EnableIf/DisableIf conditions
//

template<typename T, EnableIf<is_Chan<T>>...>
inline
T& identity(T& C) {
	if(!C.is_square()) throw std::runtime_error("not square");
	C.eye();
	return C;
}

template<typename T, EnableIf<is_Chan<T>>...>
inline
T identity(uint n) {
	T C(n, n);
	return identity(C);
}


template<typename T, EnableIf<is_Chan<T>>...>
inline
T& no_interference(T& C) {
	typedef typename T::elem_type eT;

	C.zeros();
	C.col(0).fill(eT(1));
	return C;
}

template<typename T, EnableIf<is_Chan<T>>...>
inline
T no_interference(uint n) {
	T C(n, 1);
	return no_interference(C);
}


template<typename T, EnableIf<is_Chan<T>>...>
inline
T& randu(T& C) {
	typedef typename T::elem_type eT;

	for(uint i = 0; i < C.n_rows; i++)
		C.row(i) = randu<Prob<eT>>(C.n_cols);

	return C;
}

template<typename T, EnableIf<is_Chan<T>>...>
inline
T randu(uint n) {
	T C(n, n);
	return randu(C);
}

template<typename T, EnableIf<is_Chan<T>>...>
inline
T randu(uint n, uint m) {
	T C(n, m);
	return randu(C);
}


template<typename T, EnableIf<is_Chan<T>>...>
inline
bool is_proper(const T& C) {
	typedef typename T::elem_type eT;

	for(uint i = 0; i < C.n_rows; i++)
		if(!is_proper<Prob<eT>>(C.row(i)))
			return false;

	return true;
}


template<typename T, EnableIf<is_Chan<T>>...>
inline
void check_proper(const T& C) {
	if(!is_proper<T>(C))
		throw std::runtime_error("not a proper matrix");
}


template<typename T, EnableIf<is_Chan<T>>...>
inline bool chan_equal(const T& A, const T& B) {
	if(A.n_rows != B.n_rows || A.n_cols != B.n_cols)
		return false;

	for(uint i = 0; i < A.n_rows; i++)
		for(uint j = 0; j < A.n_cols; j++)
			if(!equal(A.at(i, j), B.at(i, j)))
				return false;

	return true;
}


// Returns a channel X such that A = B X
//
template<typename T, EnableIf<is_Chan<T>>...>
inline
T factorize(const T& A, const T& B) {
	typedef typename T::elem_type eT;

	// A: M x N
	// B: M x R
	// X: R x N   unknowns
	//
	uint M = A.n_rows,
		 N = A.n_cols,
		 R = B.n_cols,
		 n_vars = R * N,		// one var for each element of X
		 n_cons = M * N + R;	// one constraint for each element of A (A[i,j] = dot(B[i,:], X[: j]), plus one constraint for row of X (sum = 1)

	if(B.n_rows != M)
		return T();

	LinearProgram<eT> lp;
	lp.A = arma::zeros<T>(n_cons, n_vars);
	lp.b.set_size(n_cons);
	lp.c = arma::zeros<T>(n_vars);			// we don't really care to optimize, so cost function = 0
	lp.sense.set_size(n_cons);
	lp.sense.fill('=');

	// Build equations for A = B X
	// We have R x N variables, that will be unfolded in a vector.
	// The varialbe X[r,n] will have variable number rN+n.
	// For each element m,n of A we have an equation A[i,j] = dot(B[i,:], X[: j])
	//
	for(uint m = 0; m < M; m++) {
		for(uint n = 0; n < N; n++) {
			uint row = m*N+n;
			lp.b(row) = A(m, n);		// sum of row is A[m,n]

			for(uint r = 0; r < R; r++)
				lp.A(row, r*N+n) = B(m, r);		// coeff B[m,r] for variable X[r,n]
		}
	}

	// equalities for summing up to 1
	//
	for(uint r = 0; r < R; r++) {
		uint row = M*N+r;
		lp.b(row) = eT(1);		// sum of row = 1

		for(uint n = 0; n < N; n++)
			lp.A(row, r*N+n) = eT(1);	// coeff 1 for variable X[r,n]
	}

	// solve program
	//
	if(!lp.solve())
		return T();

	// reconstrict channel from solution
	//
	Chan<eT> X(R, N);
	for(uint r = 0; r < R; r++)
		for(uint n = 0; n < N; n++)
			X(r, n) = lp.x(r*N+n);

	return X;
}


#endif
