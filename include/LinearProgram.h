#ifndef _QIF_LinearProgram_h_
#define _QIF_LinearProgram_h_
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


// Solve the linear program
// {min/max} dot(c,x)
// subject to A x {>=|==|<=} b
// x >= 0

template<typename eT>
class LinearProgram {
	public:
		enum class status_t { optimal, infeasible, unbounded, error };
		enum class method_t { simplex_primal, simplex_dual, interior };

		Mat<eT>
			A;				// constraints
		Col<eT>
			x,				// solution
			b,				// constants
			c;				// cost function
		Col<char>
			sense;			// sense of each constraint, can be '<', '=', '>' (<,> really mean <=,>=), default is '<'

		bool maximize = true;
		bool non_negative = 1;
		method_t method = method_t::simplex_primal;
		status_t status;

		LinearProgram() {}
		LinearProgram(const Mat<eT>& A, const Col<eT>& b, const Col<eT>& c) : A(A), b(b), c(c) { check_sizes(); }

		bool solve();
		std::string to_mps();

		inline eT optimum() { return arma::dot(x, c); }
		LinearProgram canonical_form();

		inline char get_sense(uint i) { return i < sense.n_rows ? sense.at(i) : '<'; }		// default sense is <

	protected:
		void check_sizes() { if(A.n_rows != b.n_rows || A.n_cols != c.n_rows) throw "invalid size"; }

		bool glpk();
		bool simplex();
};

#endif
