/*
This file belongs to the LIBQIF library.
A Quantitative InfoMat<eT>ion Flow C++ Toolkit Library.
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

The canonical_form and simplex methods are adapted from
https://github.com/IainNZ/RationalSimplex.jl
Original code is under the MIT licence.

The MIT License (MIT)
Copyright (c) 2014 Iain Dunning
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "LinearProgram.h"
#include "aux.h"
#include <glpk.h>
#include <string>
#include <cassert>
#include <vector>
#include <exception>

using std::string;
using std::to_string;



template<typename eT>
bool LinearProgram<eT>::solve() {
	check_sizes();

	return glpk();
}

// for rats, we use the simplex() method after transforming to canonical form
//
template<>
bool LinearProgram<rat>::solve() {
	check_sizes();

	if(method != method_t::simplex_primal)
		throw std::runtime_error("not supported");

	LinearProgram<rat> lp = canonical_form();

	bool res = lp.simplex();
	status = lp.status;

	if(res)
		x = lp.x.subvec(0, A.n_cols-1);

	return res;
}


template<typename eT>
bool LinearProgram<eT>::glpk() {
	// create problem
	glp_prob *lp = glp_create_prob();

	glp_set_obj_dir(lp, maximize ? GLP_MAX : GLP_MIN);

	// add variables
	// CAREFULL: all glp indexes are 1-based
	//
	glp_add_cols(lp, A.n_cols);
	for(uint j = 0; j < A.n_cols; j++) {
		if(non_negative)
			glp_set_col_bnds(lp, j+1, GLP_LO, 0.0, 0.0);	// x_j >= 0
		else
			glp_set_col_bnds(lp, j+1, GLP_FR, 0.0, 0.0);	// x_j is free

		glp_set_obj_coef(lp, j+1, c.at(j));					// coefficient in the cost functoin
	}

	// add constraints. glpk uses a "sparse" way of entering the rows, using
	// flat arrays ia, ja, ar. ia[z],ja[z] are the indexes of A for the value to
	// set, and ar[z] = A[ ia[z], ja[z] ]
	// we add entries only for non-zero elements, it's much faster!
	//
	glp_add_rows(lp, A.n_rows);

	int size = 0;
	for(uint i = 0; i < A.n_rows; i++)
		for(uint j = 0; j < A.n_cols; j++)
			if(!equal(A(i, j), eT(0)))
				size++;

	std::vector<int>	ia(size+1),
						ja(size+1);
	std::vector<double> ar(size+1);

	int index = 1;
	for(uint i = 0; i < A.n_rows; i++) {
		char sense_i = sense.n_rows > i ? sense.at(i) : '<';	// default sense is <
		if(sense_i == '<')
			glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, b.at(i));	// row_i dot x <= b(i)
		else if(sense_i == '>')
			glp_set_row_bnds(lp, i+1, GLP_LO, b.at(i), 0.0);	// row_i dot x >= b(i)
		else
			glp_set_row_bnds(lp, i+1, GLP_FX, b.at(i), 0.0);	// row_i dot x >= b(i)

		for(uint j = 0; j < A.n_cols; j++) {
			if(equal(A(i, j), eT(0))) continue;
			ia[index] = i+1;
			ja[index] = j+1;
			ar[index] = A(i, j);
			index++;
		}
	}

	glp_load_matrix(lp, size, &ia[0], &ja[0], &ar[0]);

	// solve
	if(method == method_t::simplex_primal || method == method_t::simplex_dual) {
		glp_smcp opt;
		glp_init_smcp(&opt);
		opt.msg_lev = GLP_OFF;		// no terminal output

		glp_simplex(lp, &opt);

		int glp_status = glp_get_status(lp);
		status =
			glp_status == GLP_OPT	? status_t::optimal :
			glp_status == GLP_NOFEAS? status_t::infeasible :
			glp_status == GLP_UNBND	? status_t::unbounded :
									  status_t::error;

		if(status == status_t::optimal) {
			x.set_size(A.n_cols);
			for(uint j = 0; j < A.n_cols; j++)
				x.at(j) = glp_get_col_prim(lp, j+1);
		}

	} else {
		glp_iptcp opt;
		glp_init_iptcp(&opt);
		opt.msg_lev = GLP_OFF;		// no terminal output

		glp_interior(lp, &opt);

		// NOTE: glpk's interior point returns GLP_NOFEAS also for unbounded problems,
		//       not sure how we can check for boundedness
		//
		int glp_status = glp_ipt_status(lp);
		// std::cout << "interior status: " << (status == GLP_OPT ? "GLP_OPT" : status == GLP_NOFEAS  ? "GLP_NOFEAS" : status == GLP_INFEAS ? "GLP_INFEAS" : status == GLP_UNDEF ? "GLP_UNDEF" : "XXX") << "\n";
		status =
			glp_status == GLP_OPT	? status_t::optimal :
			glp_status == GLP_NOFEAS? status_t::infeasible :
									  status_t::error;

		if(status == status_t::optimal) {
			x.set_size(A.n_cols);
			for(uint j = 0; j < A.n_cols; j++)
				x.at(j) = glp_ipt_col_prim(lp, j+1);
		}
	}

	// clean
	glp_delete_prob(lp);
	glp_free_env();

	return status == status_t::optimal;
}

template<>
bool LinearProgram<rat>::glpk() {
	throw std::runtime_error("not available for rat");
}


// transform the progarm in canonical form:
//        min  dot(c,x)
// subject to  A x == b
//               x >= 0
// b must be >= 0.
//
template<typename eT>
LinearProgram<eT> LinearProgram<eT>::canonical_form() {
	uint m = A.n_rows,
		 n = A.n_cols;

	LinearProgram lp;
	lp.maximize = false;		// minimization
	lp.non_negative = 1;		// x_i >= 0
	lp.sense.set_size(m);		// all constraints
	lp.sense.fill('=');			// are equalities

	// Count number of auxiliaries we will need
	uint extra = 0;
	for(uint i = 0; i < m; i++)
		if(get_sense(i) != '=')
			extra += 1;

	lp.A = arma::zeros<Mat<eT>>(m, n + extra);
	lp.A.cols(0, n-1) = A;

	lp.c = arma::zeros<Col<eT>>(n + extra);
	lp.c.rows(0, n-1) = c;

	if(maximize)		// canonical form program is minimization
		lp.c *= eT(-1);

	// Add the auxiliaries
	uint offset = 0;
	for(uint i = 0; i < m; i++) {
		if(get_sense(i) != '=') {
			lp.A.at(i, n + offset) = eT(get_sense(i) == '<' ? 1 : -1);
			offset++;
		}
	}

	// Make sure right-hand-side is non-negative
	lp.b = b;
	for(uint i = 0; i < m; i++) {
		if(b.at(i) < eT(0)) {
			lp.A.row(i) *= eT(-1);
			lp.b.at(i) *= eT(-1);
		}
	}

	return lp;
}

// simplex
// Solve the linear program in canonical form
//        min  dot(c,x)
// subject to  A x == b
//               x >= 0
// b must be >= 0.
// 
// This is mainly to be used with rats
//
// The algorithm is the two-phase primal revised simplex method.
// In the first phase auxiliaries are created which we eliminate
// until we have a basis consisting solely of actual variables.
// This is pretty much the "textbook algorithm", and shouldn't
// be used for anything that matters. It doesn't exploit sparsity
// at all. You could use it with floating points but it wouldn't 
// work for anything except the most simple problem due to accumulated
// errors and the comparisons with zero.
//
template<typename eT>
bool LinearProgram<eT>::simplex() {
	using arma::zeros;
	using arma::ones;
	using arma::eye;
	using arma::umat;

	uint m = A.n_rows,
		 n = A.n_cols;

	assert(!maximize);
	for(uint i = 0; i < m; i++)
		assert(!less_than(b.at(i), eT(0)));

	Mat<char> is_basic	= zeros<Mat<char>>(n + m);
	umat basic			= zeros<umat>(m);				// indices of current basis
	Mat<eT> Binv		= eye<Mat<eT>>(m, m);			// inverse of basis matrix
	Row<eT> cB			= ones<Row<eT>>(m);				// costs of basic variables
	x					= zeros<Col<eT>>(n + m);		// current solution

	// Intialize phase 1
	for(uint i = 0; i < m; i++) {
		basic(i) = i + n;
		is_basic(i + n) = 1;
		x(i + n) = b(i);
	}
	bool phase_one = true;

	// Begin simplex iterations
	while(true) {
		// Calculate dual solution...
		Row<eT> pi_T = cB * Binv;

		// ... and thus the reduced costs
		int entering = -1;
		for(uint j = 0; j < n; j++) {
			if(is_basic(j)) continue;
			eT rc = (phase_one ? eT(0) : c(j)) - dot(pi_T, A.col(j));
			if(less_than(rc, eT(0))) {
				entering = j;
				break;
			}
		}

		// If we couldn't find a variable with a negative reduced cost, 
		// we terminate this phase because we are at optimality for this
		// phase - not necessarily optimal for the actual problem.
		if(entering == -1) {
			if(phase_one) {
				phase_one = false;
				// Check objective - if 0, we are OK
				for(uint j = n; j < n + m; j++) {
					if(less_than(eT(0), x(j))) {
						// It couldn't reduce objective to 0 which is equivalent
						// to saying a feasible basis with no artificials could
						// not be found
						status = status_t::infeasible;
						// std::cout << "-- infeasible\n";
						goto EXIT;
					}
				}
				// Start again in phase 2 with our nice feasible basis
				for(uint i = 0; i < m; i++) {
					cB(i) = basic(i) > n ? eT(0) : c(basic(i));
				}
				continue;
			} else {
				status = status_t::optimal;
				// std::cout << "-- optimal\n";
				goto EXIT;
			}
		}

		// Calculate how the solution will change when our new
		// variable enters the basis and increases from 0
		Col<eT> BinvAs = Binv * A.col(entering);

		// Perform a "ratio test" on each variable to determine
		// which will reach 0 first
		int leaving = -1;
		eT min_ratio = eT(0);
		for(uint i = 0; i < m; i++) {
			if(less_than(eT(0), BinvAs(i))) {
				eT ratio = x(basic(i)) / BinvAs(i);
				if(less_than(ratio, min_ratio) || leaving == -1) {
					min_ratio = ratio;
					leaving = i;
				}
			}
		}

		// If no variable will leave basis, then we have an 
		// unbounded problem.
		if(leaving == -1) {
			status = status_t::unbounded;
			// std::cout << "--unbounded\n";
			goto EXIT;
		}

		// Now we update solution...
		for(uint i = 0; i < m; i++) {
			x(basic(i)) -= min_ratio * BinvAs(i);
		}
		x(entering) = min_ratio;

		// ... and the basis inverse...
		// Our tableau is ( Binv b | Binv | BinvAs )
		// and we doing a pivot on the leaving row of BinvAs
		eT pivot_value = BinvAs(leaving);
		for(uint i = 0; i < m; i++) {  // all rows except leaving row
			if(i == static_cast<uint>(leaving)) continue;
			eT factor = BinvAs(i) / pivot_value;
			for(uint j = 0; j < m; j++)
				Binv(i, j) -= factor * Binv(leaving, j);
		}
		for(uint j = 0; j < m; j++)
			Binv(leaving, j) /= pivot_value;

		// ... and variable status flags
		is_basic(basic(leaving)) = 0;
		is_basic(entering) = 0;
		cB(leaving) = phase_one ? eT(0) : c(entering);
		basic(leaving) = entering;
	}

EXIT:
	x = x.subvec(0, n-1);		// the solution are the first n vars

	return status == status_t::optimal;
}


template<typename eT>
string LinearProgram<eT>::to_mps() {
	string s;
	s = "NAME PROG\n";

	// rows
	s += "ROWS\n";
	s += " N  OBJ\n";	// objective function
	for(uint i = 0; i < A.n_rows; i++) {
		char sense_i = sense.n_rows > i ? sense.at(i) : '<';	// default sense is <
		string s_sense = sense_i == '<' ? "L" :
						 sense_i == '>' ? "G" :
							 			  "E";
		s += " " + s_sense + " ROW" + to_string(i+1) + "\n";
	}

	// columns
	s += "COLUMNS\n";
	for(uint j = 0; j < A.n_cols; j++) {
		s += " X" + to_string(j+1) + " OBJ " + to_string(c.at(j)) + "\n";

		for(uint i = 0; i < A.n_rows; i++)
			s += " X" + to_string(j+1) + " ROW" + to_string(i+1) + " " + to_string(A.at(i, j)) + "\n";
	}

	// RHS
	s += "RHS\n";
	for(uint i = 0; i < A.n_rows; i++) {
		s += " RHS ROW" + to_string(i+1) + " " + to_string(b.at(i)) + "\n";
	}

	if(!non_negative) {
		// No BOUNDS assumes >= 0
		s += "\nBOUNDS\n";
		for(uint j = 0; j < A.n_cols; j++)
			s += " FR BND X" + to_string(j+1) + "\n";
	}

	s += "ENDATA\n";

	return s;
}

template<>
string LinearProgram<rat>::to_mps() {
	// TODO: make to_mps work for rat
	throw std::runtime_error("not supported");
}

template class LinearProgram<double>;
template class LinearProgram<float>;
template class LinearProgram<rat>;


