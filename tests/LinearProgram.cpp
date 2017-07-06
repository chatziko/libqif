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
#include "tests_aux.h"
using namespace qif::lp;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class LinearProgramTest : public BaseTest<eT> {};

TYPED_TEST_CASE_P(LinearProgramTest);


TYPED_TEST_P(LinearProgramTest, Optimal) {
	typedef TypeParam eT;
	typedef method_t m_t;
	typedef status_t s_t;
	const bool is_rat = std::is_same<eT, rat>::value;

	// the default acceptance range is too string for linear programs, we need a more permissive mrd
	eT md =	def_max_diff<eT>();
	eT mrd = def_max_rel_diff<float>();		// always use the mrd for floats

	for(m_t method : { m_t::simplex_primal, m_t::simplex_dual, m_t::interior }) {
	for(bool presolve : { false, true }) {
		if(presolve && method == m_t::interior) continue;					// interior method has no presolver
		if(is_rat && (presolve || method != m_t::simplex_primal)) continue;	// rat supports only simplex_primal/no presolve

		LinearProgram<eT> lp;
		lp.method = method;
		lp.glp_presolve = presolve;

		lp.A = format_num<eT>("1 2; 3 1");
		lp.b = format_num<eT>("1 2");
		lp.c = format_num<eT>("0.6 0.5");

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(s_t::optimal, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(46)/100, lp.optimum(), md, mrd);
		expect_mat(format_num<eT>("0.6; 0.2"), lp.x, md, mrd);

		lp.A = format_num<eT>("1 1 0; 0 1 1");
		lp.b = format_num<eT>("1 1");
		lp.c = format_num<eT>("1 2 -1");

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(s_t::optimal, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(2), lp.optimum(), md, mrd);
		expect_mat(format_num<eT>("0; 1; 0"), lp.x, md, mrd);

		lp.maximize = false;
		lp.A = format_num<eT>("3 -4; 1 2; 1 0");
		lp.b = format_num<eT>("12 4 1");
		lp.c = format_num<eT>("3 4");
		lp.sense = "< > >";

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(s_t::optimal, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(9), lp.optimum(), md, mrd);
		expect_mat(format_num<eT>("1; 1.5"), lp.x, md, mrd);

		lp.maximize = false;
		lp.A = format_num<eT>("1 2 2; 2 1 2; 2 2 1");
		lp.b = format_num<eT>("20 20 20");
		lp.c = format_num<eT>("-10 -12 -12");
		lp.sense.set_size(0);

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(s_t::optimal, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(-136), lp.optimum(), md, mrd);
		expect_mat(format_num<eT>("4; 4; 4"), lp.x, md, mrd);
	}}
}

TYPED_TEST_P(LinearProgramTest, Infeasible) {
	typedef TypeParam eT;
	typedef method_t m_t;
	typedef status_t s_t;
	const bool is_rat = std::is_same<eT, rat>::value;

	for(m_t method : { m_t::simplex_primal, m_t::simplex_dual, m_t::interior }) {
	for(bool presolve : { false, true }) {
		if(presolve && method == m_t::interior) continue;					// interior method has no presolver
		if(is_rat && (presolve || method != m_t::simplex_primal)) continue;	// rat supports only simplex_primal/no presolve

		LinearProgram<eT> lp;
		lp.method = method;
		lp.glp_presolve = presolve;

		s_t status = method == m_t::interior
				? s_t::infeasible_or_unbounded	// sometimes we just know that the problem is infeasible OR unbounded
				: s_t::infeasible;

		lp.A = format_num<eT>("1; 1");
		lp.b = format_num<eT>("3 2");
		lp.c = format_num<eT>("1");
		lp.sense = "> <";

		EXPECT_FALSE(lp.solve());
		EXPECT_EQ(status, lp.status);

		lp.A = format_num<eT>("1; -1");
		lp.b = format_num<eT>("3 -2");
		lp.c = format_num<eT>("4");
		lp.sense = "> >";

		EXPECT_FALSE(lp.solve());
		EXPECT_EQ(status, lp.status);
	}}
}

TYPED_TEST_P(LinearProgramTest, Unbounded) {
	typedef TypeParam eT;
	typedef method_t m_t;
	typedef status_t s_t;
	const bool is_rat = std::is_same<eT, rat>::value;

	for(m_t method : { m_t::simplex_primal, m_t::simplex_dual, m_t::interior }) {
	for(bool presolve : { false, true }) {
		if(presolve && method == m_t::interior) continue;					// interior method has no presolver
		if(is_rat && (presolve || method != m_t::simplex_primal)) continue;	// rat supports only simplex_primal/no presolve

		LinearProgram<eT> lp;
		lp.method = method;
		lp.glp_presolve = presolve;

		s_t status = method == m_t::simplex_primal && !presolve
				? s_t::unbounded
				: s_t::infeasible_or_unbounded;	// sometimes we just know that the problem is infeasible OR unbounded

		lp.maximize = false;
		lp.A = format_num<eT>("1");
		lp.b = format_num<eT>("2");
		lp.c = format_num<eT>("-1");
		lp.sense = ">";

		EXPECT_FALSE(lp.solve());
		EXPECT_EQ(status, lp.status);
	}}
}


REGISTER_TYPED_TEST_CASE_P(LinearProgramTest, Optimal, Infeasible, Unbounded);

INSTANTIATE_TYPED_TEST_CASE_P(LinearProgram, LinearProgramTest, AllTypes);

