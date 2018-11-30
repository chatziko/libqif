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

using namespace channel;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ChanTest : public BaseTest<T> {};
template <typename T>
class ChanTestReals : public BaseTest<T> {};

TYPED_TEST_CASE_P(ChanTest);
TYPED_TEST_CASE_P(ChanTestReals);


TYPED_TEST_P(ChanTest, Construct) {
	typedef TypeParam eT;

	const char* s = "1 0 0; 0 1 0";
	Chan<eT> C(s);

	expect_channel( s,  Chan<eT>(s)              ); // char*
	expect_channel( s,  Chan<eT>(std::string(s)) ); // std::string
	expect_channel( s,  Chan<eT>(C)              ); // copy

	// malformed channel
	//
	const char* s2 = "1 2; 3 0.5";
	Chan<eT> C2(s2);

	EXPECT_ANY_THROW( assert_proper(Chan<eT>(s2));              ); // char*
	EXPECT_ANY_THROW( assert_proper(Chan<eT>(std::string(s2))); ); // std::string
	EXPECT_ANY_THROW( assert_proper(Chan<eT>(C2));              ); // Mat
}

TYPED_TEST_P(ChanTest, Identity) {
	typedef TypeParam eT;

	Chan<eT> C;
	C = identity<eT>(0);
	expect_channel(0, 0, C);

	C = identity<eT>(3);
	expect_channel("1 0 0; 0 1 0; 0 0 1", C);
}

TYPED_TEST_P(ChanTest, Randu) {
	typedef TypeParam eT;

	Chan<eT> C(200, 200);
	randu(C);
	expect_channel(200, 200, C);

	C = randu<eT>(5);
	expect_channel(5, 5, C);

	C = randu<eT>(4, 6);
	expect_channel(4, 6, C);
}

TYPED_TEST_P(ChanTest, Factorize) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// non factorizable
	expect_channel(0, 0, factorize(t.id_10, t.noint_10));
	expect_channel(0, 0, factorize(t.id_4,  t.noint_10));

	int n = 4, m = 6;
	Chan<eT>
		B = channel::randu<eT>(n, m),
		A = B * channel::randu<eT>(m, n),
		X1 = factorize(A, B),
		Z1 = B * X1;

	expect_channel(m, n, X1);
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, A, Z1, 1e-4, 0);		// default is subgrad method, with tolerance 1e-4

	// factorize_lp
	//
	// TODO: factorize_lp is unstable under float, it fails half of the time, we should investigate
	if(std::is_same<eT, float>::value) return;

	expect_channel(0, 0, factorize_lp(t.id_10, t.noint_10));
	expect_channel(0, 0, factorize_lp(t.id_4,  t.noint_10));

	Chan<eT>
		X2 = factorize_lp(A, B),
		Z2 = B * X2;

	expect_channel(m, n, X2);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, A, Z2);
}

TYPED_TEST_P(ChanTestReals, FactorizeSubgrad) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// non factorizable
	expect_channel(0, 0, factorize_subgrad(t.id_10, t.noint_10));
	expect_channel(0, 0, factorize_subgrad(t.id_4,  t.noint_10));

	// a "fat" B matrix with cols = 1.5 * rows seems to be a good case where the initial solution of A = B X is
	// not a proper matrix, and needs to be improved by the subgradient method (when cols == rows or cols = 2 * rows it
	// seems that the solution is goot right away, with almost no iterations).
	//
	int n = 10, m = 15;
	Chan<eT>
		B = randu<eT>(n, m),
		A = B * randu<eT>(m, n),
		X = factorize_subgrad(A, B),
		Z = B * X;

	expect_channel(m, n, X);
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, A, Z, 1e-4, 0);

	// the following matrices cause the S matrix of the subgradient method to contain inf, causing X to contain -nan
	//
	A = "0.4405 0.5595;"
		"0.6588 0.3412 ";
	B = "0.9694 0.0062 0.0244;"
		"0.5312 0.1401 0.3287 ";

	X = factorize_subgrad(A, B),
	Z = B * X;

	expect_channel(3, 2, X);
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, A, Z, 1e-4, 0);
}

TYPED_TEST_P(ChanTest, LeftFactorize) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// non factorizable
	expect_channel(0, 0, left_factorize(t.id_10, t.noint_10));
	expect_channel(0, 0, left_factorize(t.id_4,  t.noint_10));

	int n = 4, m = 6;
	Chan<eT>
		B = channel::randu<eT>(n, m),
		A = channel::randu<eT>(m, n) * B,
		X1 = left_factorize(A, B),
		Z1 = X1 * B;

	expect_channel(m, n, X1);
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, A, Z1, 1e-4, 0);		// default is subgrad method, with tolerance 1e-4
}

TYPED_TEST_P(ChanTest, BayesianUpdate) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// the identity channel should produce the real prior in 2 iterations
	Prob<eT> pi;
	EXPECT_EQ(2u, bayesian_update<eT>(t.id_10, t.prand_10 * t.id_10, pi));
	expect_prob(t.prand_10, pi);

	// a non interfering channel should produce the uniform prior in 1 iteration
	pi.clear();
	EXPECT_EQ(1u, bayesian_update<eT>(t.noint_10, t.prand_10 * t.noint_10, pi));
	expect_prob(t.unif_10, pi);

	if(std::is_same<eT, double>::value) {
		// the geometric should produce the real prior in many iterations (and with limited accuracy)
		// Note: for rat this is slow (probably has to do with the huge denominators in the random elements)
		pi.clear();
		auto C = mechanism::geometric<eT>(10);
		bayesian_update<eT>(C, t.prand_10 * C, pi, eT(1e-7));
		expect_prob(t.prand_10, pi, eT(0), eT(1e-4));
	}
}


// run ChanTest for all types, ChanTestReals only for native types
//
REGISTER_TYPED_TEST_CASE_P(ChanTest, Construct, Identity, Randu, Factorize, LeftFactorize, BayesianUpdate);
REGISTER_TYPED_TEST_CASE_P(ChanTestReals, FactorizeSubgrad);

INSTANTIATE_TYPED_TEST_CASE_P(Chan, ChanTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Chan, ChanTestReals, NativeTypes);

