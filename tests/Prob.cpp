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
#include <string>
#include "gtest/gtest.h"

#include "Prob.h"
#include "aux.h"
#include "tests_aux.h"



// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ProbTest : public BaseTest<T> {};
template <typename T>
class ProbTestReals : public BaseTest<T> {};

TYPED_TEST_CASE_P(ProbTest);
TYPED_TEST_CASE_P(ProbTestReals);		// tests that run only on double/float


TYPED_TEST_P(ProbTest, Construct) {
	typedef TypeParam eT;

	std::string s = format_num<eT>("0.5 0.25 0.25");
	Prob<eT> pi(s);

	expect_prob( s, Prob<eT>(s.c_str()) ); // char*
	expect_prob( s, Prob<eT>(s)         ); // std::string
	expect_prob( s, Prob<eT>(pi)        ); // copy

	// malformed prob
	//
	std::string s2 = format_num<eT>("0.1 0.2 0.3");
	Prob<eT> pi2(s2);

	EXPECT_ANY_THROW( check_proper(Prob<eT>(s2.c_str())); ); // char*
	EXPECT_ANY_THROW( check_proper(Prob<eT>(s2));         ); // std::string
	EXPECT_ANY_THROW( check_proper(Prob<eT>(pi2));        ); // Prob
}

TYPED_TEST_P(ProbTest, Uniform) {
	typedef TypeParam eT;

	Prob<eT> pi(1);
	uniform(pi);
	expect_prob(1, pi);

	pi = uniform<Prob<eT>>(4);
	expect_prob(format_num<eT>("0.25 0.25 0.25 0.25"), pi);
}

TYPED_TEST_P(ProbTest, Randu) {
	typedef TypeParam eT;

	Prob<eT> pi(200);
	randu(pi);
	expect_prob(200, pi);

	pi = randu<Prob<eT>>(5);
	expect_prob(5, pi);
}

TYPED_TEST_P(ProbTest, Dirac) {
	typedef TypeParam eT;

	Prob<eT> pi(4);
	dirac(pi);
	expect_prob("1 0 0 0", pi);

	pi = dirac<Prob<eT>>(4, 2);
	expect_prob("0 0 1 0", pi);
}

TYPED_TEST_P(ProbTest, Total_variation) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, eT(0), total_variation(t.unif_4, t.unif_4));
	EXPECT_PRED2(equal<eT>, eT(0), total_variation(t.dirac_4, t.dirac_4));
	EXPECT_PRED2(equal<eT>, eT(3)/4, total_variation(t.unif_4, t.dirac_4));
	EXPECT_PRED2(equal<eT>, eT(3)/4, total_variation(t.dirac_4, t.unif_4));
}

TYPED_TEST_P(ProbTest, Bounded_entropy_distance) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, eT(0), bounded_entropy_distance(t.unif_4, t.unif_4));
	EXPECT_PRED2(equal<eT>, eT(0), bounded_entropy_distance(t.dirac_4, t.dirac_4));

	EXPECT_PRED2(equal<eT>, eT(1), bounded_entropy_distance(t.unif_4, t.dirac_4));
	EXPECT_PRED2(equal<eT>, eT(1), bounded_entropy_distance(t.dirac_4, t.unif_4));

	EXPECT_PRED2(equal<eT>, eT(9)/eT(14), bounded_entropy_distance(t.unif_4, t.pi5));
	EXPECT_PRED2(equal<eT>, eT(9)/eT(14), bounded_entropy_distance(t.pi5, t.unif_4));
}

TYPED_TEST_P(ProbTestReals, Multiplicative_distance) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, eT(0), multiplicative_distance(t.unif_4, t.unif_4));
	EXPECT_PRED2(equal<eT>, eT(0), multiplicative_distance(t.dirac_4, t.dirac_4));

	EXPECT_EQ(infinity<eT>(), multiplicative_distance(t.unif_4, t.dirac_4));
	EXPECT_EQ(infinity<eT>(), multiplicative_distance(t.dirac_4, t.unif_4));

	EXPECT_PRED2(equal<eT>, std::log(0.7/0.25), multiplicative_distance(t.unif_4, t.pi5));
	EXPECT_PRED2(equal<eT>, std::log(0.7/0.25), multiplicative_distance(t.pi5, t.unif_4));
}


// run the ProbTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_CASE_P(ProbTest, Construct, Uniform, Randu, Dirac, Total_variation, Bounded_entropy_distance);
REGISTER_TYPED_TEST_CASE_P(ProbTestReals, Multiplicative_distance);

INSTANTIATE_TYPED_TEST_CASE_P(Prob, ProbTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Prob, ProbTestReals, NativeTypes);

