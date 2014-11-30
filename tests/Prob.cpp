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
#include "Prob.h"
#include "tests_aux.h"
#include "gtest/gtest.h"
#include <string>
#include <type_traits>

using namespace std;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ProbTest : public ::testing::Test {};

TYPED_TEST_CASE_P(ProbTest);


TYPED_TEST_P(ProbTest, Construct) {
	typedef TypeParam eT;

	const char* s = "0.5 0.25 0.25";
	Prob<eT> pi("0.5 0.25 0.25");

	expect_prob( s,  Prob<eT>(s)              ); // char*
	expect_prob( s,  Prob<eT>(std::string(s)) ); // std::string
	expect_prob( s,  Prob<eT>(pi)             ); // copy

	// malformed prob
	// Note: "cout <<" is to avoid the compiler removing the code as unused!
	//
	const char* s2 = "0.1 0.2 0.3";
	Prob<eT> pi2(s2);

	EXPECT_ANY_THROW( check_proper(Prob<eT>(s2));              ); // char*
	EXPECT_ANY_THROW( check_proper(Prob<eT>(std::string(s2))); ); // std::string
	EXPECT_ANY_THROW( check_proper(Prob<eT>(pi2));              ); // Prob
}

TYPED_TEST_P(ProbTest, Uniform) {
	typedef TypeParam eT;

	Prob<eT> pi(1);
	uniform(pi);
	expect_prob(1, pi);

	pi = uniform<Prob<eT>>(4);
	expect_prob("0.25 0.25 0.25 0.25", pi);
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

	Prob<eT> pi1(4), pi2(4);
	dirac(pi1);
	uniform(pi2);

	EXPECT_PRED2(equal<eT>, eT(0),				total_variation(pi1, pi1));
	EXPECT_PRED2(equal<eT>, eT(0),				total_variation(pi2, pi2));
	EXPECT_PRED2(equal<eT>, (pi1(0)-pi2(0))/2,	total_variation(pi1, pi2));
	EXPECT_PRED2(equal<eT>, (pi1(0)-pi2(0))/2,	total_variation(pi2, pi1));
}

TYPED_TEST_P(ProbTest, Bounded_entropy_distance) {
	typedef TypeParam eT;

	Prob<eT> pi1(4), pi2(4);
	dirac(pi1);
	uniform(pi2);
	Prob<eT> pi3("0.1 0.1 0.1 0.7");

	EXPECT_PRED2(equal<eT>, eT(0),	bounded_entropy_distance(pi1, pi1));
	EXPECT_PRED2(equal<eT>, eT(0),	bounded_entropy_distance(pi2, pi2));

	EXPECT_PRED2(equal<eT>, eT(1),	bounded_entropy_distance(pi1, pi2));
	EXPECT_PRED2(equal<eT>, eT(1),	bounded_entropy_distance(pi2, pi1));

	EXPECT_PRED2(equal<eT>, 0.642857143,	bounded_entropy_distance(pi2, pi3));
	EXPECT_PRED2(equal<eT>, 0.642857143,	bounded_entropy_distance(pi3, pi2));
}



// run the ProbTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_CASE_P(ProbTest, Construct, Uniform, Randu, Dirac, Total_variation, Bounded_entropy_distance);

typedef ::testing::Types<double, float, urat> ProbTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Prob, ProbTest, ProbTypes);

