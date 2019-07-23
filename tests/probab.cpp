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

using namespace probab;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ProbTest : public BaseTest<T> {};

TYPED_TEST_CASE_P(ProbTest);


TYPED_TEST_P(ProbTest, Construct) {
	typedef TypeParam eT;

	std::string s = format_num<eT>("0.5 0.25 0.25");
	Prob<eT> pi = { eT(1)/2, eT(1)/4, eT(1)/4 };

	EXPECT_PRED_FORMAT2(prob_equal2<eT>, Prob<eT>(s.c_str()),  pi); // char*
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, Prob<eT>(s),          pi); // std::string
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, Prob<eT>(pi),         pi); // copy

	// malformed prob
	//
	std::string s2 = format_num<eT>("0.1 0.2 0.3");
	Prob<eT> pi2 = { eT(1)/10, eT(2)/10, eT(3)/10 };

	EXPECT_ANY_THROW( assert_proper(Prob<eT>(s2.c_str())); ); // char*
	EXPECT_ANY_THROW( assert_proper(Prob<eT>(s2));         ); // std::string
	EXPECT_ANY_THROW( assert_proper(Prob<eT>(pi2));        ); // Prob
}

TYPED_TEST_P(ProbTest, Uniform) {
	typedef TypeParam eT;

	Prob<eT> pi(1);
	uniform(pi);
	EXPECT_PRED_FORMAT2(prob_is_proper_size2<eT>, pi, 1);

	pi = uniform<eT>(4);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, format_num<eT>("0.25 0.25 0.25 0.25"));
}

TYPED_TEST_P(ProbTest, Randu) {
	typedef TypeParam eT;

	Prob<eT> pi(200);
	randu(pi);
	EXPECT_PRED_FORMAT2(prob_is_proper_size2<eT>, pi, 200);

	pi = probab::randu<eT>(5);
	EXPECT_PRED_FORMAT2(prob_is_proper_size2<eT>, pi, 5);
}

TYPED_TEST_P(ProbTest, Dirac) {
	typedef TypeParam eT;

	Prob<eT> pi(4);
	dirac(pi);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, "1 0 0 0");

	pi = probab::dirac<eT>(4, 2);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, "0 0 1 0");
}


// run the ProbTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_CASE_P(ProbTest, Construct, Uniform, Randu, Dirac);

INSTANTIATE_TYPED_TEST_CASE_P(Prob, ProbTest, AllTypes);

