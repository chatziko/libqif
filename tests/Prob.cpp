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


TYPED_TEST_P(ProbTest, constructors) {
	typedef TypeParam eT;

	const char* s = std::is_same<eT, rat>::value
		? "1/2 1/4 1/4"
		: "0.5 0.25 0.25";
	Row<eT> m(s);
	Prob<eT> p(s);

	expect_prob( s,  Prob<eT>(s)              ); // char*
	expect_prob( s,  Prob<eT>(std::string(s)) ); // std::string
	expect_prob( s,  Prob<eT>(m)              ); // Row
	expect_prob( s,  Prob<eT>(p)              ); // copy
	expect_prob( s,  Prob<eT>(Prob<eT>(s))    ); // move (Note: most likely the compiler is removing the move call completely, see http://en.cppreference.com/w/cpp/language/copy_elision)
	expect_prob( s,  Prob<eT>(Row<eT>(s))     ); // Row, move semantics

	// malformed prob
	// Note: "cout <<" is to avoid the compiler removing the code as unused!
	//
	const char* s2 = std::is_same<eT, rat>::value
		? "1/2 1/3 1/8"
		: "0.1 0.1 0.3";
	Row<eT> m2(s2);

	EXPECT_ANY_THROW( cout << Prob<eT>(s2);              ); // char*
	EXPECT_ANY_THROW( cout << Prob<eT>(std::string(s2)); ); // std::string
	EXPECT_ANY_THROW( cout << Prob<eT>(m2);              ); // Row
	EXPECT_ANY_THROW( cout << Prob<eT>(Row<eT>(s2));     ); // Row, move semantics
}

TYPED_TEST_P(ProbTest, uniform) {
	typedef TypeParam eT;

	Prob<eT> p;
	p.uniform(1);
	expect_prob(1, p);

	p.uniform(4);
	const char* s = std::is_same<eT, rat>::value
		? "1/4 1/4 1/4 1/4"
		: "0.25 0.25 0.25 0.25";
	expect_prob(s, p);
}

TYPED_TEST_P(ProbTest, randu) {
	typedef TypeParam eT;

	Prob<eT> p(200);
	p.randu();
	expect_prob(200, p);

	p.randu(5);
	expect_prob(5, p);
}

TYPED_TEST_P(ProbTest, dirac) {
	typedef TypeParam eT;

	Prob<eT> p(4);
	p.dirac(0);
	const char* s = std::is_same<eT, rat>::value
		? "1/1 0/1 0/1 0/1"
		: "1 0 0 0";
	expect_prob(s, p);

	p.dirac(2);
	s = std::is_same<eT, rat>::value
		? "0/1 0/1 1/1 0/1"
		: "0 0 1 0";
	expect_prob(s, p);
}



// run the ProbTest test-case for double, float, rat
//
REGISTER_TYPED_TEST_CASE_P(ProbTest, constructors, uniform, randu, dirac);

typedef ::testing::Types<double, float, rat> ProbTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Prob, ProbTest, ProbTypes);

