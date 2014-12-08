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

#include "Chan.h"
#include "tests_aux.h"


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ChanTest : public BaseTest<T> {};

TYPED_TEST_CASE_P(ChanTest);


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

	EXPECT_ANY_THROW( check_proper(Chan<eT>(s2));              ); // char*
	EXPECT_ANY_THROW( check_proper(Chan<eT>(std::string(s2))); ); // std::string
	EXPECT_ANY_THROW( check_proper(Chan<eT>(C2));              ); // Mat
}

TYPED_TEST_P(ChanTest, Identity) {
	typedef TypeParam eT;

	Chan<eT> C;
	C = identity<Chan<eT>>(0);
	expect_channel(0, 0, C);

	C = identity<Chan<eT>>(3);
	expect_channel("1 0 0; 0 1 0; 0 0 1", C);
}

TYPED_TEST_P(ChanTest, Randu) {
	typedef TypeParam eT;

	Chan<eT> C(200, 200);
	randu(C);
	expect_channel(200, 200, C);

	C = randu<Chan<eT>>(5);
	expect_channel(5, 5, C);

	C = randu<Chan<eT>>(4, 6);
	expect_channel(4, 6, C);
}

TYPED_TEST_P(ChanTest, Factorize) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// non factorizable
	expect_channel(0, 0, factorize(t.id_10, t.noint_10));
	expect_channel(0, 0, factorize(t.id_4,  t.noint_10));

	// TODO: Factorize is unstable under float, it fails half of the time, we should investigate
	if(std::is_same<eT, float>::value) return;

	int n = 5, m = 4;
	Chan<eT>
		B = randu<Chan<eT>>(n, n),
		A = B * randu<Chan<eT>>(n, m),
		X = factorize(A, B),
		Z = B * X;

	expect_channel(n, m, X);
	expect_channel(A, Z);
}



// run the ChanTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_CASE_P(ChanTest, Construct, Identity, Randu, Factorize);

INSTANTIATE_TYPED_TEST_CASE_P(Chan, ChanTest, AllTypes);

