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
#include "gtest/gtest.h"
#include "Shannon.h"
#include "aux.h"
#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ShannonTest : public ::testing::Test {};

TYPED_TEST_CASE_P(ShannonTest);



TYPED_TEST_P(ShannonTest, Vulnerability) {
	typedef TypeParam eT;

	Prob<eT> pi("0.5 0.5");
	Shannon<eT> shan(Chan<eT>("1 0; 0 1"));

	ASSERT_ANY_THROW( shan.vulnerability(pi); );
	ASSERT_ANY_THROW( shan.cond_vulnerability(pi); );
}

TYPED_TEST_P(ShannonTest, Entropy) {
	typedef TypeParam eT;

	Prob<eT> pi;
	Shannon<eT> shan;

	pi = uniform<Prob<eT>>(2);
	EXPECT_PRED2(equal<eT>, 1, shan.entropy(pi));

	pi = uniform<Prob<eT>>(10);
	EXPECT_PRED2(equal<eT>, log2(10), shan.entropy(pi));

	pi = Prob<eT>("1 0 0 0");
	EXPECT_PRED2(equal<eT>, 0, shan.entropy(pi));

	pi = Prob<eT>("0.2 0.8");
	EXPECT_PRED2(equal<eT>, 0.721928094887362, shan.entropy(pi));
}

TYPED_TEST_P(ShannonTest, Cond_entropy) {
	typedef TypeParam eT;

	Shannon<eT> shan;

	shan.C = identity<Chan<eT>>(2);
	EXPECT_PRED2(equal<eT>, 0, shan.cond_entropy(uniform<Prob<eT>>(2)));
	EXPECT_PRED2(equal<eT>, 0, shan.cond_entropy(dirac<Prob<eT>>(2)));
	EXPECT_PRED2(equal<eT>, 0, shan.cond_entropy(Prob<eT>("0.2 0.8")));

	shan.C = identity<Chan<eT>>(10);
	EXPECT_PRED2(equal<eT>, 0, shan.cond_entropy(uniform<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 0, shan.cond_entropy(dirac<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 0, shan.cond_entropy(Prob<eT>("0.2 0.8 0 0 0 0 0 0 0 0")));

	no_interference(shan.C);
	EXPECT_PRED2(equal<eT>, log2(10), shan.cond_entropy(uniform<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 0, shan.cond_entropy(dirac<Prob<eT>>(10)));

	Prob<eT> pi = randu<Prob<eT>>(10);
	EXPECT_PRED2(equal<eT>, shan.entropy(pi), shan.cond_entropy(pi));

	shan.C = Chan<eT>("0.8 0.2; 0.3 0.7");
	pi = "0.25 0.75";
	EXPECT_PRED2(equal<eT>, 0.669020059980807, shan.cond_entropy(pi));

	shan.C = identity<Chan<eT>>(10);
	ASSERT_ANY_THROW(shan.cond_entropy(uniform<Prob<eT>>(2)););
}

TYPED_TEST_P(ShannonTest, Capacity) {
	typedef TypeParam eT;

	Shannon<eT> shan;

	shan.C = identity<Chan<eT>>(2);
	EXPECT_PRED2(equal<eT>, 1, shan.capacity());

	shan.C = identity<Chan<eT>>(10);
	EXPECT_PRED2(equal<eT>, log2(10), shan.capacity());

	shan.C = no_interference<Chan<eT>>(10);
	EXPECT_PRED2(equal<eT>, 0, shan.capacity());

	shan.C = Chan<eT>("0.8 0.2; 0.3 0.7");
	EXPECT_PRED2(equal<eT>, 0.19123721482206, shan.capacity());

	// symmetric
	shan.C = Chan<eT>(
		".3 .2 .5;"
		".5 .3 .2;"
		".2 .5 .3;"
	);
	double cap = log2(shan.C.n_cols) - shan.entropy(shan.C.row(0));
	EXPECT_PRED2(equal<eT>, cap, shan.capacity());

	// weakly symmetric
	shan.C = Chan<eT>(
		"0.333333333 0.166666667 0.5;"
		"0.333333333 0.5         0.166666667;"
	);
	cap = log2(shan.C.n_cols) - shan.entropy(shan.C.row(0));
	EXPECT_PRED2(equal<eT>, cap, shan.capacity());
}


// run the ChanTest test-case for double, float
//
REGISTER_TYPED_TEST_CASE_P(ShannonTest, Vulnerability, Entropy, Cond_entropy, Capacity);

INSTANTIATE_TYPED_TEST_CASE_P(Shannon, ShannonTest, NativeTypes);

