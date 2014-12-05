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
#include "MinEntropy.h"
#include "aux.h"
#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MinEntropyTest : public BaseTest<eT> {};
template <typename eT>
class MinEntropyTestReals : public BaseTest<eT> {};

TYPED_TEST_CASE_P(MinEntropyTest);
TYPED_TEST_CASE_P(MinEntropyTestReals);		// tests that run only on double/float


TYPED_TEST_P(MinEntropyTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, eT(5)/10, MinEntropy<eT>().vulnerability(t.unif_2));
	EXPECT_PRED2(equal<eT>, eT(1)/10, MinEntropy<eT>().vulnerability(t.unif_10));
	EXPECT_PRED2(equal<eT>, 1,            MinEntropy<eT>().vulnerability(t.dirac_4));
	EXPECT_PRED2(equal<eT>, eT(8)/10, MinEntropy<eT>().vulnerability(t.pi1));
}

TYPED_TEST_P(MinEntropyTest, Cond_vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, 1,        MinEntropy<eT>(t.id_2).cond_vulnerability(t.unif_2));
	EXPECT_PRED2(equal<eT>, 1,        MinEntropy<eT>(t.id_2).cond_vulnerability(t.dirac_2));
	EXPECT_PRED2(equal<eT>, 1,        MinEntropy<eT>(t.id_2).cond_vulnerability(t.pi1));
	
	EXPECT_PRED2(equal<eT>, 1,        MinEntropy<eT>(t.id_10).cond_vulnerability(t.unif_10));
	EXPECT_PRED2(equal<eT>, 1,        MinEntropy<eT>(t.id_10).cond_vulnerability(t.dirac_10));
	EXPECT_PRED2(equal<eT>, 1,        MinEntropy<eT>(t.id_10).cond_vulnerability(t.pi2));
	
	EXPECT_PRED2(equal<eT>, eT(1)/10, MinEntropy<eT>(t.noint_10).cond_vulnerability(t.unif_10));
	EXPECT_PRED2(equal<eT>, 1,        MinEntropy<eT>(t.noint_10).cond_vulnerability(t.dirac_10));

	EXPECT_PRED2(equal<eT>, MinEntropy<eT>().vulnerability(t.pi2), MinEntropy<eT>(t.noint_10).cond_vulnerability(t.pi2));

	EXPECT_PRED2(equal<eT>, MinEntropy<eT>().vulnerability(t.pi3), MinEntropy<eT>(t.c1).cond_vulnerability(t.pi3));	// no change in vulnerability
	EXPECT_PRED2(equal<eT>, eT(31)/40, MinEntropy<eT>(t.c1).cond_vulnerability(t.pi4));

	ASSERT_ANY_THROW(MinEntropy<eT>(t.id_10).cond_vulnerability(t.unif_2));
}

TYPED_TEST_P(MinEntropyTest, Max_mult_leakage) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, 2,         MinEntropy<eT>(t.id_2).max_mult_leakage());
	EXPECT_PRED2(equal<eT>, 10,        MinEntropy<eT>(t.id_10).max_mult_leakage());
	EXPECT_PRED2(equal<eT>, 1,         MinEntropy<eT>(t.noint_10).max_mult_leakage());
	EXPECT_PRED2(equal<eT>, eT(15)/10, MinEntropy<eT>(t.c1).max_mult_leakage());

	EXPECT_PRED2(equal<eT>, MinEntropy<eT>(t.crand_10).mult_leakage(t.unif_10), MinEntropy<eT>(t.crand_10).max_mult_leakage()); // max_mult_leakage is given for uniform prior
}


TYPED_TEST_P(MinEntropyTestReals, Entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, -qif::log2(0.5), MinEntropy<eT>().entropy(t.unif_2));
	EXPECT_PRED2(equal<eT>, -qif::log2(0.1), MinEntropy<eT>().entropy(t.unif_10));
	EXPECT_PRED2(equal<eT>, 0,               MinEntropy<eT>().entropy(t.dirac_4));
	EXPECT_PRED2(equal<eT>, -qif::log2(0.8), MinEntropy<eT>().entropy(t.pi1));
}

TYPED_TEST_P(MinEntropyTestReals, Cond_entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, 0,   MinEntropy<eT>(t.id_2).cond_entropy(t.unif_2));
	EXPECT_PRED2(equal<eT>, 0,   MinEntropy<eT>(t.id_2).cond_entropy(t.dirac_2));
	EXPECT_PRED2(equal<eT>, 0,   MinEntropy<eT>(t.id_2).cond_entropy(t.pi1));
	
	EXPECT_PRED2(equal<eT>, 0,   MinEntropy<eT>(t.id_10).cond_entropy(t.unif_10));
	EXPECT_PRED2(equal<eT>, 0,   MinEntropy<eT>(t.id_10).cond_entropy(t.dirac_10));
	EXPECT_PRED2(equal<eT>, 0,   MinEntropy<eT>(t.id_10).cond_entropy(t.pi2));
	
	EXPECT_PRED2(equal<eT>, -qif::log2(0.1), MinEntropy<eT>(t.noint_10).cond_entropy(t.unif_10));
	EXPECT_PRED2(equal<eT>, 0,               MinEntropy<eT>(t.noint_10).cond_entropy(t.dirac_10));

	EXPECT_PRED2(equal<eT>, MinEntropy<eT>().entropy(t.pi2), MinEntropy<eT>(t.noint_10).cond_entropy(t.pi2));

	EXPECT_PRED2(equal<eT>, MinEntropy<eT>().entropy(t.pi3), MinEntropy<eT>(t.c1).cond_entropy(t.pi3)); // no change in entropy
	EXPECT_PRED2(equal<eT>, -qif::log2(0.775),               MinEntropy<eT>(t.c1).cond_entropy(t.pi4));

	ASSERT_ANY_THROW(MinEntropy<eT>(t.id_10).cond_entropy(t.unif_2));
}

TYPED_TEST_P(MinEntropyTestReals, Capacity) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal<eT>, 1,              MinEntropy<eT>(t.id_2).capacity());
	EXPECT_PRED2(equal<eT>, qif::log2(10),  MinEntropy<eT>(t.id_10).capacity());
	EXPECT_PRED2(equal<eT>, 0,              MinEntropy<eT>(t.noint_10).capacity());
	EXPECT_PRED2(equal<eT>, qif::log2(1.5), MinEntropy<eT>(t.c1).capacity());

	EXPECT_PRED2(equal<eT>, MinEntropy<eT>(t.crand_100).leakage(t.unif_100), MinEntropy<eT>(t.crand_100).capacity()); // capacity is given for uniform prior
}

// run the MinEntropyTest test-case for all types, and the MinEntropyTestReals only for double/float
//
REGISTER_TYPED_TEST_CASE_P(MinEntropyTest, Vulnerability, Cond_vulnerability, Max_mult_leakage);
REGISTER_TYPED_TEST_CASE_P(MinEntropyTestReals, Entropy, Cond_entropy, Capacity);

INSTANTIATE_TYPED_TEST_CASE_P(MinEntropy, MinEntropyTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(MinEntropyReals, MinEntropyTestReals, NativeTypes);

