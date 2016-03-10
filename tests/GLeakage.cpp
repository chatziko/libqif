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

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class GLeakageTest : public BaseTest<eT> {};
template <typename eT>
class GLeakageTestReals : public BaseTest<eT> {};

TYPED_TEST_CASE_P(GLeakageTest);
TYPED_TEST_CASE_P(GLeakageTestReals);		// tests that run only on double/float


// TODO: test more gain functions. Currently these tests are copies from MinEntropy and test only the id gain function


TYPED_TEST_P(GLeakageTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, eT(5)/10, GLeakage<eT>(t.id_2,  t.id_2).vulnerability(t.unif_2));
	EXPECT_PRED2(equal2<eT>, eT(1)/10, GLeakage<eT>(t.id_10, t.id_10).vulnerability(t.unif_10));
	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.id_4,  t.id_4).vulnerability(t.dirac_4));
	EXPECT_PRED2(equal2<eT>, eT(8)/10, GLeakage<eT>(t.id_2,  t.id_2).vulnerability(t.pi1));
}

TYPED_TEST_P(GLeakageTest, Cond_vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.id_2,     t.id_2).cond_vulnerability(t.unif_2));
	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.id_2,     t.id_2).cond_vulnerability(t.dirac_2));
	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.id_2,     t.id_2).cond_vulnerability(t.pi1));
	
	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.id_10,    t.id_10).cond_vulnerability(t.unif_10));
	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.id_10,    t.id_10).cond_vulnerability(t.dirac_10));
	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.id_10,    t.id_10).cond_vulnerability(t.pi2));
	
	EXPECT_PRED2(equal2<eT>, eT(1)/10, GLeakage<eT>(t.noint_10, t.id_10).cond_vulnerability(t.unif_10));
	EXPECT_PRED2(equal2<eT>, 1,        GLeakage<eT>(t.noint_10, t.id_10).cond_vulnerability(t.dirac_10));

	EXPECT_PRED2(equal2<eT>, GLeakage<eT>(t.noint_10, t.id_10).vulnerability(t.pi2), GLeakage<eT>(t.noint_10, t.id_10).cond_vulnerability(t.pi2));

	EXPECT_PRED2(equal2<eT>, GLeakage<eT>(t.c1, t.id_2).vulnerability(t.pi3), GLeakage<eT>(t.c1, t.id_2).cond_vulnerability(t.pi3));	// no change in vulnerability
	EXPECT_PRED2(equal2<eT>, eT(31)/40, GLeakage<eT>(t.c1, t.id_2).cond_vulnerability(t.pi4));

	ASSERT_ANY_THROW(GLeakage<eT>(t.id_10, t.id_10).cond_vulnerability(t.unif_2));
}


TYPED_TEST_P(GLeakageTestReals, Entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, -internal::log2(0.5), GLeakage<eT>(t.id_2,  t.id_2).entropy(t.unif_2));
	EXPECT_PRED2(equal2<eT>, -internal::log2(0.1), GLeakage<eT>(t.id_10, t.id_10).entropy(t.unif_10));
	EXPECT_PRED2(equal2<eT>, 0,               GLeakage<eT>(t.id_4,  t.id_4).entropy(t.dirac_4));
	EXPECT_PRED2(equal2<eT>, -internal::log2(0.8), GLeakage<eT>(t.id_2,  t.id_2).entropy(t.pi1));
}

TYPED_TEST_P(GLeakageTestReals, Cond_entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, 0,   GLeakage<eT>(t.id_2, t.id_2).cond_entropy(t.unif_2));
	EXPECT_PRED2(equal2<eT>, 0,   GLeakage<eT>(t.id_2, t.id_2).cond_entropy(t.dirac_2));
	EXPECT_PRED2(equal2<eT>, 0,   GLeakage<eT>(t.id_2, t.id_2).cond_entropy(t.pi1));
	
	EXPECT_PRED2(equal2<eT>, 0,   GLeakage<eT>(t.id_10, t.id_10).cond_entropy(t.unif_10));
	EXPECT_PRED2(equal2<eT>, 0,   GLeakage<eT>(t.id_10, t.id_10).cond_entropy(t.dirac_10));
	EXPECT_PRED2(equal2<eT>, 0,   GLeakage<eT>(t.id_10, t.id_10).cond_entropy(t.pi2));
	
	EXPECT_PRED2(equal2<eT>, -internal::log2(0.1), GLeakage<eT>(t.noint_10, t.id_10).cond_entropy(t.unif_10));
	EXPECT_PRED2(equal2<eT>, 0,               GLeakage<eT>(t.noint_10, t.id_10).cond_entropy(t.dirac_10));

	EXPECT_PRED2(equal2<eT>, GLeakage<eT>(t.noint_10, t.id_10).entropy(t.pi2), GLeakage<eT>(t.noint_10, t.id_10).cond_entropy(t.pi2));

	EXPECT_PRED2(equal2<eT>, GLeakage<eT>(t.c1, t.id_2).entropy(t.pi3), GLeakage<eT>(t.c1, t.id_2).cond_entropy(t.pi3)); // no change in entropy
	EXPECT_PRED2(equal2<eT>, -internal::log2(0.775), GLeakage<eT>(t.c1, t.id_2).cond_entropy(t.pi4));

	ASSERT_ANY_THROW(GLeakage<eT>(t.id_10, t.id_2).cond_entropy(t.unif_2));
}

// run the GLeakageTest test-case for all types, and the GLeakageTestReals only for double/float
//
REGISTER_TYPED_TEST_CASE_P(GLeakageTest, Vulnerability, Cond_vulnerability);
REGISTER_TYPED_TEST_CASE_P(GLeakageTestReals, Entropy, Cond_entropy);

INSTANTIATE_TYPED_TEST_CASE_P(GLeakage, GLeakageTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(GLeakageReals, GLeakageTestReals, NativeTypes);

