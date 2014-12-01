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
template <typename T>
class MinEntropyTest : public ::testing::Test {};

TYPED_TEST_CASE_P(MinEntropyTest);



TYPED_TEST_P(MinEntropyTest, Vulnerability) {
	typedef TypeParam eT;

	Prob<eT> pi;
	MinEntropy<eT> minent;

	pi = uniform<Prob<eT>>(2);
	EXPECT_PRED2(equal<eT>, 0.5, minent.vulnerability(pi));

	pi = uniform<Prob<eT>>(10);
	EXPECT_PRED2(equal<eT>, 0.1, minent.vulnerability(pi));

	pi = Prob<eT>("1 0 0 0");
	EXPECT_PRED2(equal<eT>, 1, minent.vulnerability(pi));

	pi = Prob<eT>("0.2 0.8");
	EXPECT_PRED2(equal<eT>, 0.8, minent.vulnerability(pi));
}

TYPED_TEST_P(MinEntropyTest, Entropy) {
	typedef TypeParam eT;

	Prob<eT> pi;
	MinEntropy<eT> minent;

	pi = uniform<Prob<eT>>(2);
	EXPECT_PRED2(equal<eT>, -qif::log2(0.5), minent.entropy(pi));

	pi = uniform<Prob<eT>>(10);
	EXPECT_PRED2(equal<eT>, -qif::log2(0.1), minent.entropy(pi));

	pi = Prob<eT>("1 0 0 0");
	EXPECT_PRED2(equal<eT>, 0, minent.entropy(pi));

	pi = Prob<eT>("0.2 0.8");
	EXPECT_PRED2(equal<eT>, -qif::log2(0.8), minent.entropy(pi));
}

TYPED_TEST_P(MinEntropyTest, Cond_vulnerability) {
	typedef TypeParam eT;

	MinEntropy<eT> minent;

	minent.C = identity<Chan<eT>>(2);
	EXPECT_PRED2(equal<eT>, 1, minent.cond_vulnerability(uniform<Prob<eT>>(2)));
	EXPECT_PRED2(equal<eT>, 1, minent.cond_vulnerability(dirac<Prob<eT>>(2)));
	EXPECT_PRED2(equal<eT>, 1, minent.cond_vulnerability(Prob<eT>("0.2 0.8")));

	minent.C = identity<Chan<eT>>(10);
	EXPECT_PRED2(equal<eT>, 1, minent.cond_vulnerability(uniform<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 1, minent.cond_vulnerability(dirac<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 1, minent.cond_vulnerability(Prob<eT>("0.2 0.8 0 0 0 0 0 0 0 0")));

	no_interference(minent.C);
	EXPECT_PRED2(equal<eT>, 0.1, minent.cond_vulnerability(uniform<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 1, minent.cond_vulnerability(dirac<Prob<eT>>(10)));

	Prob<eT> pi = randu<Prob<eT>>(10);
	EXPECT_PRED2(equal<eT>, minent.vulnerability(pi), minent.cond_vulnerability(pi));

	minent.C = Chan<eT>("0.8 0.2; 0.3 0.7");
	pi = "0.25 0.75";
	EXPECT_PRED2(equal<eT>, minent.vulnerability(pi), minent.cond_vulnerability(pi));	// no change in vulnerability
	pi = "0.75 0.25";
	EXPECT_PRED2(equal<eT>, 0.775, minent.cond_vulnerability(pi));

	minent.C = identity<Chan<eT>>(10);
	ASSERT_ANY_THROW(minent.cond_vulnerability(uniform<Prob<eT>>(2)););
}

TYPED_TEST_P(MinEntropyTest, Cond_entropy) {
	typedef TypeParam eT;

	MinEntropy<eT> minent;

	minent.C = identity<Chan<eT>>(2);
	EXPECT_PRED2(equal<eT>, 0, minent.cond_entropy(uniform<Prob<eT>>(2)));
	EXPECT_PRED2(equal<eT>, 0, minent.cond_entropy(dirac<Prob<eT>>(2)));
	EXPECT_PRED2(equal<eT>, 0, minent.cond_entropy(Prob<eT>("0.2 0.8")));

	minent.C = identity<Chan<eT>>(10);
	EXPECT_PRED2(equal<eT>, 0, minent.cond_entropy(uniform<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 0, minent.cond_entropy(dirac<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 0, minent.cond_entropy(Prob<eT>("0.2 0.8 0 0 0 0 0 0 0 0")));

	no_interference(minent.C);
	EXPECT_PRED2(equal<eT>, qif::log2(10), minent.cond_entropy(uniform<Prob<eT>>(10)));
	EXPECT_PRED2(equal<eT>, 0, minent.cond_entropy(dirac<Prob<eT>>(10)));

	Prob<eT> pi = randu<Prob<eT>>(10);
	EXPECT_PRED2(equal<eT>, minent.entropy(pi), minent.cond_entropy(pi));

	minent.C = Chan<eT>("0.8 0.2; 0.3 0.7");
	pi = "0.25 0.75";
	EXPECT_PRED2(equal<eT>, minent.entropy(pi), minent.cond_entropy(pi));	// no change in entropy
	pi = "0.75 0.25";
	EXPECT_PRED2(equal<eT>, -qif::log2(0.775), minent.cond_entropy(pi));

	minent.C = identity<Chan<eT>>(10);
	ASSERT_ANY_THROW(minent.cond_entropy(uniform<Prob<eT>>(2)););
}

TYPED_TEST_P(MinEntropyTest, Capacity) {
	typedef TypeParam eT;

	MinEntropy<eT> minent;

	minent.C = identity<Chan<eT>>(2);
	EXPECT_PRED2(equal<eT>, 1, minent.capacity());

	minent.C = identity<Chan<eT>>(10);
	EXPECT_PRED2(equal<eT>, qif::log2(10), minent.capacity());

	minent.C = no_interference<Chan<eT>>(10);
	EXPECT_PRED2(equal<eT>, 0, minent.capacity());

	minent.C = Chan<eT>("0.8 0.2; 0.3 0.7");
	EXPECT_PRED2(equal<eT>, qif::log2(1.5), minent.capacity());

	minent.C = randu<Chan<eT>>(100);
	EXPECT_PRED2(equal<eT>, minent.leakage(uniform<Prob<eT>>(100)), minent.capacity());		// capacity is given for uniform prior
}


// run the ChanTest test-case for double, float
//
REGISTER_TYPED_TEST_CASE_P(MinEntropyTest, Vulnerability, Entropy, Cond_vulnerability, Cond_entropy, Capacity);

INSTANTIATE_TYPED_TEST_CASE_P(MinEntropy, MinEntropyTest, NativeTypes);

