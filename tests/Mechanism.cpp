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
#include "Mechanism.h"
#include "aux.h"
#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MechanismTest : public BaseTest<eT> {};

TYPED_TEST_CASE_P(MechanismTest);


TYPED_TEST_P(MechanismTest, Is_private) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_TRUE(Mechanism<eT>(t.noint_10).is_private(0));

//	EXPECT_PRED2(equal2<eT>, eT(1)/10, Mechanism<eT>().vulnerability(t.unif_10));
//	EXPECT_PRED2(equal2<eT>, 1,            Mechanism<eT>().vulnerability(t.dirac_4));
//	EXPECT_PRED2(equal2<eT>, eT(8)/10, Mechanism<eT>().vulnerability(t.pi1));
}

TYPED_TEST_P(MechanismTest, Smallest_epsilon) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, eT(0), Mechanism<eT>(t.noint_10).smallest_epsilon());

//	EXPECT_PRED2(equal2<eT>, eT(1)/10, Mechanism<eT>().vulnerability(t.unif_10));
//	EXPECT_PRED2(equal2<eT>, 1,            Mechanism<eT>().vulnerability(t.dirac_4));
//	EXPECT_PRED2(equal2<eT>, eT(8)/10, Mechanism<eT>().vulnerability(t.pi1));
}


REGISTER_TYPED_TEST_CASE_P(MechanismTest, Is_private, Smallest_epsilon);

INSTANTIATE_TYPED_TEST_CASE_P(Mechanism, MechanismTest, NativeTypes);

