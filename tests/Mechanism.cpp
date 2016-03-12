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
class MechanismTest : public BaseTest<eT> {};

TYPED_TEST_CASE_P(MechanismTest);


TYPED_TEST_P(MechanismTest, Is_private) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Mechanism<eT> mech;

	mech.C = t.noint_10;
	EXPECT_TRUE(mech.is_private(0));

	mech.C = t.id_10;
	EXPECT_FALSE(mech.is_private(10000));

	mech.C = t.c1;
	EXPECT_TRUE (mech.is_private(std::log(7.0/2)));
	EXPECT_FALSE(mech.is_private(std::log(6.9/2)));
}

TYPED_TEST_P(MechanismTest, Smallest_epsilon) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Mechanism<eT> mech;

	mech.C = t.noint_10;
	EXPECT_PRED2(equal2<eT>, eT(0), mech.smallest_epsilon());

	mech.C = t.id_10;
	EXPECT_PRED2(equal2<eT>, infinity<eT>(), mech.smallest_epsilon());

	mech.C = t.c1;
	EXPECT_PRED2(equal2<eT>, std::log(7.0/2), mech.smallest_epsilon());
}

TYPED_TEST_P(MechanismTest, Reals) {
	typedef TypeParam eT;

	uint size = 10;
	eT step = 1.1;
	eT epsilon = 0.9;

	auto d = metric::scale(metric::euclidean<eT, uint>(), step);

	Mechanism<eT> geom = mechanisms::geometric<eT>(size, step, epsilon);
	Mechanism<eT> expon = mechanisms::exponential<eT>(size, d, epsilon);
	Mechanism<eT> tc = mechanisms::tight_constraints<eT>(size, d, epsilon);

	expect_channel(size, size, geom.C);
	EXPECT_TRUE(geom.is_private(epsilon));
	EXPECT_FALSE(geom.is_private(epsilon - 0.01));
	EXPECT_PRED2(equal2<eT>, epsilon, geom.smallest_epsilon());

	expect_channel(size, size, expon.C);
	EXPECT_TRUE(expon.is_private(epsilon));
	EXPECT_PRED2(equal2<eT>, 0.64197307180467134, expon.smallest_epsilon());

	expect_channel(geom.C, tc.C);
}

TYPED_TEST_P(MechanismTest, Discrete) {
	typedef TypeParam eT;

	uint size = 5;
	eT epsilon = 0.9;

	auto d = metric::discrete<eT, uint>();

	Mechanism<eT> tc = mechanisms::tight_constraints<eT>(size, d, epsilon);
	Mechanism<eT> expon = mechanisms::exponential<eT>(size, d, 2*epsilon);

	expect_channel(size, size, tc.C);
	EXPECT_TRUE(tc.is_private(epsilon));
	EXPECT_FALSE(tc.is_private(epsilon - 0.01));
	EXPECT_PRED2(equal2<eT>, epsilon, tc.smallest_epsilon());

	// the exponential mechanism with 2*epsilon should be the same as the tight constraints
	EXPECT_PRED2(chan_equal2<eT>, expon.C, tc.C);
}

TYPED_TEST_P(MechanismTest, Grid) {
	typedef TypeParam eT;

	uint width = 3,
		 height = 3,
		 size = width * height;
	eT step = 1.1,
	   epsilon = 0.9;

	auto d = metric::scale(metric::grid<eT, Point<eT>>(width), step);

	Mechanism<eT> laplace = mechanisms::planar_laplace_grid<eT>(width, height, step, epsilon);
	Mechanism<eT> tc = mechanisms::tight_constraints<eT>(size, d, epsilon);
	Mechanism<eT> expon = mechanisms::exponential<eT>(size, d, epsilon);

	expect_channel(size, size, tc.C);
	EXPECT_TRUE(tc.is_private(epsilon));
	EXPECT_FALSE(tc.is_private(epsilon - 0.01));
	EXPECT_PRED2(equal2<eT>, epsilon, tc.smallest_epsilon());

	// for the planar laplace, the accuracy that we get through numeric integration is not that great
	EXPECT_PRED2(channel::is_proper<eT>, laplace.C, eT(1e-3));
	EXPECT_TRUE(laplace.is_private(epsilon));
	EXPECT_PRED4(equal<eT>, 0.8844, laplace.smallest_epsilon(), 0, eT(1e-3));

	expect_channel(size, size, expon.C);
	EXPECT_TRUE(expon.is_private(epsilon));
	EXPECT_PRED2(equal2<eT>, 0.58945591528726249, expon.smallest_epsilon());
}



REGISTER_TYPED_TEST_CASE_P(MechanismTest, Is_private, Smallest_epsilon, Reals, Discrete, Grid);

INSTANTIATE_TYPED_TEST_CASE_P(Mechanism, MechanismTest, NativeTypes);

