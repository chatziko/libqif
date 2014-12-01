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

using namespace std;

TEST(Shannon, vulnerability) {
	prob pi("0.5 0.5");
	Shannon shan(chan("1 0; 0 1"));

	ASSERT_ANY_THROW( shan.vulnerability(pi); );
	ASSERT_ANY_THROW( shan.cond_vulnerability(pi); );
}

TEST(Shannon, entropy) {
	prob pi;
	Shannon shan;

	pi = uniform<prob>(2);
	EXPECT_EQ(1, shan.entropy(pi));

	pi = uniform<prob>(10);
	EXPECT_FLOAT_EQ(log2(10), shan.entropy(pi));

	pi = prob("1 0 0 0");
	EXPECT_EQ(0, shan.entropy(pi));

	pi = prob("0.2 0.8");
	EXPECT_FLOAT_EQ(0.721928094887362, shan.entropy(pi));
}

TEST(Shannon, cond_entropy) {
	Shannon shan;

	shan.C = identity<chan>(2);
	EXPECT_EQ(0, shan.cond_entropy(uniform<prob>(2)));
	EXPECT_EQ(0, shan.cond_entropy(dirac<prob>(2)));
	EXPECT_EQ(0, shan.cond_entropy(prob("0.2 0.8")));

	shan.C = identity<chan>(10);
	EXPECT_EQ(0, shan.cond_entropy(uniform<prob>(10)));
	EXPECT_EQ(0, shan.cond_entropy(dirac<prob>(10)));
	EXPECT_EQ(0, shan.cond_entropy(prob("0.2 0.8 0 0 0 0 0 0 0 0")));

	no_interference(shan.C);
	EXPECT_FLOAT_EQ(log2(10), shan.cond_entropy(uniform<prob>(10)));
	EXPECT_EQ(0, shan.cond_entropy(dirac<prob>(10)));

	prob pi = randu<prob>(10);
	EXPECT_FLOAT_EQ(shan.entropy(pi), shan.cond_entropy(pi));

	shan.C = chan("0.8 0.2; 0.3 0.7");
	pi = "0.25 0.75";
	EXPECT_FLOAT_EQ(0.669020059980807, shan.cond_entropy(pi));

	shan.C = identity<chan>(10);
	ASSERT_ANY_THROW(shan.cond_entropy(uniform<prob>(2)););
}

TEST(Shannon, capacity) {
	Shannon shan;

	shan.C = identity<chan>(2);
	EXPECT_EQ(1, shan.capacity());

	shan.C = identity<chan>(10);
	EXPECT_FLOAT_EQ(log2(10), shan.capacity());

	shan.C = no_interference<chan>(10);
	EXPECT_EQ(0, shan.capacity());

	shan.C = chan("0.8 0.2; 0.3 0.7");
	EXPECT_FLOAT_EQ(0.19123721482206, shan.capacity());

	// symmetric
	shan.C = chan(
		".3 .2 .5;"
		".5 .3 .2;"
		".2 .5 .3;"
	);
	double cap = log2(shan.C.n_cols) - shan.entropy(shan.C.row(0));
	EXPECT_FLOAT_EQ(cap, shan.capacity());

	// weakly symmetric
	shan.C = chan(
		"0.333333333 0.166666667 0.5;"
		"0.333333333 0.5         0.166666667;"
	);
	cap = log2(shan.C.n_cols) - shan.entropy(shan.C.row(0));
	EXPECT_FLOAT_EQ(cap, shan.capacity());
}
