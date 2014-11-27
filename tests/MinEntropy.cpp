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
#include <string>

using namespace std;

TEST(MinEntropy_QIF_functions, incorrect_X_size) {
	string new_channel_elements = "1 0;0 1";
	chan new_channel = chan(new_channel_elements);
	MinEntropy minEntropy = MinEntropy(new_channel);
	string new_prior_elements = "0.33 0.33 0.34";
	prob prior = prob(new_prior_elements);
	ASSERT_ANY_THROW(minEntropy.vulnerability(prior););
	ASSERT_ANY_THROW(minEntropy.cond_vulnerability(prior););
	ASSERT_ANY_THROW(minEntropy.entropy(prior););
	ASSERT_ANY_THROW(minEntropy.cond_entropy(prior););
	ASSERT_ANY_THROW(minEntropy.leakage(prior););
}

/* Untested functions:
MinEntropy(chan c);
~MinEntropy();
double vulnerability(prob pi);
double cond_vulnerability(prob pi);
double leakage(prob pi);
double entropy(prob pi);
double cond_entropy(prob pi);
double capacity();
*/
