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
	Channel new_channel = Channel (new_channel_elements);
	MinEntropy minEntropy = MinEntropy(new_channel);
	string new_prior_elements ="0.33 0.33 0.34";
	Prob prior = Prob(new_prior_elements);
	ASSERT_ANY_THROW(double r= minEntropy.vulnerability(prior););
	ASSERT_ANY_THROW(double r= minEntropy.cond_vulnerability(prior););
	ASSERT_ANY_THROW(double r= minEntropy.entropy(prior););
	ASSERT_ANY_THROW(double r= minEntropy.cond_entropy(prior););
	ASSERT_ANY_THROW(double r= minEntropy.leakage(prior););
}

/* Untested functions:
MinEntropy(Channel c);
~MinEntropy();
double vulnerability(Prob pi);
double cond_vulnerability(Prob pi);
double leakage(Prob pi);
double entropy(Prob pi);
double cond_entropy(Prob pi);
double capacity();
*/