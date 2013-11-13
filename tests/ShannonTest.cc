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
#include <string>

using namespace std;

TEST(S_vulnerability, not_supported) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	Shannon shannon = Shannon(new_channel);
	string new_prior_elements ="0.5 0.5";
	Prob prior = Prob(new_prior_elements);
	ASSERT_ANY_THROW(double r= shannon.vulnerability(prior));
}

TEST(S_cond_vulnerability, not_supported) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	Shannon shannon = Shannon(new_channel);
	string new_prior_elements ="0.5 0.5";
	Prob prior = Prob(new_prior_elements);
	ASSERT_ANY_THROW(double r= shannon.cond_vulnerability(prior););
}

TEST(Shannon_QIF_functions, incorrect_X_size) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	Shannon shannon = Shannon(new_channel);
	string new_prior_elements ="0.33 0.33 0.34";
	Prob prior = Prob(new_prior_elements);
	ASSERT_ANY_THROW(double r= shannon.entropy(prior););
	ASSERT_ANY_THROW(double r= shannon.cond_entropy(prior););
	ASSERT_ANY_THROW(double r= shannon.leakage(prior););
}

/* Untested functions:
Shannon(Channel c);
~Shannon();
double leakage(Prob pi);
double entropy(Prob pi);
double cond_entropy(Prob pi);
double capacity();
*/