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
#include "GLeakage.h"
#include <string>

using namespace std;

TEST(GLeakage, incorrect_X_size) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	string new_gain_elements = "1 0 0;0 1 0";
	Gain new_gain = Gain (new_gain_elements);
	ASSERT_ANY_THROW(GLeakage gleakage = GLeakage(new_channel,new_gain););
}

TEST(GLeakage_QIF_functions, incorrect_X_size) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	string new_gain_elements = "1 0;0 1";
	Gain new_gain = Gain (new_gain_elements);
	GLeakage gleakage = GLeakage(new_channel,new_gain);
	string new_prior_elements ="0.33 0.33 0.34";
	Prob prior = Prob(new_prior_elements);
	ASSERT_ANY_THROW(double r= gleakage.vulnerability(prior););
	ASSERT_ANY_THROW(double r= gleakage.cond_vulnerability(prior););
	ASSERT_ANY_THROW(double r= gleakage.entropy(prior););
	ASSERT_ANY_THROW(double r= gleakage.cond_entropy(prior););
	ASSERT_ANY_THROW(double r= gleakage.leakage(prior););
}

TEST(GLeakage_plotter_functions, plot_without_selected_engine) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	string new_gain_elements = "1 0;0 1";
	Gain new_gain = Gain (new_gain_elements);
	GLeakage gleakage = GLeakage(new_channel,new_gain);
	ASSERT_ANY_THROW(gleakage.plot2d_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot2d_cond_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot2d_leakage(););
	ASSERT_ANY_THROW(gleakage.plot2d_entropy(););
	ASSERT_ANY_THROW(gleakage.plot2d_cond_entropy(););

	//------------------------------------------------------
	new_channel_elements = "1 0 0;0 1 0;0 0 1";
	new_channel = Channel (new_channel_elements);
	new_gain_elements = "1 0 0;0 1 0;0 0 1";
	new_gain = Gain (new_gain_elements);
	gleakage = GLeakage(new_channel,new_gain);
	ASSERT_ANY_THROW(gleakage.plot3d_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot3d_cond_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot3d_leakage(););
	ASSERT_ANY_THROW(gleakage.plot3d_entropy(););
	ASSERT_ANY_THROW(gleakage.plot3d_cond_entropy(););
}

TEST(GLeakage_plotter_functions, incorrect_X_size) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	string new_gain_elements = "1 0;0 1";
	Gain new_gain = Gain (new_gain_elements);
	GLeakage gleakage = GLeakage(new_channel,new_gain);
	gleakage.change_to_scilab();
	ASSERT_ANY_THROW(gleakage.plot3d_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot3d_cond_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot3d_leakage(););
	ASSERT_ANY_THROW(gleakage.plot3d_entropy(););
	ASSERT_ANY_THROW(gleakage.plot3d_cond_entropy(););
	//--------------------------------------------------
	new_channel_elements = "1 0 0;0 1 0;0 0 1";
	new_channel = Channel (new_channel_elements);
	new_gain_elements = "1 0 0;0 1 0;0 0 1";
	new_gain = Gain (new_gain_elements);
	gleakage = GLeakage(new_channel,new_gain);
	gleakage.change_to_scilab();
	ASSERT_ANY_THROW(gleakage.plot2d_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot2d_cond_vulnerability(););
	ASSERT_ANY_THROW(gleakage.plot2d_leakage(););
	ASSERT_ANY_THROW(gleakage.plot2d_entropy(););
	ASSERT_ANY_THROW(gleakage.plot2d_cond_entropy(););
}

TEST(G_capacity, not_supported) {
	string new_channel_elements = "1 0;0 1";
	Channel new_channel = Channel (new_channel_elements);
	string new_gain_elements = "1 0;0 1";
	Gain new_gain = Gain (new_gain_elements);
	GLeakage gleakage = GLeakage(new_channel,new_gain);
	ASSERT_ANY_THROW(gleakage.capacity(););
}

/* Untested functions:
double vulnerability(Prob pi);
double cond_vulnerability(Prob pi);
double leakage(Prob pi);
double entropy(Prob pi);
double cond_entropy(Prob pi);
double capacity();
void * compare_over_prior(Channel other_channel);
void * compare_over_gain(Channel other_channel,Prob prior);
*/