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
#include "Gain.h"
#include "gtest/gtest.h"
#include <string>

using namespace std;

TEST(Gain, CorrectSizeAndElements) {
	string new_gain_elements = "1 0 0; 0 1 0";
	Gain new_gain = Gain (new_gain_elements);
	EXPECT_EQ(3, new_gain.inputs_number());
	EXPECT_EQ(2, new_gain.guesses_number());
	for (int i = 0; i < new_gain.inputs_number(); ++i) {
		for (int j = 0; j < new_gain.guesses_number(); ++j) {
			if(i==j){
				EXPECT_EQ(1, new_gain.at(i,j));
			}else{
				EXPECT_EQ(0, new_gain.at(i,j));
			}
		}
	}
}

TEST(Gain, NoCorrectElements) {
	string new_gain_elements = "1 2;3 0.5 2";
	ASSERT_ANY_THROW(Gain new_gain = Gain (new_gain_elements););
}

TEST(G_new_id_function, Negative) {
	ASSERT_ANY_THROW(Gain new_gain = new_gain.new_id_function(-1););
}

TEST(G_new_id_function, Zero) {
	Gain new_gain = new_gain.new_id_function(0);
	EXPECT_EQ(0, new_gain.inputs_number());
	EXPECT_EQ(0, new_gain.guesses_number());
}

TEST(G_new_id_function, Positive) {
	Gain new_gain = new_gain.new_id_function(3);
	EXPECT_EQ(3, new_gain.inputs_number());
	EXPECT_EQ(3, new_gain.guesses_number());
	for (int i = 0; i < new_gain.inputs_number(); ++i) {
		for (int j = 0; j < new_gain.guesses_number(); ++j) {
			if(i==j){
				EXPECT_EQ(1, new_gain.at(i,j));
			}else{
				EXPECT_EQ(0, new_gain.at(i,j));
			}
		}
	}
}

/* Untested functions:
~Gain();
Gain zeros ( int inputs, int guesses );
Gain ones ( int inputs, int guesses );
Gain clone();
bool is_symmetric();
bool is_partial_symmetric();
bool is_equal_to ( const Gain& other );
vec get_row ( int index );
vec get_column ( int index );
void set_row ( int index,vec new_row_elements );
void set_row ( int index,char* new_row_elements );
*/