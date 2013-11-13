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
#include "Channel.h"
#include "gtest/gtest.h"
#include <string>

using namespace std;

TEST(Channel, CorrectSizeAndElements) {
	string new_channel_elements = "1 0 0;0 1 0";
	Channel new_channel = Channel (new_channel_elements);
	EXPECT_EQ(2, new_channel.inputs_number());
	EXPECT_EQ(3, new_channel.outputs_number());
	for (int i = 0; i < new_channel.inputs_number(); ++i) {
		for (int j = 0; j < new_channel.outputs_number(); ++j) {
			if(i==j){
				EXPECT_EQ(1, new_channel.at(i,j));
			}else{
				EXPECT_EQ(0, new_channel.at(i,j));
			}
		}
	}
}

TEST(Channel, NoCorrectElements) {
	string new_channel_elements = "1 2;3 0.5";
	ASSERT_ANY_THROW(Channel new_channel = Channel(new_channel_elements););
}

TEST(new_id_channel, Negative) {
	ASSERT_ANY_THROW(Channel new_channel = new_channel.new_id_channel(-1););
}

TEST(new_id_channel, Zero) {
	Channel new_channel = new_channel.new_id_channel(0);
	EXPECT_EQ(0, new_channel.inputs_number());
	EXPECT_EQ(0, new_channel.outputs_number());
}

TEST(new_id_channel, Positive) {
	Channel new_channel = new_channel.new_id_channel(3);
	EXPECT_EQ(3, new_channel.inputs_number());
	EXPECT_EQ(3, new_channel.outputs_number());
	for (int i = 0; i < new_channel.inputs_number(); ++i) {
		for (int j = 0; j < new_channel.outputs_number(); ++j) {
			if(i==j){
				EXPECT_EQ(1, new_channel.at(i,j));
			}else{
				EXPECT_EQ(0, new_channel.at(i,j));
			}
		}
	}
}

/* Untested functions:
Channel (std::string& new_channel_elements);
~Channel();
Channel clone();
bool is_symmetric();
bool is_partial_symmetric();
bool is_equal_to ( const Channel& other );
vec get_row ( int index );
vec get_column ( int index );
void set_row ( int index,vec new_row_elements );
void set_row ( int index,char* new_row_elements );		
*/