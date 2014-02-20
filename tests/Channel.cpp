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

TEST(chan, CorrectSizeAndElements) {
	string new_channel_elements = "1 0 0;0 1 0";
	chan new_channel = chan(new_channel_elements);
	EXPECT_EQ(2, new_channel.n_rows);
	EXPECT_EQ(3, new_channel.n_cols);
	for(int i = 0; i < new_channel.n_rows; ++i) {
		for(int j = 0; j < new_channel.n_cols; ++j) {
			if(i == j) {
				EXPECT_EQ(1, new_channel.at(i, j));
			} else {
				EXPECT_EQ(0, new_channel.at(i, j));
			}
		}
	}
}

TEST(chan, NoCorrectElements) {
	string new_channel_elements = "1 2;3 0.5";
	ASSERT_ANY_THROW(chan new_channel = chan(new_channel_elements););
}

TEST(identity, Zero) {
	chan new_channel;
	new_channel.identity(0);
	EXPECT_EQ(0, new_channel.n_rows);
	EXPECT_EQ(0, new_channel.n_cols);
}

TEST(identity, Positive) {
	chan new_channel;
	new_channel.identity(3);
	EXPECT_EQ(3, new_channel.n_rows);
	EXPECT_EQ(3, new_channel.n_cols);

	for(int i = 0; i < new_channel.n_rows; ++i)
		for(int j = 0; j < new_channel.n_cols; ++j)
			EXPECT_EQ(i == j ? 1 : 0, new_channel.at(i, j));
}

/* Untested functions:
chan (std::string& new_channel_elements);
~chan();
chan clone();
bool is_symmetric();
bool is_partial_symmetric();
bool is_equal_to ( const chan& other );
*/
