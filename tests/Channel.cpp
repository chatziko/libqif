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
#include "tests_aux.h"
#include "gtest/gtest.h"
#include <string>

using namespace std;

TEST(chan, constructors) {
	const char* s = "1 0 0; 0 1 0";
	mat  m(s);
	chan c(s);

	expect_channel( "", chan()               ); // empty
	expect_channel( s,  chan(s)              ); // char*
	expect_channel( s,  chan(std::string(s)) ); // std::string
	expect_channel( s,  chan(m)              ); // mat
	expect_channel( s,  chan(c)              ); // copy
	expect_channel( s,  chan(chan(s))        ); // move (Note: most likely the compiler is removing the move call completely, see http://en.cppreference.com/w/cpp/language/copy_elision)
	expect_channel( s,  chan(mat(s))         ); // mat, move semantics

	// malformed channel
	// Note: "cout <<" is to avoid the compiler removing the code as unused!
	//
	const char* s2 = "1 2;3 0.5";
	mat m2(s2);

	EXPECT_ANY_THROW( cout << chan(s2);              ); // char*
	EXPECT_ANY_THROW( cout << chan(std::string(s2)); ); // std::string
	EXPECT_ANY_THROW( cout << chan(m2);              ); // mat
	EXPECT_ANY_THROW( cout << chan(mat(s2));         ); // mat, move semantics
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

	for(uint i = 0; i < new_channel.n_rows; ++i)
		for(uint j = 0; j < new_channel.n_cols; ++j)
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
