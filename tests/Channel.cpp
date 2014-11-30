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
#include <type_traits>


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ChannelTest : public ::testing::Test {};

TYPED_TEST_CASE_P(ChannelTest);


TYPED_TEST_P(ChannelTest, constructors) {
	typedef TypeParam eT;

	const char* s = "1 0 0; 0 1 0";
	Mat<eT> m(s);
	Channel<eT> c(s);

	expect_channel( "", Channel<eT>()               ); // empty
	expect_channel( s,  Channel<eT>(s)              ); // char*
	expect_channel( s,  Channel<eT>(std::string(s)) ); // std::string
	expect_channel( s,  Channel<eT>(m)              ); // Mat
	expect_channel( s,  Channel<eT>(c)              ); // copy
	expect_channel( s,  Channel<eT>(Channel<eT>(s)) ); // move (Note: most likely the compiler is removing the move call completely, see http://en.cppreference.com/w/cpp/language/copy_elision)
	expect_channel( s,  Channel<eT>(Mat<eT>(s))     ); // Mat, move semantics

	// malformed channel
	// Note: "cout <<" is to avoid the compiler removing the code as unused!
	//
	const char* s2 = "1 2; 3 0.5";
	Mat<eT> m2(s2);

	EXPECT_ANY_THROW( std::cout << Channel<eT>(s2);              ); // char*
	EXPECT_ANY_THROW( std::cout << Channel<eT>(std::string(s2)); ); // std::string
	EXPECT_ANY_THROW( std::cout << Channel<eT>(m2);              ); // Mat
	EXPECT_ANY_THROW( std::cout << Channel<eT>(Mat<eT>(s2));     ); // Mat, move semantics
}

TYPED_TEST_P(ChannelTest, identity) {
	typedef TypeParam eT;

	Channel<eT> c;
	c.identity(0);
	expect_channel(0, 0, c);

	c.identity(3);
	expect_channel("1 0 0; 0 1 0; 0 0 1", c);
}

TYPED_TEST_P(ChannelTest, randu) {
	typedef TypeParam eT;

	Channel<eT> c(200, 200);
	c.randu();
	expect_channel(200, 200, c);

	c.randu(5);
	expect_channel(5, 5, c);

	c.randu(4, 6);
	expect_channel(4, 6, c);
}



// run the ChannelTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_CASE_P(ChannelTest, constructors, identity, randu);

typedef ::testing::Types<double, float, urat> ChannelTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Channel, ChannelTest, ChannelTypes);

