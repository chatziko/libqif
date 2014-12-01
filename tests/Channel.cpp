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
#include <string>
#include "gtest/gtest.h"

#include "Channel.h"
#include "tests_aux.h"


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ChannelTest : public ::testing::Test {};

TYPED_TEST_CASE_P(ChannelTest);


TYPED_TEST_P(ChannelTest, Construct) {
	typedef TypeParam eT;

	const char* s = "1 0 0; 0 1 0";
	Channel<eT> C(s);

	expect_channel( s,  Channel<eT>(s)              ); // char*
	expect_channel( s,  Channel<eT>(std::string(s)) ); // std::string
	expect_channel( s,  Channel<eT>(C)              ); // copy

	// malformed channel
	//
	const char* s2 = "1 2; 3 0.5";
	Channel<eT> C2(s2);

	EXPECT_ANY_THROW( check_proper(Channel<eT>(s2));              ); // char*
	EXPECT_ANY_THROW( check_proper(Channel<eT>(std::string(s2))); ); // std::string
	EXPECT_ANY_THROW( check_proper(Channel<eT>(C2));              ); // Mat
}

TYPED_TEST_P(ChannelTest, Identity) {
	typedef TypeParam eT;

	Channel<eT> C;
	C = identity<Channel<eT>>(0);
	expect_channel(0, 0, C);

	C = identity<Channel<eT>>(3);
	expect_channel("1 0 0; 0 1 0; 0 0 1", C);
}

TYPED_TEST_P(ChannelTest, Randu) {
	typedef TypeParam eT;

	Channel<eT> C(200, 200);
	randu(C);
	expect_channel(200, 200, C);

	C = randu<Channel<eT>>(5);
	expect_channel(5, 5, C);

	C = randu<Channel<eT>>(4, 6);
	expect_channel(4, 6, C);
}



// run the ChannelTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_CASE_P(ChannelTest, Construct, Identity, Randu);

typedef ::testing::Types<double, float, urat> ChannelTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Channel, ChannelTest, ChannelTypes);

