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
#include "Mechanism.h"
#include "gtest/gtest.h"
#include <string>

using namespace std;

TEST(Mechanism, CorrectSizeAndElements) {
	string new_channel_elements = "1 0 0;0 1 0";
	string graph_elements = "1 2;1 3;2 3";
    Graph graph = Graph(3,graph_elements);
	Mechanism m=Mechanism(new_channel_elements,graph);
	EXPECT_EQ(2, m.inputs_number());
	EXPECT_EQ(3, m.outputs_number());
}

TEST(Mechanism, NoCorrectElements) {
	string new_channel_elements = "1 2;3 0.5";
	string graph_elements = "1 2; 1 3; 2 3";
    Graph graph = Graph(3,graph_elements);
	ASSERT_ANY_THROW(Mechanism m= Mechanism(new_channel_elements,graph););
}

/* Untested functions:
~Mechanism();
bool is_diffenrential_private(double epsilon);
*/