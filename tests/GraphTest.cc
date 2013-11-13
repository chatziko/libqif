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
#include "Graph.h"
#include "gtest/gtest.h"
#include <string>

using namespace std;

TEST(Graph, NoIntegerNodes) {
	string new_graph_elements = "1 2;3 0.5";
	ASSERT_ANY_THROW(Graph graph = Graph(3,new_graph_elements););
}

TEST(Graph, NegativeNodes) {
	string new_graph_elements = "-1 2;3 2";
	ASSERT_ANY_THROW(Graph graph = Graph(3,new_graph_elements););
}

TEST(Graph, NodeZero) {
	string new_graph_elements = "0 2;3 2";
	ASSERT_ANY_THROW(Graph graph = Graph(3,new_graph_elements););
}
 
/*
TEST(Graph, NoCorrectElementsDistribution) {
	string new_graph_elements = "1 2 3;3 2";
	ASSERT_ANY_THROW(Graph graph = Graph(3,new_graph_elements););
}
*/

TEST(Graph, NoCorrectElementsRespectToSize) {
	string new_graph_elements = "1 22;3 2";
	ASSERT_ANY_THROW(Graph graph = Graph(3,new_graph_elements););
}
