#ifndef _QIF_Graph_h_
#define _QIF_Graph_h_
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
#include <armadillo>
#include <vector>
#include <utility>
#include "types.h"

class Graph {
	public:
		Graph(uint vertex_num, std::string& edges);

		Graph(uint vertex_num, std::vector< std::pair<int, int> >& edges);

//		~Graph();

		uint vertex_number();
		bool is_an_edge(uint v1, uint v2);
		uint get_distance(uint v1, uint v2);

	protected:
		uint V;
		mat adjacency;
		mat distances;
};
#endif
