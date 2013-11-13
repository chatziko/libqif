#include "Mechanism.h"
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
Mechanism::Mechanism(StringType& new_channel_elements,Graph& new_graph): Channel(new_channel_elements)
{
	graph = &new_graph;
}

//Mechanism::~Mechanism()
//{
	//Its not implemented yet
//}

bool Mechanism::is_differential_private(DoubleType epsilon)
{
/* Algorithm
for each (x, x') in edges(G):
    for each y in Y:
        if C[x,y] > e^epsilon * C[x',y]:
            return false
return true
*/
//for each (x, x') in edges(G):
	for (int x = 0;x<graph->vertex_number();++x){
		for (int x2 = 0; x2 < graph->vertex_number(); ++x2)
		{
			if(graph->is_an_edge(x,x2))
			{
				//for each y in Y:
				for (int y = 0; y < matrix.n_cols; ++y)
				{
					//if C[x,y] > e^epsilon * C[x',y]:
					if (matrix.at(x,y) > exp(epsilon) * matrix.at(x2,y))
					{
						return false;
					}
				}
			}
		}
	}
	return true;
}