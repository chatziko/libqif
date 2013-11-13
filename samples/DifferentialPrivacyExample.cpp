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
#include <iostream>
#include <string>

#include "Mechanism.h"

int main()
{
	std::cout << "Using QIF Library Example" << std::endl;
   
	//the following line is an example of how to write the graph edges   
    std::string graph_elements = "1 2; 1 3; 2 3";
    //this graph have 3 vertex.
    Graph graph = Graph(3,graph_elements);

    std::string channel_elements ="1 0 0; 0 1 0; 0 0 1";
    Mechanism mechanism= Mechanism(channel_elements,graph);

    //this example asks with epsilon=0.05
    bool result= mechanism.is_differential_private(0.05);
}