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
#include "GLeakage.h"
#include <string> 

int main()
{
    std::string channel_elements ="1 0 0; 0 1 0; 0 0 1";
    Channel C= Channel(channel_elements);

    std::string gain_elements ="1 0 0; 0 1 0; 0 0 1";
    Gain g=Gain(gain_elements);

    GLeakage gl = GLeakage(C,g);

    gl.change_to_scilab();
    gl.plot3d_leakage();

    std::string vector_elements = "0.333 0.333 0.334";
    Prob p1= Prob(vector_elements);

    double lgl=gl.leakage(p1);
}