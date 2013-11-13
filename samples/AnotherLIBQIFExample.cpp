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

#include "Shannon.h"
#include "MinEntropy.h"
#include "Guessing.h"
#include "GLeakage.h"
int main()
{
	
    std::cout << "Using QIF Library Example" << std::endl;
       
    //Creating the channel matrix
    std::string channel_elements ="1 0 0; 0 1 0; 0 0 1";
    Channel C= Channel(channel_elements);
    std::cout << C.str << std::endl;

    //Creating the gain function matrix reusing the channel elements
    Gain g=Gain(channel_elements);
    
    //Creating the probability distribution
    std::string vector_elements ="0.3333 0.3333 0.3334";
    Prob p1= Prob(vector_elements);

    std::cout << "Calculating the GLeakage" << std::endl;
    //Calculating the GLeakage
    GLeakage gl= GLeakage(C,g);
    double Lg=gl.leakage(p1);
    std::cout << "Lg " << Lg << std::endl;
    std::cout << "Calculating ends" << std::endl;

    Shannon sl = Shannon(C);
    MinEntropy ml = MinEntropy(C);
    Guessing gul = Guessing(C);
    //Using the Plotter 
    gl.change_to_scilab();
    gl.plot3d_leakage();

    sl.change_to_scilab();
    sl.plot3d_entropy();
}
