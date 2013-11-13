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
#include "GLeakage.h"

int main()
{
    std::cout << "Using LIBQIF Library Example" << std::endl;
    
    std::string channel_elements = "0.3 0.7; 0.7 0.3; 0.3 0.7";
    Channel C= Channel(channel_elements);
    
    std::string function_elements = "1 0 0; 0 1 0; 0 0 1";
    Gain g=Gain(function_elements);
    
    //Creating the probability distribution vectors
    std::string vector1_elements = "0.33333 0.33333 0.33334";
    std::string vector2_elements = "0 0.5 0.5";
    std::string vector3_elements = "0.5 0.5 0";
    std::string vector4_elements = "0.25 0.5 0.25";    
    
    Prob p1= Prob(vector1_elements);
    Prob p2= Prob(vector2_elements);
    Prob p3= Prob(vector3_elements);    
    Prob p4= Prob(vector4_elements);
    
    //GLeakage
    GLeakage gl= GLeakage(C,g);

    //calculating measures    
    double Lg1=gl.leakage(p1);
    double Lg2=gl.leakage(p2);
    double Lg3=gl.leakage(p3);
    double Lg4=gl.leakage(p4);
    
    std::cout << "Lg p1" << Lg1 << std::endl;
    std::cout << "Lg p2" << Lg2 << std::endl;
    std::cout << "Lg p3" << Lg3 << std::endl;
    std::cout << "Lg p4" << Lg4 << std::endl;
}