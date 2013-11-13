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
    std::cout << "Using LIBQIF Library looking for properties" << std::endl;
    
    std::string random = "0.3 0.7; 0.7 0.3; 0.3 0.7";
    std::string balanced = "0.25 0.5 0.25; 0 1 0; 0.5 0 0.5";
    std::string metrics = "1 0.5 0; 0.5 1 0.5; 0 0.5 1";
    std::string id = "1 0 0; 0 1 0; 0 0 1";
    std::string k_tries = "1 1 0; 1 0 1; 0 1 1";
 	   
    Channel C_rand= Channel(random);
    Channel C_balanced= Channel(balanced);
    Channel C_id= Channel(id);
    
    Gain g_id=Gain(id);
    Gain g_metrics=Gain(metrics);
    Gain g_k_tries=Gain(k_tries);
        
    //GLeakage
    GLeakage gl1= GLeakage(C_rand,g_id);
    GLeakage gl2= GLeakage(C_balanced,g_id);
    GLeakage gl3= GLeakage(C_id,g_id);
    
    GLeakage gl4= GLeakage(C_rand,g_metrics);
    GLeakage gl5= GLeakage(C_balanced,g_metrics);
    GLeakage gl6= GLeakage(C_id,g_metrics);
    
    GLeakage gl7= GLeakage(C_rand,g_k_tries);
    GLeakage gl8= GLeakage(C_balanced,g_k_tries);
    GLeakage gl9= GLeakage(C_id,g_k_tries);

    gl1.change_to_scilab();
    gl2.change_to_scilab();
    gl3.change_to_scilab();
    
    gl4.change_to_scilab();
    
    gl5.change_to_scilab();
    gl6.change_to_scilab();
    gl7.change_to_scilab();
    gl8.change_to_scilab();
    gl9.change_to_scilab();

    //ploting

    gl1.plot3d_leakage();
    gl2.plot3d_leakage();
    gl3.plot3d_leakage();
    
    gl4.plot3d_leakage(); //<-----
    
    gl5.plot3d_leakage();
    gl6.plot3d_leakage();
    gl7.plot3d_leakage();
    gl8.plot3d_leakage();
    gl9.plot3d_leakage();
	
}