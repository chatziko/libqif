#include "Shannon.h"
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
Shannon::Shannon(Channel& channel)
{
	C=&channel;	
}
	
//Shannon::~Shannon()
//{
//	C->~Channel();	
//}

//-------------- declaring the theoric algoritmhs implementation
DoubleType Shannon::vulnerability(Prob& pi)
{
	throw 1; //It is not supported
}
		
DoubleType Shannon::cond_vulnerability(Prob& pi)
{
	throw 1; //It is not supported
}			

DoubleType Shannon::entropy(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	double sum_x=0;
	for(int x=0;x<C->inputs_number();x++){
		sum_x+= pi.at(x) * (log(pi.at(x)) / log(2));
	}
	return -sum_x;
	// - sum x pi(x) log2 p(x)
	//log2 p(x) = log p(x) / log (2)
}		

DoubleType Shannon::cond_entropy(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	double sum_y;
	double sum_x=0;
	for(int x=0; x<C->inputs_number();x++){
		sum_y=0;
		for(int y=0;y<C->outputs_number();y++){
			sum_y+= C->at(x,y) * (log (C->at(x,y)) / log(2));
		}
		sum_x+= pi.at(x) * sum_y;
	}
	return -sum_x;
	//- sum x p(x) * (sum y C[x,y] log2 C[x,y])
	//log2 C[x,y] = log C[x,y] / log (2)
}		

DoubleType Shannon::leakage(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	return(entropy(pi)-cond_entropy(pi));
}	

DoubleType Shannon::capacity()
{
	//implements the Blahut-Arimoto Algorithm
	return 0;
}
