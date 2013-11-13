#include "Guessing.h"
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
Guessing::Guessing(Channel& channel)
{
	C=&channel;
}
	
//Guessing::~Guessing()
//{
	
//}

//-------------- declaring the theoric algoritmhs implementation
DoubleType Guessing::vulnerability(Prob& pi)
{
	throw 1; //It is not supported
}
		
DoubleType Guessing::cond_vulnerability(Prob& pi)
{
	throw 1; //It is not supported
}		

DoubleType Guessing::leakage(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	return (entropy(pi)-cond_entropy(pi));
}

//internal function to implement entropy and conditional entropy
DoubleType G(Prob& pi){
	DoubleType sum=0;
	for (int x = 0; x < pi.size(); ++x)
	{
		sum += x * pi.at(x);
	}
	return sum;
}

DoubleType Guessing::entropy(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	//Sort pi
	VectorType prob_vector = arma::vec(pi.str);
	sort(prob_vector,1);
	pi = Prob(prob_vector);
	//Call G
	return G(pi);
}		

DoubleType Guessing::cond_entropy(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	DoubleType result=0;
	//for all y
	for(int y=0;y<C->outputs_number();y++){
		//create the vector vy = pi(1)* C[1,y] ... pi(x)* C[x,y]
		VectorType new_vector = arma::vec(C->inputs_number());
		for(int x=0;x<C->inputs_number();x++){
			new_vector[x]= pi.at(x) * C->at(x,y);
		}
		//Sort vy
		sort(new_vector,1); //ascendig
		//call G
		Prob vy= Prob(new_vector);
		result +=G(vy);
	}
	return result;
}		

DoubleType Guessing::capacity()
{
	throw 1; //It is not supported
}
