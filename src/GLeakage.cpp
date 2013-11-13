#include "GLeakage.h"
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
GLeakage::GLeakage(Channel& channel,Gain& gain_function)
{
	if(channel.inputs_number()!=gain_function.inputs_number())
	{
		throw 1; // X must be equal for both
	}
	C=&channel;
	g=&gain_function;
}	

//GLeakage::~GLeakage()
//{
//	C->~Channel();
//	g->~Gain();
//}	
/*
void * GLeakage::compare_over_prior(Channel& other_channel)
{
	
}
	
void * GLeakage::compare_over_gain(Channel& other_channel,Prob& prior)
{
	
}
*/
//-------------- declaring the theoric algoritmhs implementation
DoubleType GLeakage::vulnerability(Prob& pi)
{

	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	//the names w x and y are from the formulas.
	double sum_x;
	double max_w=0;
	for(int w=0;w<g->guesses_number();++w)
	{
		sum_x=0;
		for(int x=0;x<g->inputs_number();++x)
		{
			sum_x+=pi.at(x) * g->at(w,x);
		}
		if(sum_x>max_w){
			max_w=sum_x;
		}
	}
	return max_w;
}
		
DoubleType GLeakage::cond_vulnerability(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	//the names w x and y are from the formulas.
	double sum_x;
	double max_w;
	double sum_y=0;
		
	for(int y=0;y<C->outputs_number();++y)
	{
		max_w=0;
		for(int w=0;w<g->guesses_number();++w)
		{
			sum_x=0;
			for(int x=0;x<g->inputs_number();++x)
			{
				sum_x+=pi.at(x) * g->at(w,x) * C->at(x,y);
			}
			if(sum_x>max_w){
				max_w=sum_x;
			}
		}
		sum_y+=max_w;
	}
	return sum_y;
}	
	
DoubleType GLeakage::leakage(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	return log (cond_vulnerability(pi)/vulnerability(pi));
}	

DoubleType GLeakage::additive_leakage(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	return (entropy(pi)-cond_entropy(pi));
}
	
DoubleType GLeakage::entropy(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	return -log(vulnerability(pi));
}		

DoubleType GLeakage::cond_entropy(Prob& pi)
{
	if(C->inputs_number()!=pi.size())
	{
		throw 1; // X must be equal for both
	}
	return -log(cond_vulnerability(pi));
}		

DoubleType GLeakage::capacity()
{
	throw 1; //It is not supported
}