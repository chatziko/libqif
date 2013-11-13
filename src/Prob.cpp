#include "Prob.h"
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

Prob::Prob(StringType& new_vector)
{
	prob_vector=arma::vec(new_vector);
	str=new_vector;
	if(!this->rep_ok()){throw 1;}
}

Prob::Prob(VectorType& new_vector)
{
	prob_vector=arma::vec(new_vector);
	if(!this->rep_ok()){throw 1;}
}

/*
Prob::~Prob()
{
	prob_vector.~vec();
}
*/

int Prob::size()
{
	return prob_vector.size();
}

DoubleType Prob::at(IntType index)
{
	return prob_vector.at(index); //OJO : VER SI ES NECESARIO CONTROLAR EL RANGO.
}