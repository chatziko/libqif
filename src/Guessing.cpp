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
Guessing::Guessing(chan& channel) {
	C = &channel;
}

//Guessing::~Guessing()
//{

//}

//-------------- declaring the theoric algoritmhs implementation
double Guessing::vulnerability(prob& pi) {
	throw 1; //It is not supported
}

double Guessing::cond_vulnerability(prob& pi) {
	throw 1; //It is not supported
}

double Guessing::leakage(prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	return (entropy(pi) - cond_entropy(pi));
}

//internal function to implement entropy and conditional entropy
double G(prob& pi) {
	//TODO: sort doesn't work cause armadillo types are messed
	//sort(pi);

	double sum = 0;
	for(uint x = 0; x < pi.size(); ++x) {
		sum += x * pi.at(x);
	}
	return sum;
}

double Guessing::entropy(prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	//Call G
	prob clone = pi;
	return G(clone);
}

double Guessing::cond_entropy(prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	double result = 0;
	//for all y
	for(uint y = 0; y < C->n_cols; y++) {
		//create the vector vy = pi(1)* C[1,y] ... pi(x)* C[x,y]
		prob vy(C->n_rows);
		for(uint x = 0; x < C->n_rows; x++) {
			vy.at(x) = pi.at(x) * C->at(x, y);
		}
		//call G
		result += G(vy);
	}
	return result;
}

double Guessing::capacity() {
	throw 1; //It is not supported
}
