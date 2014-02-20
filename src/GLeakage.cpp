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
GLeakage::GLeakage(chan& channel, Gain& gain_function) {
	if(channel.n_rows != gain_function.n_cols) {
		throw 1;
	}
	C = &channel;
	g = &gain_function;
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
double GLeakage::vulnerability(Prob& pi) {

	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	//the names w x and y are from the formulas.
	double sum_x;
	double max_w = 0;
	for(int w = 0; w < g->n_rows; ++w) {
		sum_x = 0;
		for(int x = 0; x < g->n_cols; ++x) {
			sum_x += pi.at(x) * g->at(w, x);
		}
		if(sum_x > max_w) {
			max_w = sum_x;
		}
	}
	return max_w;
}

double GLeakage::cond_vulnerability(Prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	//the names w x and y are from the formulas.
	double sum_x;
	double max_w;
	double sum_y = 0;

	for(int y = 0; y < C->n_cols; ++y) {
		max_w = 0;
		for(int w = 0; w < g->n_rows; ++w) {
			sum_x = 0;
			for(int x = 0; x < g->n_cols; ++x) {
				sum_x += pi.at(x) * g->at(w, x) * C->at(x, y);
			}
			if(sum_x > max_w) {
				max_w = sum_x;
			}
		}
		sum_y += max_w;
	}
	return sum_y;
}

double GLeakage::leakage(Prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	return log(cond_vulnerability(pi) / vulnerability(pi));
}

double GLeakage::additive_leakage(Prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	return (entropy(pi) - cond_entropy(pi));
}

double GLeakage::entropy(Prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	return -log(vulnerability(pi));
}

double GLeakage::cond_entropy(Prob& pi) {
	if(C->n_rows != pi.size()) {
		throw 1; // X must be equal for both
	}
	return -log(cond_vulnerability(pi));
}

double GLeakage::capacity() {
	throw 1; //It is not supported
}
