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

#include "Channel.h"
#include "aux.h"


Channel::Channel(std::string& s) :
	mat(s) {

	if(!this->is_proper())
		throw 1;
}

Channel::Channel(mat& m) :
	mat(m) {

	if(!this->is_proper())
		throw 1;
}

Channel::Channel(mat&& m) :
	mat(m) {

	if(!this->is_proper())
		throw 1;
}

Channel& Channel::randu() {
	mat::randu();

	//cout <<  arma::sum(*this, 1);
	//each_col() /= sums;

	return *this;
}

bool Channel::is_symmetric() {
	if(this->n_cols != this->n_rows) return false;

	uint i, j;
	for(i = 0; i < this->n_cols; i++)
		for(j = i + 1; j < this->n_rows; j++)
			if(this->at(i, j) != this->at(j, i)) return false;

	return true;
}

bool Channel::is_proper() {
	for(uint i = 0; i < n_rows; i++) {
		double sum = 0;
		for(uint j = 0; j < n_cols; j++) {
			// elements should be non-negative
			double elem = at(i, j);
			if(less_than(elem, 0))
				return false;

			sum += elem;
		}

		// rows should add up to 1
		if(!equal(sum, 1))
			return false;
	}
	return true;
}

bool Channel::all(std::function<bool(double)> f) {
	for(double x : *this)
		if(!f(x))
			return false;
	return true;
}

bool Channel::any(std::function<bool(double)> f) {
	for(double x : *this)
		if(f(x))
			return true;
	return false;
}

bool Channel::is_zero() {
	return all([](double x){
		return equal(x, 0);
	});
}

