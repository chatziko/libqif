#ifndef _QIF_Types_h_
#define _QIF_Types_h_
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
#include <armadillo>

#include "rat.h"
#include "aux.h"


using arma::Mat;
using arma::mat;
using arma::Row;
using arma::Col;

// helper aliases for SFINAE, see http://loungecpp.wikidot.com/tips-and-tricks:enable-if-for-c-11
enum class enabled {}; // just a type that can be used as a template parameter and is as inocuous as possible
template <typename T> using Invoke = typename T::type;
template <typename Condition> using EnableIf  = Invoke<std::enable_if< Condition::value, enabled>>;

template <typename T> using eT = typename T::elem_type;

template<typename eT> using Chan    = Mat<eT>;
template<typename T>  using is_Chan = arma::is_Mat_only<T>;

template<typename eT> using Prob    = Row<eT>;
template<typename T>  using is_Prob = arma::is_Row<T>;

typedef uint32_t uint;

typedef Chan<double> chan;
typedef Chan<float> fchan;
typedef Chan<rat>   rchan;

typedef Row<double>  prob;
typedef Row<float>  fprob;
typedef Row<rat>    rprob;

typedef Mat<rat>     rmat;
typedef Col<rat>  rcolvec;
typedef Row<rat>  rrowvec;


template<typename eT>
struct Point {
	eT x, y;
	Point() {}
	Point(eT x, eT y) : x(x), y(y) {}

	bool operator==(const Point<eT>& rhs) const { return equal(this->x, rhs.x) && equal(this->y, rhs.y); }
};

template<typename T>
struct is_Point { static const bool value = false; };
template<typename eT>
struct is_Point<Point<eT>> { typedef eT elem_type; static const bool value = true;  };

typedef Point<double> point;
typedef Point<float> fpoint;
typedef Point<rat>   rpoint;


#endif
