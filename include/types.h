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

using arma::Mat;
using arma::mat;
using arma::vec;
using arma::Row;

// helper aliases for SFINAE, see http://loungecpp.wikidot.com/tips-and-tricks:enable-if-for-c-11
enum class enabled {}; // just a type that can be used as a template parameter and is as inocuous as possible
template <typename T> using Invoke = typename T::type;
template <typename Condition> using EnableIf  = Invoke<std::enable_if< Condition::value, enabled>>;
template <typename Condition> using DisableIf = Invoke<std::enable_if<!Condition::value, enabled>>;


template<typename IntType> class rational; // forward

template<typename eT> using Chan    = Mat<eT>;
template<typename T>  using is_Chan = arma::is_Mat_only<T>;

template<typename eT> using Prob    = Row<eT>;
template<typename T>  using is_Prob = arma::is_Row<T>;

typedef uint32_t uint;

template<typename IntType> using Rational = rational<IntType>;
typedef Rational<uintmax_t> urat;

typedef Chan<double>  chan;
typedef Chan<float>  fchan;
typedef Chan<urat>   rchan;

typedef Row<double>  prob;
typedef Row<float>  fprob;
typedef Row<urat>   rprob;


template<typename T>
struct is_Rational
  { static const bool value = false; };
template<typename IntType>
struct is_Rational< Rational<IntType> >
  { static const bool value = true; };


#endif
