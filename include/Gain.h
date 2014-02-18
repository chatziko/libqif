#ifndef _QIF_Gain_h_
#define _QIF_Gain_h_
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
#include "types.h"

/*! \class Gain
 *  \brief A matrix class.
 *
 * Some details about the class
 */
class Gain :
	public arma::mat {

	public:
		// inherit the constructors from parent (C++11 feature)
		using arma::mat::Mat;

		// ! A normal constructor member taking 1 argument.
		/* !
		\param new_gain_elements is an array of double argument wich contains the new elements of the gain matrix ordered by row by row.
		\pre Correct size: the length of each array on new_gain_elements must be the same.
		\sa ~Channel() new_id_function (int size)
		*/
		//Gain (double** new_gain_function );

		//! A normal constructor member taking 1 argument.
		/*!
		\param new_gain_elements an string argument wich contains the new elements of the gain matrix ordered by row by row.
		\pre Correct size: the length of each row on new_gain_elements must be the same.
		\sa ~Channel() new_id_channel (int size)
		*/

		static Gain identity(UIntType size);

};
#endif
