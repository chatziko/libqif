#ifndef _QIF_Channel_h_
#define _QIF_Channel_h_
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
#include <string>

#include "types.h"
#include "aux.h"

/*! \class Channel
 *  \brief A channel matrix class.
 *
 *  A channel must satisfy that the sum each row is 1 and each element is greater than or equal to 0.
 */
class Channel :
	public arma::mat {

	public:
		// inherit the constructors from parent (C++11 feature)
		using arma::mat::Mat;

		// ! A normal constructor member taking 1 argument.
		/* !
		\param new_channel_elements is an array of double argument wich contains the new elements of the channel ordered by row by row.
		\pre Correct size: the length of each array on new_channel_elements must be the same.
		\pre Correct elements: each element of new_channel_elements must be a number between 0 and 1.
		\pre Probability distribution channel: the sum each outputs_number elements must be 1.
		\sa ~Channel() new_id_channel (int size)
		*/
		//Channel (double** new_channel_elements );

		//! A normal constructor member taking 1 argument.
		/*!
		\param new_channel_elements an array of double argument wich contains the new elements of the channel ordered by row by row.
		\pre Correct size: the length of each row on new_channel_elements must be the same.
		\pre Correct elements: each element of new_channel_elements must be a number between 0 and 1.
		\pre Probability distribution channel: the sum each outputs_number elements must be 1.
		\sa ~Channel() new_id_channel (int size)
		*/
		Channel() {};

		Channel(std::string&);

		Channel(mat&);

		Channel(mat&&);

		//! A normal destroyer member.
		/*!
		\sa Channel()
		*/
		//~Channel();

		//! Generates a new identity channel with the size specified in the argument.
		/*!
		\param size an integer argument which corresponds to the number of rows and columns of the channel.
		\return A new identity channel.
		\sa ~Channel()
		*/
		inline Channel& identity()       { if(!is_square()) throw 1; eye(); return *this; }
		inline Channel& identity(uint n) { eye(n, n); return *this; }

		Channel& randu();

		//! Checks if the channel is symmetric.
		/*!
		\return Return True iff the channel matrix is symmetric.
		\sa is_partial_simmetric()
		*/
		bool is_symmetric();

		bool is_proper();

		bool all(std::function<bool(double)>);
		bool any(std::function<bool(double)>);

		bool is_zero();
};


namespace arma {
	template<>
	struct is_Mat_only< Channel > :
		is_Mat_only< mat > {};

	template<>
	struct is_Mat_only< const Channel > :
		is_Mat_only< const mat > {};

	template<>
	struct is_Mat< Channel > :
		is_Mat_only< mat > {};

	template<>
	struct is_Mat< const Channel > :
		is_Mat_only< const mat > {};

	template<>
	struct Proxy< Channel > :
		Proxy< mat > {
		using Proxy< mat >::Proxy;
	};

	template<>
	struct Proxy< const Channel > :
		Proxy< const mat > {
		using Proxy< const mat >::Proxy;
	};
}

#endif
