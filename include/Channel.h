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
#include <string>

#include "types.h"
#include "aux.h"

/*! \class Channel
 *  \brief A channel matrix class.
 *
 *  A channel must satisfy that the sum each row is 1 and each element is greater than or equal to 0.
 */

template<typename eT>
class Channel :
	public Mat<eT> {

	public:
		// inherit the constructors from parent (C++11 feature)
		using Mat<eT>::Mat;

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
		inline Channel() {}

		inline Channel(const char*        s) : Mat<eT>(s) { force_proper(); }
		inline Channel(const std::string& s) : Mat<eT>(s) { force_proper(); }

		inline Channel(const Channel<eT>& c) : Mat<eT>(c) {}
		inline Channel(Channel<eT>&&      c) : Mat<eT>(c) {}

		inline Channel(const Mat<eT>&     m) : Mat<eT>(m) { force_proper(); }
		inline Channel(Mat<eT>&&          m) : Mat<eT>(m) { force_proper(); }

		inline const Channel& operator=(const Channel&  c) { Mat<eT>::operator=(c); return *this; }
		inline const Channel& operator=(const Channel&& c) { Mat<eT>::operator=(c); return *this; }

		inline const Channel& operator=(const Mat<eT>&  c) { Mat<eT>::operator=(c); force_proper(); return *this; }
		inline const Channel& operator=(const Mat<eT>&& c) { Mat<eT>::operator=(c); force_proper(); return *this; }

		//! A normal destroyer member.
		/*!
		\sa Channel()
		*/
		//~Channel() { std::cout << "destroy\n"; }

		//! Generates a new identity channel with the size specified in the argument.
		/*!
		\param size an integer argument which corresponds to the number of rows and columns of the channel.
		\return A new identity channel.
		\sa ~Channel()
		*/
		// Note: inherited methods must be called as this->method since no declaration is available
		//
		inline const Channel<eT>& identity()			{ if(!this->is_square()) throw 1; this->eye(); return *this; }
		inline const Channel<eT>& identity(uint n)		{ this->eye(n, n); return *this; }

		inline const Channel<eT>& randu()				{ Mat<eT>::randu();     normalize(); return *this; }
		inline const Channel<eT>& randu(uint s)			{ Mat<eT>::randu(s);    normalize(); return *this; }
		inline const Channel<eT>& randu(uint r, uint c)	{ Mat<eT>::randu(r, c); normalize(); return *this; }

		inline const Channel<eT>& normalize()			{ this->each_col() /= sum(*this, 1); return *this; }

		//! Checks if the channel is symmetric.
		/*!
		\return Return True iff the channel matrix is symmetric.
		\sa is_partial_simmetric()
		*/
		bool is_symmetric() const;

		bool is_proper() const;
		inline void force_proper() const { if(!is_proper()) throw 1; }

		bool all(std::function<bool(const eT&)>) const;
		bool any(std::function<bool(const eT&)>) const;

		bool is_zero() const;
};


namespace arma {
	// make armadillo support boost::rational<eT> as element type
	template<typename eT>
	struct is_supported_elem_type< boost::rational<eT> > {
		static const bool value = true;
	};

	template<typename eT>
	struct is_Mat_only< Channel<eT> > :
		is_Mat_only< Mat<eT> > {};

	template<typename eT>
	struct is_Mat_only< const Channel<eT> > :
		is_Mat_only< const Mat<eT> > {};

	template<typename eT>
	struct is_Mat< Channel<eT> > :
		is_Mat_only< Mat<eT> > {};

	template<typename eT>
	struct is_Mat< const Channel<eT> > :
		is_Mat_only< const Mat<eT> > {};

	template<typename eT>
	struct Proxy< Channel<eT> > :
		Proxy< Mat<eT> > {
		using Proxy< Mat<eT> >::Proxy;
	};

	template<typename eT>
	struct Proxy< const Channel<eT> > :
		Proxy< const Mat<eT> > {
		using Proxy< const Mat<eT> >::Proxy;
	};
}

#endif
