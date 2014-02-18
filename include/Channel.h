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
		Channel();

		Channel(StringType&);

		Channel(MatrixType&);

		Channel(MatrixType&&);

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
		static Channel identity(UIntType size);

		//! Checks if the channel is symmetric.
		/*!
		\return Return True iff the channel matrix is symmetric.
		\sa is_partial_simmetric()
		*/
		bool is_symmetric();

	protected:
		/*! \brief This method checks the invariant representation of the class.
		 *
		 *
		 *  A channel must satisfy that the sum each row is 1 and each element is greater than or equal to 0.
		 */
		bool rep_ok() {
			int x = this->n_rows;
			int y = this->n_cols;
			bool result = true; //flag used to control.
			for(int row = 0; row < x; ++row) {
				arma::rowvec current_vector = this->row(row);
				arma::vec::iterator c = current_vector.begin();
				arma::vec::iterator d = current_vector.end();

				int j = 0; //index used to check if should change the row.
				double sum = 0; //sumation used to check.

				for(arma::vec::iterator i = c; i != d && result; i++) {
					sum += (*i);
					result = result && (*i) >= 0;   // all the elements are greater than or equal to 0.
					j++;

					if(j == y) {
						result = result && sum == 1; //the sum of each row is 1.
						j = 0;
						sum = 0;
					}
				}
			}
			return result;
		}
};
#endif
