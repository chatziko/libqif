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
class Channel
{
	public:	
		StringType str;
		
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
		Channel(StringType& new_channel_elements);
		
		Channel(MatrixType& new_channel_elements);
		
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
		Channel new_id_channel(IntType size);
		
		//! Create another new channel with the same content.
		/*!
		\return A new channel with the same content of this.
		\sa ~Channel()
		*/
		Channel clone();
		
		//! Provides the number of channel rows.
		/*!
		\return The number of channel rows(number of inputs).
		*/
		int inputs_number();
		
		//! Provides the number of channel columns.
		/*!
		\return The number of channel columns(number of outputs).
		*/
		int outputs_number();
		
		//! Checks if the channel is simmetrimatrix.
		/*!
		\return Return True iff the channel matrix is simmetrimatrix.
		\sa is_partial_simmetric()
		*/
		bool is_symmetric();
		
		//! Checks if the channel is partial simmetrimatrix.
		/*!
		\return Return True iff the channel matrix is partial simmetrimatrix.
		\sa is_simmetric()
		*/
		//bool is_partial_symmetric();
		
		//! Checks if the argument is equal to the current channel.
		/*!
		\param other is a Channel.
		\return Return True iff the argument channel matrix is equal to this.
		*/
		//bool is_equal_to(const Channel& other);
		
		//! Return the row at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the channel.
		\pre The index argument must be between 0 and inputs_number()-1.
		\return The row at the position required.
		\sa inputs_number(),get_column (index)
		*/
		VectorType get_row(IntType index);
		
		//! Return the column at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the column position on the channel.
		\pre The index argument must be between 0 and outputs_number()-1.
		\return The columun at the position required.
		\sa outputs_number(), get_row (index)
		*/
		VectorType get_column(IntType index);
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the channel.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		//void set_row(IntType index,VectorType& new_row_elements);
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the channel.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		//void set_row(IntType index,StringType& new_row_elements);
		
		DoubleType at(IntType index_x,IntType index_y);
		
	protected:
		MatrixType matrix; /*!< This is a matrix defined in the Armadillo Library */
		/*! \brief This method checks the invariant representation of the class.
		 *
		 *
		 *  A channel must satisfy that the sum each row is 1 and each element is greater than or equal to 0.
		 */
		bool rep_ok()
		{
			int x = matrix.n_rows; //number of matrix rows
			int y = matrix.n_cols; //number of matrix columns
			bool result = true; //flag used to control.
			for (int row = 0; row < x; ++row)
			{
				arma::rowvec current_vector= matrix.row(row);
				arma::vec::iterator c = current_vector.begin();
				arma::vec::iterator d = current_vector.end();

				int j = 0; //index used to check if should change the row.
				double sum = 0; //sumation used to check.
				
				for ( arma::vec::iterator i = c; i != d && result; i++ )
				{
					sum += ( *i );
					result = result && ( *i ) >= 0; // all the elements are greater than or equal to 0.
					j++;

					if ( j == y )
					{
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
