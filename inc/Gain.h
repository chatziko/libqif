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
class Gain
{
	public:
		StringType str;
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
		Gain(StringType& new_gain_elements);
		
		//! A normal destroyer member.
		/*!
		\sa Gain()
		*/
		//~Gain();
		
		//! Provides the number of channel columns.
		/*!
		\return The number of channel columns(number of inputs).
		*/
		int inputs_number();
		
		//! Provides the number of channel rows.
		/*!
		\return The number of channel rows(number of guesses).
		*/
		int guesses_number();
		
//		Gain zeros(IntType inputs,IntType guesses);
		
//		Gain ones(IntType inputs,IntType guesses);
		
		Gain new_id_function(IntType size);
		
		//! Create another new gain function with the same content.
		/*!
		\return A new gain function with the same content of this.
		\sa ~Gain()
		*/
		Gain clone();
		
		//! Checks if the gain matrix is simmetric.
		/*!
		\return Return True iff the channel matrix is simmetric.
		\sa is_partial_simmetric()
		*/
		bool is_symmetric();
		
		//! Checks if the gain matrix is partial simmetric.
		/*!
		\return Return True iff the channel matrix is partial simmetric.
		\sa is_simmetric()
		*/
		//bool is_partial_symmetric();
		
		//! Checks if the argument is equal to the current channel.
		/*!
		\param other is a Gain.
		\return Return True iff the argument gain matrix is equal to this.
		*/
		//bool is_equal_to(const Gain& other);
		
		//! Return the row at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the gain matrix.
		\pre The index argument must be between 0 and inputs_number()-1.
		\return The row at the position required.
		\sa inputs_number(),get_column (index)
		*/
		VectorType get_row(IntType index);
		
		//! Return the column at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the column position on the gain matrix.
		\pre The index argument must be between 0 and outputs_number()-1.
		\return The columun at the position required.
		\sa outputs_number(), get_row (index)
		*/
		VectorType get_column(IntType index);
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the gain matrix.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		//void set_row(IntType index,VectorType& new_row_elements);
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the gain matrix.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		//void set_row(IntType index,StringType& new_row_elements);
		
		DoubleType at(IntType index_x,IntType index_y);
		
	protected:
		MatrixType matrix;
		
	private:
		Gain(MatrixType& new_gain_elements);
};
#endif
