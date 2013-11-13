#ifndef _QIF_Prob_h_
#define _QIF_Prob_h_
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
/*! \class Prob
 *  \brief A probability distribution vector class.
 *
 * This class satisfies that the sum of all the elements is equal to 1 and each element is greater than or equal to 0.
 */
class Prob
{
	public:
		StringType str;
		//! A normal constructor member taking 1 argument.
		/*!
		\param elements_number is an integer which represents the number of elements of the new vector.
		\param new_vector is a pointer to a doubles array which the elements of the new probability vector.
		\pre Correct size: the elements_number argument must coincide with the length of the array new_vector_elements.
		\pre Correct elements: each element of new_vector_elements must be a number between 0 and 1.
		\pre Probability distribution vector: the sum of all the elements must be 1.
		\sa ~Probability_vector()
		*/
		//Prob (double* new_vector_elements );

		//! A normal constructor member taking 1 argument.
		/*!
		\param elements_number is an integer which represents the number of elements of the new vector.
		\param new_vector is a pointer to a doubles array which the elements of the new probability vector.
		\pre Correct size: the elements_number argument must coincide with the length of the array new_vector_elements.
		\pre Correct elements: each element of new_vector_elements must be a number between 0 and 1.
		\pre Probability distribution vector: the sum of all the elements must be 1.
		\sa ~Probability_vector()
		*/
		Prob(StringType& new_vector_elements);

		Prob(VectorType& vector_elements);
		
		//! normal destroyer member.
		/*!
		\sa Probability_vector()
		*/
		//~Prob ();
		
		//! A function which returns the length of the vector.
		/*!
		\return The test results
		*/
		IntType size();
		
		//! A function wich takes an index position and returns the choosen element.
		/*!
		\param inputs an integer argument which corresponds to the choosen position.
		\pre The index argument must be a number between 0 and size()-1.
		\return The element at position index.
		\sa size()
		*/
		DoubleType at(IntType index);
		
	protected:
		VectorType prob_vector;/*!< This is a vector defined in the Armadillo Library */

		/*! \brief This method checks the invariant representation of the class.
		 *
		 *
		 *  A probability vector must satisfy that the sum of all elements is 1 and each element is greater than or equal to 0.
		 */
		bool rep_ok()
		{
			arma::vec::iterator c = prob_vector.begin();
			arma::vec::iterator d = prob_vector.end();

			double sum = 0; //sumation used to check.
			bool result = true; //flag used to control.

			for ( arma::vec::iterator i = c; i != d && result; ++i )
			{
				sum += ( *i );
				result = result && ( *i ) >= 0; // all the elements are greater than or equal to 0.
			}
			result = result && sum == 1; //the sum of each row is 1.
			return result;
		}
};
#endif
