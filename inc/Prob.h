#ifndef _QIF_Prob_h_
#define _QIF_Prob_h_

#include <armadillo>
using namespace arma;

/*! \class Prob
 *  \brief A probability distribution vector class.
 *
 * This class satisfies that the sum of all the elements is equal to 1 and each element is greater than or equal to 0.
 */
class Prob
{
	public:
		//! A normal constructor member taking 1 argument.
		/*!
		\param elements_number is an integer which represents the number of elements of the new vector.
		\param new_vector is a pointer to a doubles array which the elements of the new probability vector.
		\pre Correct size: the elements_number argument must coincide with the length of the array new_vector_elements.
		\pre Correct elements: each element of new_vector_elements must be a number between 0 and 1.
		\pre Probability distribution vector: the sum of all the elements must be 1.
		\sa ~Probability_vector()
		*/
		Prob (double* new_vector_elements );

		//! A normal constructor member taking 1 argument.
		/*!
		\param elements_number is an integer which represents the number of elements of the new vector.
		\param new_vector is a pointer to a doubles array which the elements of the new probability vector.
		\pre Correct size: the elements_number argument must coincide with the length of the array new_vector_elements.
		\pre Correct elements: each element of new_vector_elements must be a number between 0 and 1.
		\pre Probability distribution vector: the sum of all the elements must be 1.
		\sa ~Probability_vector()
		*/
		Prob (char* new_vector_elements );
		
		//! normal destroyer member.
		/*!
		\sa Probability_vector()
		*/
		~Prob ();
		
		//! A function which returns the length of the vector.
		/*!
		\return The test results
		*/
		int size();
		
		//! A function wich takes an index position and returns the choosen element.
		/*!
		\param inputs an integer argument which corresponds to the choosen position.
		\pre The index argument must be a number between 0 and size()-1.
		\return The element at position index.
		\sa size()
		*/
		double at ( int index );
		
	protected:
		vec V;/*!< This is a vector defined in the Armadillo Library */

		/*! \brief This method checks the invariant representation of the class.
		 *
		 *
		 *  A probability vector must satisfy that the sum of all elements is 1 and each element is greater than or equal to 0.
		 */
		bool rep_ok()
		{
			vec::iterator c = V.begin();
			vec::iterator d = V.end();

			double sum = 0; //sumation used to check.
			bool result = true; //flag used to control.

			for ( mat::iterator i = c; i != d && result; ++i )
			{
				sum += ( *i );
				result = result && ( *i ) >= 0; // all the elements are greater than or equal to 0.
			}
			result = result && sum == 1; //the sum of each row is 1.
			return result;
		}
};
#endif
