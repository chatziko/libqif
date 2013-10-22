#ifndef _QIF_Gain_h_
#define _QIF_Gain_h_

#include <armadillo>
using namespace arma;

/*! \class Gain
 *  \brief A matrix class.
 *
 * Some details about the class
 */

class Gain
{
	public:
		
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
		Gain (char* new_gain_elements );
		
		//! A normal destroyer member.
		/*!
		\sa Gain()
		*/
		~Gain();
		
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
		
		Gain zeros ( int inputs, int guesses );
		
		Gain ones ( int inputs, int guesses );
		
		Gain new_id_function ( int size);
		
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
		bool is_partial_symmetric();
		
		//! Checks if the argument is equal to the current channel.
		/*!
		\param other is a Gain.
		\return Return True iff the argument gain matrix is equal to this.
		*/
		bool is_equal_to ( const Gain& other );
		
		//! Return the row at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the gain matrix.
		\pre The index argument must be between 0 and inputs_number()-1.
		\return The row at the position required.
		\sa inputs_number(),get_column (index)
		*/
		vec get_row ( int index );
		
		//! Return the column at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the column position on the gain matrix.
		\pre The index argument must be between 0 and outputs_number()-1.
		\return The columun at the position required.
		\sa outputs_number(), get_row (index)
		*/
		vec get_column ( int index );
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the gain matrix.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		void set_row ( int index,vec new_row_elements );
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the gain matrix.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		void set_row ( int index,char* new_row_elements );
		
		double at ( int index_x,int index_y );
		
	protected:
		mat C;
		
	private:
		Gain (mat new_gain_elements );
};
#endif
