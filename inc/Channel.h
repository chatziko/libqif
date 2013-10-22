#ifndef _QIF_Channel_h_
#define _QIF_Channel_h_

#include <armadillo>
#include <string>
using namespace arma;

/*! \class Channel
 *  \brief A channel matrix class.
 *
 *  A channel must satisfy that the sum each row is 1 and each element is greater than or equal to 0.
 */

class Channel
{
	public:	
		
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
		Channel (std::string& new_channel_elements);
		Channel(mat new_channel_elements);
		
		//! A normal destroyer member.
		/*!
		\sa Channel()
		*/
		~Channel();
		
		//! Generates a new identity channel with the size specified in the argument.
		/*!
		\param size an integer argument which corresponds to the number of rows and columns of the channel.
		\return A new identity channel.
		\sa ~Channel()
		*/
		Channel new_id_channel ( int size );
		
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
		
		//! Checks if the channel is simmetric.
		/*!
		\return Return True iff the channel matrix is simmetric.
		\sa is_partial_simmetric()
		*/
		bool is_symmetric();
		
		//! Checks if the channel is partial simmetric.
		/*!
		\return Return True iff the channel matrix is partial simmetric.
		\sa is_simmetric()
		*/
		bool is_partial_symmetric();
		
		//! Checks if the argument is equal to the current channel.
		/*!
		\param other is a Channel.
		\return Return True iff the argument channel matrix is equal to this.
		*/
		bool is_equal_to ( const Channel& other );
		
		//! Return the row at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the channel.
		\pre The index argument must be between 0 and inputs_number()-1.
		\return The row at the position required.
		\sa inputs_number(),get_column (index)
		*/
		vec get_row ( int index );
		
		//! Return the column at the position the index argument.
		/*!
		\param index is an integer argument which corresponds to the column position on the channel.
		\pre The index argument must be between 0 and outputs_number()-1.
		\return The columun at the position required.
		\sa outputs_number(), get_row (index)
		*/
		vec get_column ( int index );
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the channel.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		void set_row ( int index,vec new_row_elements );
		
		//! Changes the column at the position of the index argument.
		/*!
		\param index is an integer argument which corresponds to the row position on the channel.
		\pre The index argument must be between 0 and inputs_number()-1.
		\pre The number of elements in new_row_elements must be outputs_number().
		*/
		void set_row ( int index,char* new_row_elements );
		
		double at ( int index_x,int index_y );
		
	protected:
		mat C; /*!< This is a matrix defined in the Armadillo Library */
		
		/*! \brief This method checks the invariant representation of the class.
		 *
		 *
		 *  A channel must satisfy that the sum each row is 1 and each element is greater than or equal to 0.
		 */
		bool rep_ok()
		{
			mat::iterator c = C.begin();
			mat::iterator d = C.end();

			int y = C.n_cols - 1; //number of matrix columns
			int j = 0; //index used to check if should change the row.
			double sum = 0; //sumation used to check.
			bool result = true; //flag used to control.

			for ( mat::iterator i = c; i != d && result; ++i )
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

			return result;
		}
	
	
};
#endif
