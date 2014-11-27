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

#include "types.h"
#include "Channel.h"

/*! \class Prob
 *  \brief A probability distribution vector class.
 *
 * This class satisfies that the sum of all the elements is equal to 1 and each element is greater than or equal to 0.
 */
template<typename eT>
class Prob :
	public Row<eT> {

	public:
		// inherit the constructors from parent (C++11 feature)
		using Row<eT>::Row;

		inline Prob() {}

		inline Prob(const char*        s) : Row<eT>(s) { force_proper(); }
		inline Prob(const std::string& s) : Row<eT>(s) { force_proper(); }

		inline Prob(const Prob<eT>& c) : Row<eT>(c) {}
		inline Prob(Prob<eT>&&      c) : Row<eT>(c) {}

		inline Prob(const Row<eT>&     m) : Row<eT>(m) { force_proper(); }
		inline Prob(Row<eT>&&          m) : Row<eT>(m) { force_proper(); }

		inline const Prob& operator=(const Prob&  c) { Row<eT>::operator=(c); return *this; }
		inline const Prob& operator=(const Prob&& c) { Row<eT>::operator=(c); return *this; }

		inline const Prob& operator=(const Row<eT>&  c) { Row<eT>::operator=(c); force_proper(); return *this; }
		inline const Prob& operator=(const Row<eT>&& c) { Row<eT>::operator=(c); force_proper(); return *this; }

		inline const Prob<eT>& uniform() {
			Channel<eT> c(1, this->n_cols);
			c.randu();
			this->row(0) = c.row(0);
			return *this;
		}
		inline const Prob<eT>& uniform(uint s)		{ this->set_size(s); return this->uniform(); }

		bool is_proper() const;
		inline void force_proper() const { if(!is_proper()) throw 1; }
};
#endif
