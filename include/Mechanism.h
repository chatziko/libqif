#ifndef _QIF_Mechanism_h_
#define _QIF_Mechanism_h_
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
#include "Chan.h"
#include "Graph.h"
#include "Metric.h"
#include <math.h>

template<typename eT>
class Mechanism {
	public:
		Chan<eT> C;
		const Metric<eT, uint> d = {metrics::euclidean<eT, uint>()};

		Mechanism<eT>()																{};
		Mechanism<eT>(const Chan<eT>& C)							: C(C)			{};
		Mechanism<eT>(const Chan<eT>& C, const Metric<eT, uint> d)	: C(C), d(d)	{};

		bool is_private(eT epsilon = eT(1));
		eT smallest_epsilon();
};

#endif
