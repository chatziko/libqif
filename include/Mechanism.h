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
#include "Metric.h"
#include <math.h>

template<typename eT>
class Mechanism {
	public:
		Chan<eT> C;
		Metric<eT, uint> d = {metrics::euclidean<eT, uint>()};

//		Mechanism<eT>()																{};
//		Mechanism<eT>(const Chan<eT>& C)							: C(C)			{};
//		Mechanism<eT>(const Chan<eT>& C, const Metric<eT, uint> d)	: C(C), d(d)	{};

		bool is_private(eT epsilon = eT(1));
		eT smallest_epsilon();
};

namespace mechanisms {

	template<typename eT>
	Mechanism<eT>
	geometric(uint n, eT step = eT(1), eT epsilon = eT(1)) {
		Mechanism<eT> mech;

		eT c = std::exp(step * epsilon),
		   lambda_f = c / (c + eT(1)),				// lambda for the first/last cols
		   lambda_m = (c - eT(1)) / (c + eT(1));	// lambda for the middle colums

		mech.d = metrics::scale(metrics::euclidean<eT, uint>(), step);
		mech.C.set_size(n, n);

		for(uint j = 0; j < n; j++) {
			eT lambda = (j == 0 || j == n-1 ? lambda_f : lambda_m);
			for(uint i = 0; i < n; i++)
				mech.C(i, j) = lambda * std::exp(- epsilon * mech.d(i, j));
		}

		return mech;
	}

	template<typename eT>
	Mechanism<eT>
	exponential(uint n, Metric<eT, uint> d, eT epsilon = eT(1)) {
		Mechanism<eT> mech;
		mech.d = d;

		auto& C = mech.C;
		C.set_size(n, n);
		C.diag().fill(eT(1));

		for(uint i = 0; i < n; i++) {
			for(uint j = i + 1; j < n; j++)
				C(i, j) = C(j, i) = std::exp(- epsilon/2 * d(i, j));
			C.row(i) /= accu(C.row(i));
		}

		return mech;
	}

	template<typename eT>
	Mechanism<eT>
	tight_constraints(uint n, Metric<eT, uint> d, eT epsilon = eT(1)) {
		Mechanism<eT> mech;
		mech.d = d;

		// build phi, it will be later transformed to a channel matrix
		auto& phi = mech.C;
		phi.set_size(n, n);
		phi.diag().fill(eT(1));

		for(uint i = 0; i < n; i++)
			for(uint j = i + 1; j < n; j++)
				phi(i, j) = phi(j, i) = std::exp(- epsilon * d(i, j));

		// invert and create diagonal
		//
		// Note: in the perl code, scaling the ones(n) vector somehow improves numerical stability
		//       (the non-negativity of diag). We should investigate if this still happens
		// eT scaler = n;
		// Row<eT> diag = scaler * (1/scaler * arma::ones<Row<eT>>(n) * phi.i());
		//
		Row<eT> diag = arma::ones<Row<eT>>(n) * phi.i();

		for(uint i = 0; i < n; i++)
			if(!less_than_or_eq(eT(0), diag(i)))
				throw std::runtime_error("negative");

		// the channel matrix = phi after multiplying all rows with diag
		phi.each_row() %= diag;

		return mech;
	}
}

#endif
