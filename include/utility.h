#ifndef _QIF_utility_h_
#define _QIF_utility_h_

#include "Chan.h"
#include "Metric.h"


namespace utility {

	template<typename eT>
	eT
	expected_distance(const Chan<eT>& C, const Prob<eT>& pi, const Mat<eT>& Dist) {
		return arma::dot(pi, arma::sum(C % Dist, 1));
	}

	template<typename eT>
	eT
	expected_distance(const Chan<eT>& C, const Prob<eT>& pi, Metric<eT, uint> dist) {
		eT sum(0);
		for(uint i = 0; i < C.n_rows; i++) {
			eT sum2(0);
			for(uint j = 0; j < C.n_cols; j++)
				sum2 += C(i,j) * dist(i, j);
			sum += pi(i) * sum2;
		}
		return sum;
	}
}

#endif
