namespace utility {

	template<typename eT>
	eT
	expected_distance(const Mat<eT>& Dist, const Prob<eT>& pi, const Chan<eT>& C) {
		return arma::dot(pi, arma::sum(C % Dist, 1));
	}

	template<typename eT>
	eT
	expected_distance(Metric<eT, uint> dist, const Prob<eT>& pi, const Chan<eT>& C) {
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
