
namespace mechanism {

template<typename eT>
Chan<eT>
geometric(uint n, Metric<eT, uint> d = metric::euclidean<eT,uint>()) {
	Chan<eT> C(n, n);

	eT c = std::exp(d(0,1)),
	   lambda_f = c / (c + eT(1)),				// lambda for the first/last cols
	   lambda_m = (c - eT(1)) / (c + eT(1));	// lambda for the middle colums

	for(uint j = 0; j < n; j++) {
		eT lambda = (j == 0 || j == n-1 ? lambda_f : lambda_m);
		for(uint i = 0; i < n; i++)
			C(i, j) = lambda * std::exp(-d(i, j));
	}

	return C;
}

template<typename eT>
Chan<eT>
exponential(uint n, Metric<eT, uint> d) {
	Chan<eT> C(n, n);
	C.diag().fill(eT(1));

	for(uint i = 0; i < n; i++) {
		for(uint j = i + 1; j < n; j++)
			C(i, j) = C(j, i) = std::exp(-d(i, j)/2);
		C.row(i) /= accu(C.row(i));
	}

	return C;
}

template<typename eT>
Chan<eT>
tight_constraints(uint n, Metric<eT, uint> d) {
	// build phi, it will be later transformed to a channel matrix
	Chan<eT> phi(n, n);
	phi.diag().fill(eT(1));

	for(uint i = 0; i < n; i++)
		for(uint j = i + 1; j < n; j++)
			phi(i, j) = phi(j, i) = std::exp(-d(i, j));

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

	return phi;
}

template<typename eT>
bool is_private(const Chan<eT>& C, Metric<eT, uint> d) {
	auto mtv = metric::mult_total_variation<eT, Prob<eT>>();

	for(uint i = 0; i < C.n_rows; i++) {
		for(uint j = i+1; j < C.n_rows; j++) {
			// non-adjacent elements are redundant to check
			if(!d.is_adjacent(i, j)) continue;

			eT mp = mtv(C.row(i), C.row(j));

			if(!less_than_or_eq(mp, d(i, j)))
				return false;
		}
	}
	return true;
}

template<typename eT>
eT smallest_epsilon(const Chan<eT>& C, Metric<eT, uint> d) {
	auto mtv = metric::mult_total_variation<eT, Prob<eT>>();

	eT res(0);
	for(uint i = 0; i < C.n_rows; i++) {
		for(uint j = i+1; j < C.n_rows; j++) {
			// non-adjacent elements are redundant to check
			if(!d.is_adjacent(i, j)) continue;

			eT ratio = mtv(C.row(i), C.row(j)) / d(i, j);
			if(less_than(res, ratio))
				res = ratio;
		}
	}
	return res;
}

} // namespace mechanism
