
namespace mechanism {

template<typename eT>
Mech<eT>
geometric(uint n, eT step = eT(1), eT epsilon = eT(1)) {
	Mech<eT> mech;

	eT c = std::exp(step * epsilon),
	   lambda_f = c / (c + eT(1)),				// lambda for the first/last cols
	   lambda_m = (c - eT(1)) / (c + eT(1));	// lambda for the middle colums

	mech.d = step * metric::euclidean<eT, uint>();
	mech.C.set_size(n, n);

	for(uint j = 0; j < n; j++) {
		eT lambda = (j == 0 || j == n-1 ? lambda_f : lambda_m);
		for(uint i = 0; i < n; i++)
			mech.C(i, j) = lambda * std::exp(- epsilon * mech.d(i, j));
	}

	return mech;
}

template<typename eT>
Mech<eT>
exponential(uint n, Metric<eT, uint> d, eT epsilon = eT(1)) {
	Mech<eT> mech;
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
Mech<eT>
tight_constraints(uint n, Metric<eT, uint> d, eT epsilon = eT(1)) {
	Mech<eT> mech;
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

template<typename eT>
bool is_private(const Mech<eT>& m, eT epsilon) {
	auto mtv = metric::mult_total_variation<eT, Prob<eT>>();

	for(uint i = 0; i < m.C.n_rows; i++) {
		for(uint j = i+1; j < m.C.n_rows; j++) {
			eT mp = mtv(m.C.row(i), m.C.row(j));

			if(!less_than_or_eq(mp, epsilon * m.d(i, j)))
				return false;
		}
	}
	return true;
}

template<typename eT>
eT smallest_epsilon(const Mech<eT>& m) {
	auto mtv = metric::mult_total_variation<eT, Prob<eT>>();

	eT res(0);
	for(uint i = 0; i < m.C.n_rows; i++) {
		for(uint j = i+1; j < m.C.n_rows; j++) {
			eT ratio = mtv(m.C.row(i), m.C.row(j)) / m.d(i, j);
			if(less_than(res, ratio))
				res = ratio;
		}
	}
	return res;
}

} // namespace mechanism
