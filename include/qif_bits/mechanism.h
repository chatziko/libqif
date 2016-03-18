
namespace mechanism {

// build a matrix D where Dij = exp(-d(i,j))
// used in the construction of most mechanisms below
//
template<typename eT>
Mat<eT>
distance_matrix(uint n_rows, uint n_cols, Metric<eT, uint> d) {
	Mat<eT> D(n_rows, n_cols);
	D.diag().fill(eT(1));

	// exploit symmetry, but careful about non-square matrices
	//
	for(uint i = 0; i < n_rows; i++)
		for(uint j = (i < n_cols ? i + 1 : 0); j < n_cols; j++) {
			D(i, j) = std::exp(-d(i, j));
			if(j < n_rows && i < n_cols)
				D(j, i) = D(i, j);
		}
	return D;
}

template<typename eT>
Chan<eT>
geometric(uint n_rows, Metric<eT, uint> d = metric::euclidean<eT,uint>(), uint n_cols = 0) {
	if(n_cols == 0) n_cols = n_rows;
	if(n_rows > n_cols) throw std::runtime_error("n_cols should be at least as big as n_rows");

	eT c = std::exp(d(0,1)),
	   lambda_f = c / (c + eT(1)),				// lambda for the first/last cols
	   lambda_m = (c - eT(1)) / (c + eT(1));	// lambda for the middle colums

	Chan<eT> C = distance_matrix(n_rows, n_cols, d);
	C.col(0)            *= lambda_f;
	C.col(n_cols-1)     *= lambda_f;
	C.cols(1, n_cols-2) *= lambda_m;

	return C;
}

template<typename eT>
Chan<eT>
exponential(uint n_rows, Metric<eT, uint> d, uint n_cols = 0) {
	if(n_cols == 0) n_cols = n_rows;

	Chan<eT> C = distance_matrix(n_rows, n_cols, eT(1)/2 * d);

	for(uint i = 0; i < n_rows; i++)
		C.row(i) /= accu(C.row(i));
	return C;
}

template<typename eT>
Chan<eT>
tight_constraints(uint n, Metric<eT, uint> d) {
	// invert distance matrix phi and create diagonal
	//
	// Note: in the perl code, scaling the ones(n) vector somehow improves numerical stability
	//       (the non-negativity of diag). We should investigate if this still happens
	// eT scaler = n;
	// Row<eT> diag = scaler * (1/scaler * arma::ones<Row<eT>>(n) * phi.i());
	//
	Chan<eT> phi = distance_matrix(n, n, d);
	Row<eT> diag = arma::ones<Row<eT>>(n) * phi.i();

	for(uint i = 0; i < n; i++)
		if(!less_than_or_eq(eT(0), diag(i)))
			return Chan<eT>();

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
