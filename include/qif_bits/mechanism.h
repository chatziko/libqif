
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

// d should be a scaled version of metric::euclidean<eT,uint>()
template<typename eT>
Chan<eT>
geometric(uint n_rows, Metric<eT, uint> d = metric::euclidean<eT,uint>(), uint n_cols = 0) {
	if(n_cols == 0) n_cols = n_rows;
	if(n_rows > n_cols) throw std::runtime_error("n_cols should be at least as big as n_rows");
	if(n_rows < 2)      throw std::runtime_error("n_rows should be at least 2");

	eT c = std::exp(d(0,1)),
	   lambda_f = c / (c + eT(1)),				// lambda for the first/last cols
	   lambda_m = (c - eT(1)) / (c + eT(1));	// lambda for the middle colums

	Chan<eT> C = distance_matrix(n_rows, n_cols, d);
	C.col(0)        	    *= lambda_f;
	C.col(n_cols-1)			*= lambda_f;
	if(n_cols > 2)
		C.cols(1, n_cols-2)	*= lambda_m;

	return C;
}

template<typename eT>
Chan<eT>
exponential(uint n_rows, Metric<eT, uint> d, uint n_cols = 0) {
	if(n_cols == 0) n_cols = n_rows;

	Chan<eT> C = distance_matrix(n_rows, n_cols, eT(1)/2 * d);
	C.each_col() /= sum(C, 1);	// normalize

	return C;
}

template<typename eT>
Chan<eT>
tight_constraints(uint n, Metric<eT, uint> d) {
	// Note: in the perl code, scaling the ones(n) vector somehow improves numerical stability
	//       (the non-negativity of diag). We should investigate if this still happens
	// eT scaler = n;
	// Row<eT> diag = scaler * (1/scaler * arma::ones<Row<eT>>(n) * phi.i());
	//
	// alternative way by inverting phi
	// Col<eT> diag = phi.i() * arma::ones<Col<eT>>(n);
	//
	// create diagonal by solving phi X = 1
	Chan<eT> phi = distance_matrix(n, n, d);
	Col<eT> diag;
	arma::solve(diag, phi, arma::ones<Col<eT>>(n));

	// the mechanism exists if the system has a non-negative solution
	if(diag.is_empty())
		return Chan<eT>();
	for(auto e : diag)
		if(!less_than_or_eq(eT(0), e))
			return Chan<eT>();

	// the channel matrix = phi after multiplying all rows with diag
	phi.each_row() %= diag.t();

	return phi;
}

// tight constraints, with the extra constraint that the coefficients of some colums are required
// to be the same. These constraints are given in form of a vector cols. If cols[i] == j it means
// that the coefficient of column i should be equal to column j. Set cols[i] = i to have an
// unconstraint column.
// Eg: (0, 0, 2, 2) means that there are 2 free coefficients (cols 0, 2) and two constrained (1, 3)
//
template<typename eT>
Chan<eT>
tight_constraints(arma::uvec cols, Metric<eT, uint> d) {
	uint n = cols.n_elem;
	Chan<eT> phi = distance_matrix(n, n, d);

	arma::uvec unique = arma::unique(cols);
	uint nu = unique.n_elem;

	// create a remapping channel based on the remap of cols
	Chan<eT> remap(n, nu, arma::fill::zeros);
	for(uint i = 0; i < n; i++) {
		uint target = arma::find(unique == cols(i)).eval()(0);
		remap(i, target) = eT(1);
	}

	// find the diag for the remapped phi
	Mat<eT> r_phi = phi * remap;
	Col<eT> r_diag;
	arma::solve(r_diag, r_phi, arma::ones<Col<eT>>(r_phi.n_rows));

	// the mechanism exists if the system has a non-negative solution
	if(r_diag.is_empty())
		return Chan<eT>();
	for(auto e : r_diag)
		if(!less_than_or_eq(eT(0), e))
			return Chan<eT>();

	// put the elements that were actually computed in the diag
	Col<eT> diag(n);
	for(uint i = 0; i < n; i++) {
		uint target = arma::find(unique == cols(i)).eval()(0);
		diag(i) = r_diag(target);
	}

	// the channel matrix = phi after multiplying all rows with diag
	phi.each_row() %= diag.t();

	return phi;
}

template<typename eT>
bool is_private(const Chan<eT>& C, Metric<eT, uint> d) {
	auto mtv = metric::mult_total_variation<eT, Prob<eT>>();

	for(uint i = 0; i < C.n_rows; i++) {
		for(uint j = i+1; j < C.n_rows; j++) {
			// chainable elements are redundant to check
			if(d.chainable(i, j)) continue;

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
			// chainable elements are redundant to check
			if(d.chainable(i, j)) continue;

			eT ratio = mtv(C.row(i), C.row(j)) / d(i, j);
			if(less_than(res, ratio))
				res = ratio;
		}
	}
	return res;
}

template<typename eT>
eT d_vulnerability(Metric<eT, uint> d, const Prob<eT>& pi) {

	eT res(0);
	for(uint i = 0; i < pi.n_cols; i++) {
		for(uint j = i+1; j < pi.n_cols; j++) {
			// chainable elements are redundant to check
			if(d.chainable(i, j)) continue;

			eT ratio = std::abs(std::log(pi(i)) - std::log(pi(j))) / d(i, j);
			if(less_than(res, ratio))
				res = ratio;
		}
	}
	return res;
}

} // namespace mechanism
