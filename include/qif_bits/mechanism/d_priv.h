
namespace mechanism::d_priv {

// build a matrix D where Dij = exp(-d(i,j))
// used in the construction of most mechanisms below
//
template<typename eT>
Mat<eT>
distance_matrix(uint n_rows, uint n_cols, Metric<eT, uint> d) {
	Mat<eT> D(n_rows, n_cols);

	// don't assume symmetry or anything about d (no big gain, and there are cases when assymetric distances are useful)
	//
	for(uint i = 0; i < n_rows; i++) {
		for(uint j = 0; j < n_cols; j++) {
			eT expon = -d(i, j);			// separate variable needed for rats
			D(i, j) = qif::exp(expon);
		}
	}
	return D;
}

// first_y,first_x are the number corresponding to the first column/row
//
template<typename eT>
Chan<eT>
geometric(uint n_rows, eT epsilon = 1.0, uint n_cols = 0, int first_x = 0, int first_y = 0) {
	if(n_cols == 0) n_cols = n_rows;
	if(n_cols < 2)      throw std::runtime_error("n_cols should be at least 2");

	// shift x's and y's so that first_x is 0. This does not affect the matrix but makes formulas easier
	first_y -= first_x;
	first_x = 0;

	// construct a metric between inputs and _outputs_
	Metric<eT,uint> d = [&](uint x, uint y) -> eT {
		return epsilon * abs((signed)(y + first_y - x));
	};

	eT c = qif::exp(epsilon),
	   lambda_f = c / (c + eT(1)),				// lambda for the first/last cols
	   lambda_m = (c - eT(1)) / (c + eT(1));	// lambda for the middle colums

	// the standard formula requires the "peaks" for each x to be included in the colums.
	// If the are not, we add extra columns and manually truncate them afterwards.
	uint trunc_left = std::max(first_y, 0),										// we left-truncate when first_y > 0
		 trunc_right = std::max((signed)(n_rows - first_y - n_cols), 0);		// we right-truncate when first_y+n_cols < n_rows
	n_cols += trunc_left + trunc_right;
	first_y -= trunc_left;

	Chan<eT> C = distance_matrix(n_rows, n_cols, d);
	C.col(0)        	    *= lambda_f;
	C.col(n_cols-1)			*= lambda_f;
	if(n_cols > 2)
		C.cols(1, n_cols-2)	*= lambda_m;

	// if we need to truncate, move the probabilities to the first/last col and remove the truncated columns
	if(trunc_left) {
		C.col(trunc_left) += arma::sum(C.cols(0, trunc_left-1), 1);
		C.shed_cols(0, trunc_left-1);
		n_cols -= trunc_left;
		first_y += trunc_left;
	}
	if(trunc_right) {
		C.col(n_cols-trunc_right-1) += arma::sum(C.cols(n_cols-trunc_right, n_cols-1), 1);
		C.shed_cols(n_cols-trunc_right, n_cols-1);
		n_cols -= trunc_right;
	}

	return C;
}

// d should be a metric between inputs and outputs
//
template<typename eT>
Chan<eT>
exponential(uint n_rows, Metric<eT, uint> d, uint n_cols = 0) {
	if(n_cols == 0) n_cols = n_rows;

	Chan<eT> C = distance_matrix(n_rows, n_cols, eT(1)/2 * d);
	channel::normalize(C);

	return C;
}

template<typename eT>
Chan<eT>
randomized_response(uint n_rows, eT epsilon = 1.0, uint n_cols = 0) {
	if(n_cols == 0) n_cols = n_rows;

	// Essentially the exponential for the discrete metric.
	// But due to the symmetry we don't need to half the epsilon!
	//
	auto d = epsilon * metric::discrete<eT,uint>();
	Chan<eT> C = distance_matrix(n_rows, n_cols, d);
	channel::normalize(C);

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

// A mechanism C such that for all x,x', mult_total_var(C.row(x), C.row(x')) = d(x,x') (exactly)
//
template<typename eT>
Chan<eT>
exact_distance(uint n, Metric<eT, uint> d) {
	Chan<eT> C(n, n+1);

	for(uint i = 0; i < n; i++) {
	for(uint j = 0; j < n; j++) {
		C(i,j) = eT(1)/n * qif::exp( -d(j,i)-1 );
	}}

	C.col(n).fill(0);
	C.col(n) = 1 - arma::sum(C, 1);

	return C;
}

} // namespace mechanism::d_priv
