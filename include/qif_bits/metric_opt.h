// metric optimization problems

namespace metric::opt {

// Computes the l1-diameter of the set of C's rows.
// Returns also the rows that produce the diameter.
//
template<typename eT = eT_def>
std::tuple<eT,uint,uint> l1_diameter(const Chan<eT>& C, std::string method = "direct") {
	eT diam(0);
	uint res_x1, res_x2;

	if(method == "direct") {
		auto l1 = qif::metric::l1<eT, Prob<eT>>();

		for(uint x1 = 0; x1 < C.n_rows; x1++) {
			for(uint x2 = x1+1; x2 < C.n_rows; x2++) {
				eT l1_cur = l1(C.row(x1), C.row(x2));

				if(less_than(diam, l1_cur)) {
					diam = l1_cur;
					res_x1 = x1;
					res_x2 = x2;
				}
			}
		}
	
	} else if(method == "linf") {
		// Computes the diameter by embeding C's rows in Rinf^d where d = 2^m (m is n_cols).
		// The embedding phi(x) in R^d of a vector x in R^m has one coordinate for every bitstring b in {0,1}^m.
		// The value of that coordinate is cdot(x, coeff) where coeff_i = (-1)^{b_i}.
		//
		// In Linf, the diamater is simply given by max_i (max_x x_i - min_x x_i). So for each coordinate b, we need to
		// compute { phi(x)_b }_x, then compute the max - min, and finally keep the maximum of those.
		//
		eT one(1);

		arma::Col<eT> coeff(C.n_cols);	// we only keep the coeff vector (b is implicit)
		coeff.fill(one);				// start with b = 00..0, coeff = (1,...,1)

		// we have a loop for every value of b. The loop finishes when b becomes 11..1, in which case
		// the update loop (which stops at the first 0) will finish all steps and i will become C.n_cols
		uint i = 0;
		while(i < C.n_cols) {
			// with a single multiplication we get cdot(coeff, row) for all rows.
			// we need to compute max-min, and update diam if the new value is bigger
			auto trans = C * coeff;
			if(eT d = arma::range(trans); less_than(diam, d)) {
				diam = d;
				res_x1 = arma::index_min(trans);
				res_x2 = arma::index_max(trans);
			}

			// compute the coeff for the next value of b. For b we would need to change 0s to 1s and vice versa, until we find the first 0 (which became 1).
			// In coeff, this translates changing 1s to -1s and vice versa, until we find the first -1 (which became 1).
			for(i = 0; i < C.n_cols && (coeff(i) *= -one) > eT(0); i++)
				;
		}

	} else {
		throw std::runtime_error("invalid method: " + method);
	}

	return std::tuple(diam, res_x1, res_x2);
}

// Computes the (unique) vector q that minimizes the max l2-distance from the rows of C.
// The vector is returned in q, the radius is the return value of the function.
//
template<typename eT = eT_def>
inline
eT min_l2_enclosing_ball(const Chan<eT>& C, Prob<eT>& q) {
	uint M = C.n_rows,
		 N = C.n_cols;

	// miniball doc says to always use double
	typedef Seb::Point<double> Point;
	typedef Seb::Smallest_enclosing_ball<double> Miniball;

	// prepare vectors
	std::vector<Point> S;
	for(uint x = 0; x < M; x++)
		S.push_back(Point(  N, arma::conv_to<std::vector<double>>::from(C.row(x)).begin() ));

	// solve
	Miniball mb(N, S);

	auto center_it = mb.center_begin();
	q.set_size(N);
	for(uint y = 0; y < N; y++)
		q(y) = center_it[y];

	return mb.radius();
}

} // namespace metric::opt