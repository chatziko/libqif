// metric optimization problems

namespace metric::optimize {

// Computes the l1-diameter of the set of C's rows.
// Returns also the rows that produce the diameter.
//
template<typename eT = eT_def>
std::tuple<eT,uint,uint> l1_diameter(const Chan<eT>& C, std::string method = "direct") {
	eT diam(0);
	uint res_x1 = 0, res_x2 = 0;

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

	return { diam, res_x1, res_x2 };
}

// Computes the max l1-distance between the rows of A and B.
// Returns also the rows that produce the max distance.
//
template<typename eT = eT_def>
std::tuple<eT,uint,uint> l1_max_distance(const Chan<eT>& A, const Chan<eT>& B, std::string method = "direct") {
	eT maxdist(0);
	uint res_x1, res_x2;

	if(method == "direct") {
		auto l1 = qif::metric::l1<eT, Prob<eT>>();

		for(uint x1 = 0; x1 < A.n_rows; x1++) {
			for(uint x2 = 0; x2 < B.n_rows; x2++) {
				eT l1_cur = l1(A.row(x1), B.row(x2));

				if(less_than(maxdist, l1_cur)) {
					maxdist = l1_cur;
					res_x1 = x1;
					res_x2 = x2;
				}
			}
		}

	} else {
		// TODO: implement linf method like for the diameter
		throw std::runtime_error("invalid method: " + method);
	}

	return { maxdist, res_x1, res_x2 };
}

// Computes the (unique) vector q that minimizes the max l2-distance from the rows of C.
// The radius and the vector q are returned.
//
template<typename eT = eT_def>
inline
std::pair<eT,Prob<eT>> l2_min_enclosing_ball(const Chan<eT>& C) {
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
	Prob<eT> q(N);
	for(uint y = 0; y < N; y++)
		q(y) = center_it[y];

	return { mb.radius(), q };
}

// Computes the probability distribution q that minimizes the max l1-distance from those in C.
// Note: this is different than the l1-SEB in the whole R^n. The result is not guaranteed to be in the
//       convex hull of C's rows, so we have to explicitly constraint q to be in the probability simplex.
//
// The radius and the vector q are returned.
// in_conv_hull forces solution to be within the convex hull of C's rows (only for the "lp" method)
//
template<typename eT = eT_def>
inline
std::pair<eT,Prob<eT>> simplex_l1_min_enclosing_ball(const Chan<eT>& C, std::string method = "lp", bool in_conv_hull = false) {
	uint M = C.n_rows,
		 N = C.n_cols;
	eT inf = infinity<eT>();

	// q: N variables
	lp::LinearProgram<eT> lp;
	auto vars = lp.make_vars(N, eT(0), eT(1));

	// equalities for summing up to 1
	//
	auto con = lp.make_con(eT(1), eT(1));
	for(uint y = 0; y < N; y++)
		lp.set_con_coeff(con, vars[y], eT(1));


	if(method == "lp") {
		// vars diff_vars[x,y] >= |C[x,y] - q[y]| for each x,y
		//
		auto diff_vars = lp.make_vars(M, N, eT(0), eT(2));

		for(uint x = 0; x < M; x++) {
			for(uint y = 0; y < N; y++) {
				// diff_vars[x,y] >=   C[x,y] - q[y]
				auto con = lp.make_con(C(x,y), inf);
				lp.set_con_coeff(con, diff_vars[x][y], eT(1));
				lp.set_con_coeff(con, vars[y], eT(1));

				// diff_vars[x,y] >= - C[x,y] + q[y]
				con = lp.make_con(-C(x,y), inf);
				lp.set_con_coeff(con, diff_vars[x][y], eT(1));
				lp.set_con_coeff(con, vars[y], eT(-1));
			}
		}

		// We set z >= sum_y diff_vars[x,y] >= sum_y |C[x,y] - q[y]| for each x
		//
		auto z = lp.make_var(eT(0), eT(2));

		for(uint x = 0; x < M; x++) {
			auto con = lp.make_con(eT(0), inf);
			lp.set_con_coeff(con, z, eT(1));

			for(uint y = 0; y < N; y++)
				lp.set_con_coeff(con, diff_vars[x][y], eT(-1));
		}

		// cost function: minimize z
		//
		lp.maximize = false;
		lp.set_obj_coeff(z, eT(1));

		// force solution to be within the convex hull of C
		//
		if(in_conv_hull) {
			// convex coefficients, sum up to 1
			auto coeff = lp.make_vars(M, eT(0), eT(1));

			auto con = lp.make_con(eT(1), eT(1));
			for(uint x = 0; x < M; x++)
				lp.set_con_coeff(con, coeff[x], eT(1));

			// vars[y] = sum_x coeff[x] C[x,y]  for each y
			for(uint y = 0; y < N; y++) {
				auto con = lp.make_con(eT(0), eT(0));
				lp.set_con_coeff(con, vars[y], eT(1));

				for(uint x = 0; x < M; x++)
					lp.set_con_coeff(con, coeff[x], -C(x,y));
			}
		}

	} else if (method == "linf") {
		// Computes the l1 enclosing ball by embeding C's rows in Rinf^d where d = 2^m (m is n_cols).
		// The embedding phi(x) in R^d of a vector x in R^m has one coordinate for every bitstring b in {0,1}^m.
		// The value of that coordinate is cdot(x, coeff) where coeff_i = (-1)^{b_i}.
		//
		if(in_conv_hull)
			throw std::runtime_error("in_conv_hull not supported for the linf method");

		// Bor each bitstring b we add constraints
		//   z >=   max_x{ phi(Cx)_b } - phi(q)_b
		//   z >= - min_x{ phi(Cx)_b } + phi(q)_b
		//
		auto z = lp.make_var(eT(-2), eT(2));

		arma::Col<eT> coeff(N);			// we only keep the coeff vector (b is implicit)
		coeff.fill(eT(0));				// start with b = 00..0, coeff = (1,...,1)

		for(uint i = 0; i < N; ) {
			// with a single multiplication we get cdot(coeff, row) for all rows.
			// we compute max-min for the constraints
			auto trans = C * coeff;

			auto con1 = lp.make_con( arma::max(trans), inf);
			auto con2 = lp.make_con(-arma::min(trans), inf);

			lp.set_con_coeff(con1, z, eT(1));
			lp.set_con_coeff(con2, z, eT(1));

			for(uint y = 0; y < N; y++) {
				lp.set_con_coeff(con1, vars[y],  coeff(y));
				lp.set_con_coeff(con2, vars[y], -coeff(y));
			}

			// compute the coeff for the next value of b
			for(i = 0; i < C.n_cols && (coeff(i) *= eT(-1)) > eT(0); i++)
				;
		}

		// cost function: minimize z
		//
		lp.maximize = false;
		lp.set_obj_coeff(z, eT(1));

	} else {
		throw std::runtime_error("invalid method: " + method);
	}

	// solve program
	//
	if(!lp.solve())
		throw std::runtime_error("simplex_l1_min_enclosing_ball: lp should be always solvable");

	// reconstruct q from solution
	//
	Prob<eT> q(N);
	for(uint y = 0; y < N; y++)
		q(y) = lp.solution(vars[y]);

	return { lp.objective(), q };
}

// Computes the (euclidean) projection of vector x \in R^n onto the
// n-simplex. That is, finds the point of the simplex that is closer  to x.
// Uses the algorithm of:
//    http://www.springerlink.com/content/q1636371674m36p1/
//
template<typename eT = eT_def>
inline
Prob<eT> simplex_project(Prob<eT> x) {
	uint n = x.n_cols;
	Row<char> done(n);
	done.fill(0);

	while(true) {
		eT t = (arma::accu(x)-1)/n;
		uint n1 = 0;

		for(uint i = 0; i < x.n_cols; i++) {
			if(done(i)) continue;

			x(i) -= t;
			if(eT(0) > x(i)) {
				x(i) = eT(0);
				done(i) = 1;
				n1++;
			}
		}

		if(n1 == 0) break;		// no negative elements
		n -= n1;
	}

	return x;
}

} // namespace metric::optimize