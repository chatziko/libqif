
namespace probab {

template<typename eT = eT_def>
inline
Prob<eT>& uniform(Prob<eT>& pi) {
	// cast to uint cause rat is confused when dividing by const uint
	pi.fill( eT(1) / uint(pi.n_cols) );
	return pi;
}

template<typename eT = eT_def>
inline Prob<eT> uniform(uint n) {
	Prob<eT> pi(n);
	uniform(pi);
	return pi;		// separate return to allow move semantics!
}


template<typename eT = eT_def>
inline
Prob<eT>& dirac(Prob<eT>& pi, uint i = 0) {
	pi.zeros();
	pi.at(i) = eT(1);
	return pi;
}

template<typename eT = eT_def>
inline
Prob<eT> dirac(uint n, uint i = 0) {
	Prob<eT> pi(n);
	dirac(pi, i);
	return pi;		// separate return to allow move semantics!
}


// Generate a distribution of n elements uniformly (among all the elements of the n-1 simplex)
//
// Algorithm: Generate n-1 numbers uniformly in [0,1], add 0 and 1 to the list,
// sort, and take the difference between consecutive elements.
// This algorithm has advantages over the simple "normalize n uniform elements" tehcnique:
//
// 1. The normalizing technique is not uniform! See the url below.
//
// 2. The algorithm involves no divisions, and the resulting sum is much closer to exactly 1.0
//    Using this algorithm with float the Kantorovich tests over random dists always pass, while with
//    the normalizing algorithm there are instabilities.
//
// 3. We avoid to sum all elements, which creates a big denominator under rats
//
// Note: if the n logn complexity is a problem, there's a linear algorithm involving logs in the following url:
// http://stats.stackexchange.com/questions/14059/generate-uniformly-distributed-weights-that-sum-to-unity
//
template<typename eT = eT_def>
inline
Prob<eT>& randu(Prob<eT>& pi) {
	pi.randu();
	pi(pi.n_cols-1) = eT(1);		// add 1 to the list. We don't really need to add 0
	pi = arma::sort(pi);

	for(uint i = pi.n_cols-1; i > 0; i--)
		pi(i) -= pi(i-1);

	// naif normalize algorithm
	//pi.randu();
	//pi /= arma::accu(pi);

	return pi;
}

template<typename eT = eT_def>
inline
Prob<eT> randu(uint n) {
	Prob<eT> pi(n);
	randu(pi);
	return pi;		// separate return to allow move semantics!
}

template<typename eT = eT_def>
inline void normalize(Prob<eT>& pi) {
	pi /= arma::accu(pi);
}


// draw from pi. Allow pi to be anything iterable
//
template<typename eT = eT_def, typename T = Prob<eT>>
inline
uint draw(const T& pi) {
	eT p = rng::randu<eT>();

	eT accu(0);
	uint cur = 0;
	for(auto pi_it = pi.begin(); pi_it != pi.end(); ++pi_it) {
		if(!less_than(accu += *pi_it, p))
			return cur;
		cur++;
	}

	return pi.n_cols - 1;
}

// draw multiple samples efficiently (with a single iteration over pi)
//
template<typename eT = eT_def, typename T = Prob<eT>>
inline
Row<uint> draw(const T& pi, uint n) {
	// we need n numbers uniformly sampled in [0,1]. Sort them and keep the indexes of the sorted list in 'order'
	Row<eT> ps(n);
	ps.randu();
	arma::uvec order = arma::sort_index(ps);

	eT accu(0);
	int cur = -1;
	auto pi_it = pi.begin();	// single iteration over pi for all samples
	Row<uint> res(n);

	// sample n elements. 
	for(uint i = 0; i < n; i++) {
		// we need to visit elements in sorted order of p
		eT p = ps(order(i));

		while(less_than(accu, p) && pi_it != pi.end()) {
			accu += *pi_it;
			pi_it++;
			cur++;
		}

		// place the result in the same position in res as p was in ps.
		res(order(i)) = cur;
	}

	return res;
}

template<typename eT = eT_def>
inline
bool is_uniform(const Prob<eT>& pi, const eT& mrd = def_mrd<eT>) {
	eT v = eT(1) / (int)pi.n_cols;
	for(uint j = 0; j < pi.n_cols; j++)
		if(!equal(v, pi(j), eT(0), mrd))
			return false;

	return true;
}

template<typename eT = eT_def>
inline
bool is_proper(const Prob<eT>& pi, const eT& mrd = def_mrd<eT>) {
	eT sum(0);
	for(uint j = 0; j < pi.n_cols; j++) {
		// elements should be non-negative
		const eT& elem = pi.at(j);
		if(less_than(elem, eT(0)))
			return false;

		sum += elem;
	}

	// sum should be 1
	if(!equal(sum, eT(1), def_md<eT>, mrd))
		return false;

	return true;
}

template<typename eT = eT_def>
inline
void assert_proper(const Prob<eT>& pi) {
	if(!is_proper(pi))
		throw std::runtime_error("not a proper dist");
}

template<typename eT = eT_def>
inline bool equal(const Prob<eT>& A, const Prob<eT>& B, const eT& md = def_md<eT>, const eT& mrd = def_mrd<eT>) {
	if(A.n_cols != B.n_cols)
		return false;

	for(uint i = 0; i < A.n_cols; i++)
		if(!qif::equal(A.at(i), B.at(i), md, mrd))
			return false;

	return true;
}



// Computes in place the (euclidean) projection of vector x \in R^n onto the
// n-simplex. That is, finds the point of the simplex that is closer  to x.
// Uses the algorithm of:
//    http://www.springerlink.com/content/q1636371674m36p1/
//
template<typename eT = eT_def>
inline
Prob<eT>& project_to_simplex(Prob<eT>& x) {
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

// Computes the probability distribution q that minimizes the max l1-distance from those in C.
// The distribution is returned in q, the actual distances is the return value of the function.
// in_conv_hull forces solution to be within the convex hull of C's rows
//
template<typename eT = eT_def>
inline
eT min_l1_enclosing_ball(const Chan<eT>& C, Prob<eT>& q, bool in_conv_hull = false) {
	uint M = C.n_rows,
		 N = C.n_cols;
	eT inf = infinity<eT>();

	// q: N variables
	lp::LinearProgram<eT> lp;
	auto vars = lp.make_vars(N, eT(0), eT(1));

	// vars diff_vars[x,y] >= |C[x,y] - q[y]| for each x,y
	//
	auto diff_vars = lp.make_vars(M, N, 0, 2);

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
	auto z = lp.make_var(0, 2);

	for(uint x = 0; x < M; x++) {
		auto con = lp.make_con(0, inf);
		lp.set_con_coeff(con, z, eT(1));

		for(uint y = 0; y < N; y++)
			lp.set_con_coeff(con, diff_vars[x][y], eT(-1));
	}

	// cost function: minimize z
	//
	lp.maximize = false;
	lp.set_obj_coeff(z, eT(1));

	// equalities for summing up to 1
	//
	auto con = lp.make_con(1, 1);
	for(uint y = 0; y < N; y++)
		lp.set_con_coeff(con, vars[y], 1);

	// force solution to be within the convex hull of C
	//
	if(in_conv_hull) {
		// convex coefficients, sum up to 1
		auto coeff = lp.make_vars(M, eT(0), eT(1));

		auto con = lp.make_con(1, 1);
		for(uint x = 0; x < M; x++)
			lp.set_con_coeff(con, coeff[x], 1);

		// vars[y] = sum_x coeff[x] C[x,y]  for each y
		for(uint y = 0; y < N; y++) {
			auto con = lp.make_con(0, 0);
			lp.set_con_coeff(con, vars[y], 1);

			for(uint x = 0; x < M; x++)
				lp.set_con_coeff(con, coeff[x], -C(x,y));
		}
	}

	// solve program
	//
	if(!lp.solve())
		throw std::runtime_error("min_l1_enclosing_ball: lp should be always solvable");

	// reconstruct q from solution
	//
	q.set_size(N);
	for(uint y = 0; y < N; y++)
		q(y) = lp.solution(vars[y]);

	return lp.objective();
}

// Computes the l1 enclosing ball by embeding C's rows in Rinf^d where d = 2^m (m is n_cols).
// The embedding phi(x) in R^d of a vector x in R^m has one coordinate for every bitstring b in {0,1}^m.
// The value of that coordinate is cdot(x, coeff) where coeff_i = (-1)^{b_i}.
//
template<typename eT = eT_def>
inline
eT min_l1_enclosing_ball_via_linf(const Chan<eT>& C, Prob<eT>& q) {
	uint N = C.n_cols;
	eT one(1);
	eT inf = infinity<eT>();

	// linear program
	// q: N variables
	lp::LinearProgram<eT> lp;
	auto vars = lp.make_vars(N, eT(0), eT(1));

	// Bor each bitstring b we add constraints
	//   z >=   max_x{ phi(Cx)_b } - phi(q)_b
	//   z >= - min_x{ phi(Cx)_b } + phi(q)_b
	//
	auto z = lp.make_var(-2, 2);

	arma::Col<eT> coeff(N);			// we only keep the coeff vector (b is implicit)
	coeff.fill(one);				// start with b = 00..0, coeff = (1,...,1)

	for(uint i = 0; i < N; ) {
		// with a single multiplication we get cdot(coeff, row) for all rows.
		// we compute max-min for the constraints
		auto trans = C * coeff;

		auto con1 = lp.make_con( arma::max(trans), inf);
		auto con2 = lp.make_con(-arma::min(trans), inf);

		lp.set_con_coeff(con1, z, 1);
		lp.set_con_coeff(con2, z, 1);

		for(uint y = 0; y < N; y++) {
			lp.set_con_coeff(con1, vars[y],  coeff(y));
			lp.set_con_coeff(con2, vars[y], -coeff(y));
		}

		// compute the coeff for the next value of b
		for(i = 0; i < C.n_cols && (coeff(i) *= -one) > eT(0); i++)
			;
	}

	// cost function: minimize z
	//
	lp.maximize = false;
	lp.set_obj_coeff(z, eT(1));

	// equalities for summing up to 1
	//
	auto con = lp.make_con(1, 1);
	for(uint y = 0; y < N; y++)
		lp.set_con_coeff(con, vars[y], 1);

	// solve program
	//
	if(!lp.solve())
		throw std::runtime_error("min_l1_enclosing_ball_via_linf: lp should be always solvable");

	// reconstruct q from solution
	//
	q.set_size(N);
	for(uint y = 0; y < N; y++)
		q(y) = lp.solution(vars[y]);

	return lp.objective();
}

// Takes a probability pi on a grid of given width.
// Returns a grid representation of that probability.
//
template<typename eT = eT_def>
Mat<eT> to_grid(const Prob<eT>& pi, uint width) {
	return arma::trans(arma::reshape(pi, pi.n_elem / width, width));
}

// Vectorizes a probability on a grid.
// Element (i,j) is mapped to i+j*width.
//
template<typename eT = eT_def>
Prob<eT> from_grid(const Mat<eT>& grid) {
	return arma::vectorise(grid, 1);
}

} // namespace probab
