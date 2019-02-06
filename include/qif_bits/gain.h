
namespace g {

//		void * compare_over_prior(chan& other_channel);
//		void * compare_over_gain(chan& other_channel,Prob<eT>& prior);

template<typename eT>
inline
void check_g_size(const Mat<eT>& G, const Prob<eT>& pi) {
	if(G.n_cols != pi.n_cols)
		throw std::runtime_error("invalid prior size");
}

template<typename eT>
inline
void check_g_size(const Mat<eT>& G1, const Mat<eT>& G2) {
	if(G1.n_cols != G2.n_cols)
		throw std::runtime_error("invalid G size");
}

// max_w sum_x pi[x] G[w, x]
//
template<typename eT>
eT vulnerability(const Mat<eT>& G, const Prob<eT>& pi) {
	check_g_size(G, pi);

	return arma::max(G * trans(pi));
}

// sum_y max_w sum_x pi[x] C[x, y] G[w, x]
//
template<typename eT>
eT post_vulnerability(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	check_g_size(G, pi);
	channel::check_prior_size(pi, C);

	eT s = eT(0);
	for(uint y = 0; y < C.n_cols; y++)
		s += arma::max(G * (trans(pi) % C.col(y)));
	return s;
}

template<typename eT>
eT add_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return post_vulnerability(G, pi, C) - vulnerability(G, pi);
}

template<typename eT>
eT mult_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return post_vulnerability(G, pi, C) / vulnerability(G, pi);
}

template<typename eT>
eT mulg_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return real_ops<eT>::log2(mult_leakage(G, pi, C));
}

template<typename eT>
arma::ucolvec strategy(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	check_g_size(G, pi);
	channel::check_prior_size(pi, C);

	arma::ucolvec strategy(C.n_cols);
	for(uint y = 0; y < C.n_cols; y++)
		(G * (trans(pi) % C.col(y))).max( strategy.at(y) );

	return strategy;
}


// additive capacity for fixed pi and g ranging over 1-spanning Vg's (larger class, default) or
// 1-spanning g's (if one_spanning_g == true)
//
template<typename eT>
eT add_capacity(const Prob<eT>& pi, const Chan<eT>& C, bool one_spanning_g = false) {
	channel::check_prior_size(pi, C);

	eT res(0);

	if(one_spanning_g) {
		// 1-spanning g's (smaller class). We compute the Kantorovich distance between the
		// hypers [pi] and [pi, C], where the underlying metric between the inners is tv
		// (CSF'14, Theorem 17). Since [pi] is a point hyper, computing the transportation
		// problem is easy, since the full 1 mass of pi needs to be transferred to all
		// the posteriors. Each posterior sigma^y needs to receive mass exactly delta(y),
		// for a total cost of sum_y delta(y) * tv(pi, sigma^y)
		//
		auto tv = metric::total_variation<eT, Prob<eT>>();
		// metric::convex_separation_quasi<eT, Prob<eT>>();
		Prob<eT> outer = pi * C;

		for(uint y = 0; y < outer.n_elem; y++)
			if(!equal(outer(y), eT(0)))		// ignore 0 prob outputs
				res += outer(y) * tv(pi, channel::posterior(C, pi, y));

		
	} else {
		// For the larger class of 1-spanning Vg's, the capacity only depends on the support of pi and is
		// equal to 1 - the sum of column minima (including only rows in the support of pi).
		// The same result can also be obtained via the Kantorovich above, replacing tv with convex_separation_quasi.
		// 
		res = 1;
		for(uint y = 0; y < C.n_cols; y++) {
			eT min(1);
			for(uint x = 0; x < C.n_rows; x++)
				if(!equal(pi(x), eT(0)) && C(x,y) < min)
					min = C(x,y);
			res -= min;
		}
	}

	return res;
}


// mult leakage bound (even for negative g) coming from the miracle theorem, adjusted so that the minimum gain is exactly 0
template<typename eT>
eT mult_leakage_bound1(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	// NOTE: wrt the notation of the book  (Thm 4.7), arma::min(G, 0) is _minus_ the kappa vector we need to add to G to make it non-negative.
	// So the z here and the lambda of the book are related by lambda = z / (z-1)
	//
	eT z = arma::cdot(pi, arma::min(G, 0)) / vulnerability<eT>(G, pi); 
	// std::cout << "z:" << z << "\n";
	return bayes::mult_capacity(C) * (1-z) + z;
}

template<typename eT>
eT post_vulnerability_bound1(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return vulnerability<eT>(G, pi) * mult_leakage_bound1(G, pi, C);
}

// mult leakage bound (even for negative g) coming from the additive theorem
template<typename eT>
eT mult_leakage_bound2(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return (G.max() - G.min()) * (eT(1) - channel::sum_column_min<eT>(C)) /  vulnerability<eT>(G, pi) + eT(1);
}

// add leakage bound (even for g not bounded by 1) coming from the additive miracle theorem, adjusted so that the max gain is exactly 1
template<typename eT>
eT add_leakage_bound1(const Mat<eT>& G, const Chan<eT>& C) {
	return G.max() * (eT(1) - channel::sum_column_min<eT>(C));
}

template<typename eT>
eT post_vulnerability_bound2(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return vulnerability<eT>(G, pi) + add_leakage_bound1(G, C);
}


// returns true if the leakage of A is always >= the leakage of B (i.e. if A is refined by B)
// - if false, a counter-example G is given
//   (i.e. a gain function such that A leaks strictly less than B for the uniform prior)
// - also, the remapping channel R that minimizes the euclidean distance between AR and B is returned.
//   When the function returns true, then AR = B holds.
//
template<typename eT>
bool leakage_ge(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& G, Chan<eT>& R) {
	if(A.n_rows != B.n_rows)
		throw std::runtime_error("invalid sizes");

	uint Cr = B.n_rows;
	uint Cc = B.n_cols;
	uint Rr = A.n_cols;
	uint Rc = B.n_cols;
	uint n_vars = Cr*Cc + Rr*Rc;
	uint n_cons = Cr*Cc + Rr*Rc + Rr;

	// We want to find the point C=AR that is closest (in euclidean distance) to B.
	// If it's B itself then B refines A, otherwise G=(B-C)^T is the gain function we want.
	//
	// The (squared) euclidean distance between C and B can be written as ("." is the dot product)
	//   |C-B|_2^2 = (C-B).(C-B) = C.C -2B.C + B.B
	// B.B is constant so we need to minimize C.C -2B.C
	//
	// We have 2 types of variables:
	// - entries of C (Cr*Cc), constraint by C=AR (Cr*Cc)
	//   var C_xz has index: x*Cc+z
	// - entries of R (Rr*Rc), constraint to form a channel (Rr*Rc+Rr)
	//   var R_yz has index: Cr*Cc + y*Rc+z

	qp::QuadraticProgram<eT> qp;
	qp.l.set_size(n_cons);
	qp.u.set_size(n_cons);
	qp.c.set_size(n_vars);

	std::list<MatrixEntry<eT>> entries;		// for batch insert

	// constraints for C
	uint cons_i = 0;
	for(uint x = 0; x < Cr; x++) {
		for(uint z = 0; z < Cc; z++) {
			// We need to add the constraint: C_xz - sum_y A_xy R_yz = 0
			// We add a line to A, in which C_xz has coeff 1, each R_yz has coeff - Axy, and the bounds are 0

			qp.l(cons_i) = 0;	// lower bound
			qp.u(cons_i) = 0;	// upper bound

			entries.push_back(MatrixEntry<eT>(cons_i, x*Cc+z, 1));	// coeff of C_xz

			for(uint y = 0; y < Rr; y++)
				entries.push_back(MatrixEntry<eT>(cons_i, Cr*Cc + y*Rc+z, -A(x,y)));	// coeff of R_yz

			cons_i++;
		}
	}

	// constraints for R
	for(uint y = 0; y < Rr; y++) {
		for(uint z = 0; z < Rc; z++) {
			// We need to add the constraint: 0 <= R_yz <= 1
			// We add a line to A, in which R_yz has coeff 1, and the bounds are 0,1
			qp.l(cons_i) = 0;	// lower bound
			qp.u(cons_i) = 1;	// upper bound
			entries.push_back(MatrixEntry<eT>(cons_i, Cr*Cc + y*Rc+z, 1));

			cons_i++;
		}

		// We need to add the constraint: sum_z R_yz = 1
		// We add a line to A, in which each R_yz has coeff 1, and the bounds are both 1
		qp.l(cons_i) = 1;	// lower bound
		qp.u(cons_i) = 1;	// upper bound

		for(uint z = 0; z < Rc; z++)
			entries.push_back(MatrixEntry<eT>(cons_i, Cr*Cc + y*Rc+z, 1));

		cons_i++;
	}

	assert(cons_i == n_cons);
	fill_spmat(qp.A, n_cons, n_vars, entries);

	// cost, quadratic part: C.C
	// coeff 2 for the C variables (because of the 1/2 in the QP), 0 for the R variables
	entries.clear();
	for(uint i = 0; i < Cr*Cc; i++)
		entries.push_back(MatrixEntry<eT>(i, i, 2));
	fill_spmat(qp.P, n_vars, n_vars, entries);

	// cost, linear part: -2BC
	// coeffs -2B for the C variables, 0 for the R variables
	qp.c.fill(0);
	for(uint x = 0; x < Cr; x++)
		for(uint z = 0; z < Cc; z++)
			qp.c(x*Cc+z) = -2 * B(x,z);

	// precision of the solution (Kostas: not 100% sure here)
	// if no value is set by the user, use 1e-5 instead of OSQP's defaults
	if(qp::Defaults::osqp_eps_abs < 0.0)	qp.osqp_eps_abs = 1e-5;
	if(qp::Defaults::osqp_eps_rel < 0.0)	qp.osqp_eps_rel = 1e-5;

	// ready
	if(!qp.solve())
		throw std::runtime_error("leakage_ge: QP infeasible, this shouldn't happen");

	// add B.B to the cost function to obtained the squared distance (see the program definition above)
	eT dist = qp.optimum() + arma::dot(B, B);

	bool res = equal(dist, eT(0), eT(qp.osqp_eps_abs), eT(qp.osqp_eps_rel));
	if(res) {
		G.clear();
	} else {
		G = B.t() - arma::reshape(qp.x.head(Cc*Cr), Cc, Cr);
		G -= G.min();	// non-negative
		G /= G.max();	// and in [0,1]
	}

	R = arma::reshape(qp.x.tail(Rc*Rr), Rc, Rr).t();

	return res;
}

template<typename eT>
bool leakage_ge(const Chan<eT>& A, const Chan<eT>& B) {
	Mat<eT> G;
	Chan<eT> R;
	return leakage_ge<eT>(A, B, G, R);
}

template<typename eT>
bool leakage_ge(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& G) {
	Chan<eT> R;
	return leakage_ge<eT>(A, B, G, R);
}


////////////////// Gain function manipulation //////////////////////////


// Adding g1+g2 produces a g such that Vg = Vg1 + Vg2
//
template<typename eT>
Mat<eT> g_add(const Mat<eT>& G1, const Mat<eT>& G2) {
	check_g_size(G1, G2);

	Mat<eT> G(G1.n_rows * G2.n_rows, G1.n_cols);
	G.fill(0);

	for(uint i = 0; i < G1.n_rows; i++)
		for(uint j = 0; j < G2.n_rows; j++)
			G.row(i * G2.n_rows + j) = G1.row(i) + G2.row(j);

	return G;
}

// g such that Vg(pi) = Vg'(pi,C)
//
template<typename eT>
Mat<eT> g_from_posterior(const Mat<eT>& G, const Chan<eT>& C) {
	if(G.n_cols != C.n_rows)
		throw std::runtime_error("invalid G size");

	Mat<eT> Gres(1, G.n_cols);
	Gres.fill(0);

	for(uint y = 0; y < C.n_cols; y++) {
		Mat<eT> Gtemp = G;
		Gtemp.each_row() %= arma::trans(C.col(y));
		Gres = add(Gres, Gtemp);
	}

	return Gres;
}

} // namespace g
