
namespace refinement {

// forward declaration
template<typename eT> eT add_metric(const Prob<eT>&, const Chan<eT>&, const Chan<eT>&, Mat<eT>&);


// true if A is refined by B (i.e. A's leakage is >= B's)
//
template<typename eT = eT_def>
inline
bool refined_by(const Chan<eT>& A, const Chan<eT>& B) {
	// true if AX = B for some X
	auto X = channel::factorize(B, A);
	return !X.empty();
}


// same as refined_by(A, B), using the method of projecting B to { AR | R } (via quadratic programming)
// - if false, a counter-example G is given
//   (i.e. a gain function such that A leaks strictly less than B for the uniform prior)
// - also, the remapping channel R that minimizes the euclidean distance between AR and B is returned.
//   When the function returns true, then AR = B holds.
//
// Note1: It should be possible to do the same thing by projecting wrt the manhattan norm (see Boyd 8.1.2).
//        This can be done via LP so it will likely be faster
// Note2: refined_by is equivalent to add_metric(uni, A, B) == 0, so this also gives an LP solution, but
//        it turns out to be much slower than the projection method (although the latter is done via QP).
//
template<typename eT>
bool refined_by(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& G, Chan<eT>& R) {
	if(A.n_rows != B.n_rows)
		throw std::runtime_error("invalid sizes");

	uint Cr = B.n_rows;
	uint Cc = B.n_cols;
	uint Rr = A.n_cols;
	uint Rc = B.n_cols;

	// We want to find the point C=AR that is closest (in euclidean distance) to B.
	// If it's B itself then B refines A, otherwise G=(B-C)^T is the gain function we want.
	//
	// The (squared) euclidean distance between C and B can be written as ("." is the dot product)
	//   |C-B|_2^2 = (C-B).(C-B) = C.C -2B.C + B.B
	// B.B is constant so we need to minimize C.C -2B.C
	//
	// We have 2 types of variables:
	// - entries of C (Cr x Cc), constrained by C=AR
	// - entries of R (Rr x Rc), constrained to form a channel

	qp::QuadraticProgram<eT> qp;

	auto C = qp.make_vars(Cr, Cc);
	auto RR = qp.make_vars(Rr, Rc);			// RR to distringuish from param R

	// constraints for C
	for(uint x = 0; x < Cr; x++) {
		for(uint z = 0; z < Cc; z++) {
			// We need to add the constraint: C_xz - sum_y A_xy R_yz = 0
			// C_xz has coeff 1, each R_yz has coeff - Axy, and the bounds are 0
			//
			auto con = qp.make_con(0, 0);
			qp.set_con_coeff(con, C[x][z], 1);	// coeff of C_xz

			for(uint y = 0; y < Rr; y++)
				qp.set_con_coeff(con, RR[y][z], -A(x,y));	// coeff of R_yz
		}
	}

	// constraints for R
	for(uint y = 0; y < Rr; y++) {
		for(uint z = 0; z < Rc; z++) {
			// We need to add the constraint: 0 <= R_yz <= 1
			// R_yz has coeff 1, and the bounds are 0,1
			//
			auto con = qp.make_con(0, 1);
			qp.set_con_coeff(con, RR[y][z], 1);
		}

		// We need to add the constraint: sum_z R_yz = 1
		// each R_yz has coeff 1, and the bounds are both 1
		auto con = qp.make_con(1, 1);

		for(uint z = 0; z < Rc; z++)
			qp.set_con_coeff(con, RR[y][z], 1);
	}

	// cost, quadratic part: C.C
	// coeff 2 for the C variables (because of the 1/2 in the QP), 0 for the R variables
	for(uint x = 0; x < Cr; x++)
		for(uint z = 0; z < Cc; z++)
			qp.set_obj_coeff(C[x][z], C[x][z], 2);

	// cost, linear part: -2BC
	// coeffs -2B for the C variables, 0 for the R variables
	for(uint x = 0; x < Cr; x++)
		for(uint z = 0; z < Cc; z++)
			qp.set_obj_coeff(C[x][z], -2 * B(x,z));

	// precision of the solution (Kostas: not 100% sure here)
	// if no value is set by the user, use 1e-5 instead of OSQP's defaults
	if(qp::Defaults::osqp_eps_abs < 0.0)	qp.osqp_eps_abs = 1e-5;
	if(qp::Defaults::osqp_eps_rel < 0.0)	qp.osqp_eps_rel = 1e-5;

	// ready
	if(!qp.solve())
		throw std::runtime_error("refined_by: QP infeasible, this shouldn't happen");

	// add B.B to the cost function to obtained the squared distance (see the program definition above)
	eT dist = qp.objective() + arma::dot(B, B);

	bool res = equal(dist, eT(0), eT(qp.osqp_eps_abs), eT(qp.osqp_eps_rel));
	if(res) {
		G.clear();
	} else {
		G.set_size(Cc, Cr);
		for(uint x = 0; x < Cr; x++)
			for(uint z = 0; z < Cc; z++)
				G(z,x) = B(x,z) - qp.solution(C[x][z]);
		G -= G.min();	// non-negative
		G /= G.max();	// and in [0,1]
	}

	R.set_size(Rr, Rc);
	for(uint y = 0; y < Rr; y++)
		for(uint z = 0; z < Rc; z++)
			R(y,z) = qp.solution(RR[y][z]);

	return res;
}

// same, without R
template<typename eT>
bool refined_by(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& G) {
	Chan<eT> R;
	return refined_by<eT>(A, B, G, R);
}



// true if A is max-case refined by B (i.e. A's max-case leakage is >= B's)
//
template<typename eT = eT_def>
inline
bool max_refined_by(const Chan<eT>& A, const Chan<eT>& B) {
	// true if X A* = B* for some X, where C* is produced from C by normalizing each column and transposing
	Mat<eT> An = A.t();
	Mat<eT> Bn = B.t();

	channel::normalize(An);
	channel::normalize(Bn);

	auto X = channel::left_factorize(Bn, An);

	return !X.empty();
}


// true if A is privacy-based refined by B (i.e. A sat. d-privact => B sat. d-privacy for all d)
//
template<typename eT = eT_def>
inline
bool priv_refined_by(const Chan<eT>& A, const Chan<eT>& B) {
	typedef eT R;		// use eT as the type for the metric result

	// This holds if B is d_A-private, where d_A is the metric "induced by" A:
	// d_A(x1, x2) = mtv(A.row(x1), A.row(x2))
	//
	auto d_A = metric::compose<R,uint,Prob<eT>>(
		metric::mult_total_variation<R, Prob<eT>>(),
		[&](uint x) -> Prob<eT> { return A.row(x); }
	);

	return mechanism::is_private(B, d_A);
}



template<typename eT>
eT add_metric(const Prob<eT>& pi, const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& G) {
	if(pi.n_elem != A.n_rows || A.n_rows != B.n_rows)
		throw std::runtime_error("invalid sizes");

	Mat<eT> AB = arma::join_rows(A, B);
	uint K = A.n_rows,
		 M = A.n_cols,
		 N = B.n_cols;

	assert(AB.n_rows == K && AB.n_cols == M+N);

	// vars: gain function (M+N)xK, bounded <= 1		(first for A, then for B)
	lp::LinearProgram<eT> lp;
	lp.maximize = true;
	auto vars = lp.make_vars(M+N, K, -infinity<eT>(), eT(1));

	// Objective function assumes each column gives a diffent best w
	// maximize
	// + sum_{M <= y < M+N} sum_x pi_x B_x,y g(y,x)
	// - sum_{0 <= y < M  } sum_x pi_x A_x,y g(y,x)
	//
	for(uint y = 0; y < M+N; y++)
		for(uint x = 0; x < K; x++)
			lp.set_obj_coeff(vars[y][x], (y < M ? -1 : 1) * pi(x) * AB(x,y));

	// s.t. 
	// sum_x pi_x A_x,y g(y,x)  >=  sum_x pi_x A_x,y g(w,x)       for all 0 <= y < M, w != y
	//
	for(uint y = 0; y < M; y++) {
		for(uint w = 0; w < M+N; w++) {
			if(y == w) continue;

			auto con = lp.make_con(0, infinity<eT>());

			for(uint x = 0; x < K; x++) {
				lp.set_con_coeff(con, vars[y][x],   pi(x) * AB(x,y));
				lp.set_con_coeff(con, vars[w][x], - pi(x) * AB(x,y));
			}
		}
	}

	// sum_x pi_x A_x,y g(y,x)  >=  0       for all 0 <= y < M
	//
	for(uint y = 0; y < M; y++) {
		auto con = lp.make_con(0, infinity<eT>());

		for(uint x = 0; x < K; x++)
			lp.set_con_coeff(con, vars[y][x], pi(x) * AB(x,y));
	}

	// ready
	if(!lp.solve())
		throw std::runtime_error("add_metric: LP infeasible, this shouldn't happen");

	// reconstrict gain function
	G.set_size(M+N+1, K);
	G.row(M+N).fill(0);

	for(uint w = 0; w < M+N; w++)
		for(uint x = 0; x < K; x++)
			G(w,x) = lp.solution(vars[w][x]);

	return lp.objective();
}

// same, but not interested in G
template<typename eT>
eT add_metric(const Prob<eT>& pi, const Chan<eT>& A, const Chan<eT>& B) {
	Mat<eT> G;
	return add_metric(pi, A, B, G);
}


// bound on the additive refinement metric for 1-bounded gain functions via the Kantorovich
//
template<typename eT>
eT add_metric_bound(const Prob<eT>& pi, const Chan<eT>& A, const Chan<eT>& B) {
	if(A.n_rows != B.n_rows)
		throw std::runtime_error("invalid sizes");

	Mat<eT> innersA, innersB;
	auto outerA = channel::hyper(A, pi, innersA);
	auto outerB = channel::hyper(B, pi, innersB);

	// We need the distributions outerA, outerB to refer to the same 'space' of inners.
	// So we put the inners together in the same matrix, and add zeroes on the right/left
	// of outerA/outerB respectively.
	//
	Mat<eT> inners = arma::join_rows(innersA, innersB);

	outerA.insert_cols(outerA.n_elem, innersB.n_cols);	// add on the right as many zeroes as B's inners
	outerB.insert_cols(0,             innersA.n_cols);	// add on the left  as many zeroes as A's inners

	// metric between the colums of A
	auto q = metric::compose<eT,uint,Prob<eT>>(
		metric::convex_separation_quasi<eT, Prob<eT>>(),

		[&](uint i) -> Prob<eT> { return inners.col(i).t(); }
	);

	// Important: we need kantorovich_lp, cause kantorovich uses the FastEMD algorithm which assumes
	// that q is a metric (but our convex_separation_quasi is not!)
	// TODO: does it work if we 'metricify' q?
	//
	auto kant = metric::kantorovich_lp<eT,Prob<eT>>(q);
	return kant(outerA, outerB);
}


}