
namespace refinement {

// true if A is refined by B (i.e. A's leakage is >= B's), that is if B = AR for some channel R
// method: see channel::factorize
//
template<typename eT = eT_def>
inline
bool refined_by(const Chan<eT>& A, const Chan<eT>& B, const std::string& method = "auto") {
	return !channel::factorize(B, A, false, method).empty();
}


// same as refined_by(A, B), also giving
// - when true, the channel R such that B = AR (up to the precision of the method)
// - when false, a counter-example G: a gain function (with values in [0,1]) such that A leaks strictly less than B for
//   the uniform prior. The counter-example is checked directly, so a false answer is always correct (up to rounding
//   errors for double).
// Only for methods giving a witness (see channel::factorize).
//
// Note: refined_by is equivalent to add_metric(uni, A, B) == 0, so this also gives an LP solution, but
//       it turns out to be much slower than the projection method (although the latter is done via QP).
//
template<typename eT>
bool refined_by(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& G, Chan<eT>& R, const std::string& method = "auto") {
	if(A.n_rows != B.n_rows)
		throw std::runtime_error("invalid sizes");

	Mat<eT> H;
	R = channel::factorize(B, A, H, false, method);
	if(!R.empty()) {
		G.clear();
		return true;
	}

	// H is a witness that B != AR for all R, so G = H^T is a gain function with V_G(B) > V_G(A) for the uniform
	// prior (see factorize_is_witness). Shifting and scaling preserves this, so we make it in [0,1].
	G = H.t();
	G -= G.min();
	G /= G.max();
	return false;
}



// true if A is max-case refined by B (i.e. A's max-case leakage is >= B's)
//
template<typename eT = eT_def>
inline
bool max_refined_by(const Chan<eT>& A, const Chan<eT>& B, const std::string& method = "auto") {
	// true if X A* = B* for some X, where C* is produced from C by normalizing each column and transposing
	Mat<eT> An = A.t();
	Mat<eT> Bn = B.t();

	channel::normalize(An);
	channel::normalize(Bn);

	auto X = channel::left_factorize(Bn, An, false, method);

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

	return measure::d_privacy::is_private(B, d_A);
}



template<typename eT>
std::pair<eT,Mat<eT>> add_metric(const Prob<eT>& pi, const Chan<eT>& A, const Chan<eT>& B) {
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

			auto con = lp.make_con(eT(0), infinity<eT>());

			for(uint x = 0; x < K; x++) {
				lp.set_con_coeff(con, vars[y][x],   pi(x) * AB(x,y));
				lp.set_con_coeff(con, vars[w][x], - pi(x) * AB(x,y));
			}
		}
	}

	// sum_x pi_x A_x,y g(y,x)  >=  0       for all 0 <= y < M
	//
	for(uint y = 0; y < M; y++) {
		auto con = lp.make_con(eT(0), infinity<eT>());

		for(uint x = 0; x < K; x++)
			lp.set_con_coeff(con, vars[y][x], pi(x) * AB(x,y));
	}

	// ready
	if(!lp.solve())
		throw std::runtime_error("add_metric: LP failed with status " + lp.status + ", this shouldn't happen");

	// reconstrict gain function
	Mat<eT> G(M+N+1, K);
	G.row(M+N).fill(eT(0));

	for(uint w = 0; w < M+N; w++)
		for(uint x = 0; x < K; x++)
			G(w,x) = lp.solution(vars[w][x]);

	return { lp.objective(), G };
}


// bound on the additive refinement metric for 1-bounded gain functions via the Kantorovich
//
template<typename eT>
eT add_metric_bound(const Prob<eT>& pi, const Chan<eT>& A, const Chan<eT>& B) {
	if(A.n_rows != B.n_rows)
		throw std::runtime_error("invalid sizes");

	auto [outerA, innersA] = channel::hyper(A, pi);
	auto [outerB, innersB] = channel::hyper(B, pi);

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