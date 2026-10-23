
namespace channel {

// normalizes all rows of C to create a channel
//
template<typename eT = eT_def>
inline
void normalize(Mat<eT>& C) {
	C.each_col() /= arma::sum(C, 1);
}

template<typename eT = eT_def>
inline
Chan<eT> normalize(const Mat<eT>& C) {
	Chan<eT> CC(C);
	normalize(CC);		// will call the in-place overload
	return CC;
}

template<typename eT = eT_def>
inline
void identity(Chan<eT>& C) {
	if(!C.is_square()) throw std::runtime_error("not square");
	C.eye();
}

template<typename eT = eT_def>
inline
Chan<eT> identity(uint n) {
	Chan<eT> C(n, n);
	identity(C);
	return C;		// separate return to allow move semantics!
}


template<typename eT = eT_def>
inline
void no_interference(Chan<eT>& C) {
	C.zeros();
	C.col(0).fill(eT(1));
}

template<typename eT = eT_def>
inline
Chan<eT> no_interference(uint n, uint cols = 1) {
	Chan<eT> C(n, cols);
	no_interference(C);
	return C;		// separate return to allow move semantics!
}


template<typename eT = eT_def>
inline
void randu(Chan<eT>& C) {
	for(uint i = 0; i < C.n_rows; i++)
		C.row(i) = probab::randu<eT>(C.n_cols);
}

template<typename eT = eT_def>
inline
Chan<eT> randu(uint n, uint m = 0) {
	if(m == 0) m = n;
	Chan<eT> C(n, m);
	randu(C);
	return C;		// separate return to allow move semantics!
}

template<typename eT = eT_def>
inline
Chan<eT> deterministic(arma::ucolvec map, uint n_cols) {
	Chan<eT> C(map.n_rows, n_cols, arma::fill::zeros);
	for(uint i = 0; i < map.n_rows; i++)
		C(i, map(i)) = 1;
	return C;
}

template<typename eT = eT_def>
inline
Chan<eT> deterministic(std::function<uint(uint)> map, uint n_rows, uint n_cols) {
	Chan<eT> C(n_rows, n_cols, arma::fill::zeros);
	for(uint i = 0; i < n_rows; i++)
		C(i, map(i)) = 1;
	return C;
}


template<typename eT = eT_def>
inline
bool is_proper(const Chan<eT>& C, const eT& mrd = def_mrd<eT>) {
	for(uint i = 0; i < C.n_rows; i++)
		if(!probab::is_proper<eT>(C.row(i), mrd))
			return false;

	return true;
}


template<typename eT = eT_def>
inline
void assert_proper(const Chan<eT>& C) {
	if(!is_proper(C))
		throw std::runtime_error("not a proper matrix");
}


template<typename eT = eT_def>
inline
void check_prior_size(const Prob<eT>& pi, const Chan<eT>& C) {
	if(C.n_rows != pi.n_cols)
		throw std::runtime_error("invalid prior size");
}


template<typename eT = eT_def>
inline bool equal(const Chan<eT>& A, const Chan<eT>& B, const eT& md = def_md<eT>, const eT& mrd = def_mrd<eT>) {
	if(A.n_rows != B.n_rows || A.n_cols != B.n_cols)
		return false;

	for(uint i = 0; i < A.n_rows; i++)
		for(uint j = 0; j < A.n_cols; j++)
			if(!qif::equal(A.at(i, j), B.at(i, j), md, mrd))
				return false;

	return true;
}


// lexicographic order
template<typename eT = eT_def>
int compare_columns(const Mat<eT>& A, uint j1, uint j2) {
	for(uint i = 0; i < A.n_rows; i++) {
		if(less_than(A(i, j1), A(i, j2)))
			return -1;
		else if(less_than(A(i, j2), A(i, j1)))
			return 1;
	}
	return 0;
}


// returns the posterior for a specific output y
//
template<typename eT = eT_def>
inline
Prob<eT> posterior(const Chan<eT>& C, const Prob<eT>& pi, uint y) {
	return (arma::trans(C.col(y)) % pi) / arma::dot(C.col(y), pi);
}


// returns all posteriors produced by C and pi. The returned matrix has the same
// size as C, with each column being a posterior
//
template<typename eT = eT_def>
inline
Mat<eT> posteriors(const Chan<eT>& C, const Prob<eT>& pi = {}) {
	Mat<eT>res = C;
	if(!pi.is_empty())					// if pi is not given it is assumed to be uniform, so no need to multiply
		res.each_col() %= pi.t();		// creates the joint
	res.each_row() /= arma::sum(res);	// normalizes each column into the posterior
	return res;
}



// returns the hyper produced by C and pi
//
#undef hyper	// MSVC adds this
template<typename eT = eT_def>
inline
std::pair<Prob<eT>,Mat<eT>> hyper(const Chan<eT>& C, const Prob<eT>& pi) {
	Prob<eT> outer = pi * C;
	Mat<eT> inners = C;

	// remove zero probability columns
	for(int i = outer.n_elem - 1; i >= 0; i--)	// inverse order to avoid changing the index
		if(qif::equal(outer(i), eT(0))) {
			inners.shed_col(i);
			outer.shed_col(i);
		}

	inners.each_col() %= pi.t();	// creates the joint
	inners.each_row() /= outer;		// normalizes each column into the posterior

	// remove linearly dependent cols and zero-probability cols
	// we first sort cols in lexicographic order. Then move the outer probability
	// of equal columns, and finally delete cols of zero probability
	//
	auto sorted = arma::linspace<arma::urowvec>(0, inners.n_cols-1, inners.n_cols);
    std::sort(sorted.begin(), sorted.end(), [&inners](uint a, uint b) {
		return compare_columns(inners, a, b) == -1;		// a < b
    });

	uint first = sorted(0);		// first col of a sequence of equal ones
	for(uint i = 1; i < sorted.n_elem; i++) {
		uint col = sorted(i);
		if(compare_columns(inners, first, col) == 0) {
			outer(first) += outer(col);
			outer(col) = eT(0);
		} else
			first = col;
	}

	// delete the zero-probability cols left by the merging (inverse order to avoid changing the index)
	for(int i = outer.n_elem - 1; i >= 0; i--)
		if(outer(i) == eT(0)) {
			inners.shed_col(i);
			outer.shed_col(i);
		}

	return { outer, inners };
}


// returns the reduced form of the channel
//
template<typename eT = eT_def>
inline
Chan<eT> reduced(const Chan<eT>& C) {
	// compute via a hyper constructed on a uniform prior
	auto [outer, R] = channel::hyper(C, probab::uniform<eT>(C.n_rows));
	R.each_row() %= outer;
	normalize(R);
	return R;
}


// returns the prior estimate produced by C, given observed output distribution out
//
template<typename eT = eT_def>
inline
std::pair<Prob<eT>, uint> iterative_bayesian_update(const Chan<eT>& C, const Prob<eT>& out, const Prob<eT>& start = {}, eT max_diff = eT(1e-6), uint max_reps = 0) {
	eT almost_zero(1e-6);

	Prob<eT> pi = start.is_empty() ? probab::uniform<eT>(C.n_rows) : start;

	if(C.n_rows != pi.n_cols || C.n_cols != out.n_cols)
		throw std::runtime_error("invalid sizes");

	for(uint count = 1; ; count++) {
		// out_cur[i] == 0 implies out[i] == 0. Since we want to have out[i]/out_cur[i] = 0
		// in such cases, we change change out_cur[i] to 1 to avoid NaNs.
		arma::Row<eT> out_cur = pi * C;

		if constexpr (std::is_same<eT, rat>::value) {
			// for rat we need to do it in two steps
			arma::Row<arma::uword> indices = out_cur < almost_zero;
			out_cur.elem( find(indices) ).ones();
		} else {
			out_cur.elem( find(out_cur < almost_zero) ).ones();
		}

		Prob<eT> new_pi = pi % trans(C * trans(out / out_cur));

		pi -= new_pi;
		eT diff = qif::norm1(pi);
		pi = new_pi;

		if(diff <= max_diff || count == max_reps)
			return { pi, count };
	}
}


// Factorization: find a channel X such that A = B X (or column-stochastic X if col_stoch == true), or show that
// none exists. The available methods are:
// - solve_lp:        solves the equations A = B X as a linear program
// - project_subgrad: projects A to { B X } wrt the max-norm, by a subgradient method (double only, precision ~1e-4)
// - project_lp:      projects A to { B X } wrt the manhattan norm, by linear programming. When no X exists, it also
//                    gives a witness H (see factorize_is_witness).
// - project_qp:      projects A to { B X } wrt the euclidean norm, by quadratic programming (double only). When no X
//                    exists, it also gives a witness H (see factorize_is_witness).


// Whether H (same size as A) proves that there is no X such that A = B X, i.e. whether <H,A> > max_X <H,BX>.
// For row-stochastic X this means V_G(A) > V_G(B) for the gain function G = H^T and the uniform prior.
// The check is exact for rat, and up to rounding errors for double.
//
template<typename eT = eT_def>
bool factorize_is_witness(Mat<eT> H, const Chan<eT>& A, const Chan<eT>& B, const bool col_stoch = false) {
	// scale to [-1,1], so that the comparison's tolerance is relative to H
	eT m = std::max(eT(H.max()), eT(-H.min()));
	if(m == eT(0))
		return false;
	H /= m;

	// max_X <H,BX> = sum_r max_n (B^T H)_rn for row-stochastic X, and sum_n max_r (B^T H)_rn for column-stochastic X
	// For row-stochastic X we use V_G(A) = sum_n max_n' (H^T A)_n'n >= <H,A> as the left-hand side, which is still
	// smaller than max_X <H,BX> when A = BX (post-processing cannot increase V_G).
	//
	eT lhs = col_stoch ? eT(arma::accu(H % A))           : eT(arma::accu(arma::max(H.t() * A)));
	eT rhs = col_stoch ? eT(arma::accu(arma::max(B.t() * H))) : eT(arma::accu(arma::max(H.t() * B)));
	return less_than(rhs, lhs);
}


// solve_lp method
//
template<typename eT = eT_def>
inline
Chan<eT> factorize_solve_lp(const Chan<eT>& A, const Chan<eT>& B, const bool col_stoch = false) {
	// A: M x N
	// B: M x R
	// X: R x N   unknowns
	//
	uint M = A.n_rows,
		 N = A.n_cols,
		 R = B.n_cols;

	if(B.n_rows != M)
		return Chan<eT>();

	// Build equations for A = B X
	// We have R x N variables
	lp::LinearProgram<eT> lp;
	auto vars = lp.make_vars(R, N, eT(0), eT(1));

	// For each element m,n of A we have an equation A[i,j] = dot(B[i,:], X[: j])
	//
	for(uint m = 0; m < M; m++) {
		for(uint n = 0; n < N; n++) {
			auto con = lp.make_con(A(m, n), A(m, n));	// sum of row is A[m,n]

			// coeff B[m,r] for variable X[r,n]
			for(uint r = 0; r < R; r++)
				lp.set_con_coeff(con, vars[r][n], B(m, r));
		}
	}

	// equalities for rows/cols summing up to 1
	//
	if(col_stoch) {
		for(uint n = 0; n < N; n++) {
			auto con = lp.make_con(eT(1), eT(1));

			// coeff 1 for variable X[r,n]
			for(uint r = 0; r < R; r++)
				lp.set_con_coeff(con, vars[r][n], eT(1));
		}

	} else {
		for(uint r = 0; r < R; r++) {
			auto con = lp.make_con(eT(1), eT(1));

			// coeff 1 for variable X[r,n]
			for(uint n = 0; n < N; n++)
				lp.set_con_coeff(con, vars[r][n], eT(1));
		}
	}

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> X(R, N);
	for(uint r = 0; r < R; r++)
		for(uint n = 0; n < N; n++)
			X(r, n) = lp.solution(vars[r][n]);

	return X;
}

template<typename eT = eT_def>
inline
void _simplex_project(Chan<eT>& C, bool col_stoch) {
	if(col_stoch)
		for(uint j = 0; j < C.n_cols; j++)
			C.col(j) = metric::optimize::simplex_project((Prob<eT>)C.col(j).t()).t();
	else
		for(uint i = 0; i < C.n_rows; i++)
			C.row(i) = metric::optimize::simplex_project((Prob<eT>)C.row(i));
}


// project_subgrad method: minimizes max_{i,j} |(BX)(i,j) - A(i,j)| using a subgradient method, until either the
// minimum is (close to) 0, or a lower bound shows that it is positive.
// see: http://see.stanford.edu/materials/lsocoee364b/02-subgrad_method_notes.pdf
//
// It is usually fast, but can be very slow for some inputs (eg when A is close to, but not, factorizable). If
// max_iter > 0 and neither happens within max_iter iterations, an empty matrix is returned and *decided is set to false.
//
template<typename eT = eT_def>
inline
Chan<eT> factorize_project_subgrad(const Chan<eT>& A, const Chan<eT>& B, const bool col_stoch = false, const eT max_diff = 1e-4,
	const uint max_iter = 0, bool* decided = nullptr) {
	using arma::dot;

	if(decided)
		*decided = true;

	const bool debug = false;

	// A: M x N
	// B: M x L
	// X: L x N   unknowns
	//
	uint M = A.n_rows,
		 N = A.n_cols,
		 L = B.n_cols;
	Chan<eT> X;

	if(B.n_rows != M)
		return X;

	// Solve B * X = A, if no solution exists then A is not factorizable
	//
	auto cerr_buf = std::cerr.rdbuf(nullptr);	// temporarily disable armadillo's
	arma::solve(X, B, A);						// warnings (printed to std::cerr)
	std::cerr.rdbuf(cerr_buf);

	if(!X.n_cols) return X;
	_simplex_project(X, col_stoch);

	// G = max l2-norm of B's rows
	eT G(0);
	for(uint i = 0; i < M; i++)
		G = std::max(G, arma::norm(B.row(i), 2));

	const eT R = sqrt(2 * (col_stoch ? N : L));		// diameter of the set of row (or column) stochastic X
		  	 //RG = R * G;
	const eT inf = infinity<eT>();

	Chan<eT> Z(M, N);
	Chan<eT> S = arma::zeros<Chan<eT>>(L, N);

	uint k;
	eT min(1), bound;
	eT sum1(- R * R), sum2(0);

	for(k = 1; true; k++) {
		if(max_iter && k > max_iter) {
			X.clear();
			if(decided)
				*decided = false;
			break;
		}

		// compute
		//    f = max_{i,j} | (B*X)(i,j) - A(i,j) |
		//      = max_{i,j,sign} sign*((B*X)(i,j) - A(i,j))     (sign in {1,-1})
		//
		// f_i, f_j, sign are the ones that give the max
		//
		Z = B * X - A;
		int sign = 1;
		eT f(0);
		uint f_i = 0, f_j = 0;

		for(uint i = 0; i < M; i++) {
			for(uint j = 0; j < N; j++) {
				eT diff = std::abs(Z(i, j));
				if(diff > f) {
					f = diff;
					f_i = i;
					f_j = j;
					sign = Z(i, j) < eT(0) ? -1 : 1;
				}
			}
		}

		// update min if we found a better one. when we reach zero we're done
		//
		if(f < min) {
			min = f;

			if(qif::equal(min, eT(0), max_diff))
				break;
		}

		// update X, using the CFM method in section 8 of the lecture notes
		//
		Col<eT> g = eT(sign) * B.row(f_i).t();	// the subgradient (this is actually the j-th col of the subgradient, all other cols are 0)
		eT g_norm_sq = dot(g, g);

		eT beta = std::max(eT(0), - eT(1.5) * dot(S.col(f_j), g) / g_norm_sq);
		S *= beta;
		S.col(f_j) += g;

		// S might have arbitrarily large values, such that its norm might become inf
		// In this case we reset it to just g (no memory). Note: this sounds safe but we should check
		eT s_norm_sq = arma::dot(S, S);
		if(s_norm_sq == inf) {
			S.fill(0);
			S.col(f_j) = g;
			s_norm_sq = g_norm_sq;
		}

		// alternative method with fixed beta (Note: for larger values of beta we were getting wrong results, maybe the bound does not hold for this method?)
		// const eT beta = 0.25;		// a high memory value seems to work well
		// S *= beta;
		// S.col(f_j) += (1-beta) * g;

		eT alpha = f / s_norm_sq;
		X -= alpha * S;

		_simplex_project(X, col_stoch);

		// compute the lower bound given by the stopping criterion in section 3.4 of the lecture notes.
		// if the bound is positive then the optimal is also positive, so the channel cannot be factorized
		// Note: the step is alpha * S (not alpha * g), so its squared norm is alpha^2 |S|^2
		//
		sum1 += alpha * (2 * f - alpha * s_norm_sq);
		sum2 += 2 * alpha;
		bound = sum1 / sum2;

		if(!less_than_or_eq(bound, eT(0), max_diff)) {
			X.clear();
			break;
		}

		if(debug && k % 100 == 0)
			std::cout << "k: " << k << ", min: " << min << ", bound: " << bound << "\n";
	}

	if(debug)
		std::cout << "k: " << k << ", min: " << min << ", bound: " << bound << "\n";

	return X;
}

// project_lp method: computes the projection C = BX of A to { BX } wrt the manhattan norm (see Boyd 8.1.2). If C = A
// (exactly for rat, up to the solver's precision for double) X is returned. Otherwise an empty matrix is returned, and
// H is a witness showing that A != BX for all X (see factorize_is_witness).
//
// The projection is computed by the LP
//   min sum_mn E_mn   s.t.  -E <= A - BX <= E,  X row (or column) stochastic
// and the witness by its dual. Since |A-C|_1 = max_{|H| <= 1} <H,A-C>, the distance is
//   min_X max_{|H| <= 1} <H,A-BX> = max_{|H| <= 1} <H,A> - max_X <H,BX>
// and max_X <H,BX> = sum_r max_n (B^T H)_rn for row-stochastic X (sum_n max_r (B^T H)_rn for column-stochastic X).
// So the dual LP is
//   max <H,A> - sum t   s.t.  t_r >= (B^T H)_rn for all r,n  (t_n for column-stochastic X),  -1 <= H <= 1
// and when the optimum (the distance) is positive, H is a witness.
//
template<typename eT = eT_def>
Chan<eT> factorize_project_lp(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& H, const bool col_stoch = false) {
	// A: M x N
	// B: M x R
	// X: R x N   unknowns
	//
	uint M = A.n_rows,
		 N = A.n_cols,
		 R = B.n_cols;

	H.clear();
	if(B.n_rows != M)
		return Chan<eT>();

	const eT inf = infinity<eT>();
	Chan<eT> res(R, N);

	// primal: the projection
	{
		lp::LinearProgram<eT> lp;
		lp.maximize = false;

		auto X = lp.make_vars(R, N, eT(0), eT(1));
		auto E = lp.make_vars(M, N, eT(0), inf);

		for(uint m = 0; m < M; m++) {
			for(uint n = 0; n < N; n++) {
				lp.set_obj_coeff(E[m][n], eT(1));

				// E_mn - (BX)_mn >= -A_mn  and  E_mn + (BX)_mn >= A_mn
				auto con1 = lp.make_con(-A(m,n), inf);
				auto con2 = lp.make_con( A(m,n), inf);
				lp.set_con_coeff(con1, E[m][n], eT(1));
				lp.set_con_coeff(con2, E[m][n], eT(1));

				for(uint r = 0; r < R; r++) {
					lp.set_con_coeff(con1, X[r][n], -B(m,r));
					lp.set_con_coeff(con2, X[r][n],  B(m,r));
				}
			}
		}

		// rows (or columns) of X sum up to 1
		if(col_stoch) {
			for(uint n = 0; n < N; n++) {
				auto con = lp.make_con(eT(1), eT(1));
				for(uint r = 0; r < R; r++)
					lp.set_con_coeff(con, X[r][n], eT(1));
			}
		} else {
			for(uint r = 0; r < R; r++) {
				auto con = lp.make_con(eT(1), eT(1));
				for(uint n = 0; n < N; n++)
					lp.set_con_coeff(con, X[r][n], eT(1));
			}
		}

		if(!lp.solve())
			throw std::runtime_error("factorize_project_lp: LP failed (" + lp.status + "), this shouldn't happen");

		for(uint r = 0; r < R; r++)
			for(uint n = 0; n < N; n++)
				res(r,n) = lp.solution(X[r][n]);

		if(qif::equal(lp.objective(), eT(0)))
			return res;
	}

	// dual: the witness
	lp::LinearProgram<eT> lp;
	lp.maximize = true;

	// t is a max of entries of B^T H with |H| <= 1, so |t| <= the max column sum of B. Bounding t doesn't change the
	// solution, but helps the solver (SoPlex fails on some rat programs with free t)
	eT bmax = R > 0 ? eT(arma::max(arma::sum(B, 0))) : eT(0);
	auto HH = lp.make_vars(M, N, eT(-1), eT(1));
	auto t = lp.make_vars(col_stoch ? N : R, -bmax, bmax);

	for(uint m = 0; m < M; m++)
		for(uint n = 0; n < N; n++)
			lp.set_obj_coeff(HH[m][n], A(m,n));
	for(auto v : t)
		lp.set_obj_coeff(v, eT(-1));

	// t_r - (B^T H)_rn >= 0  (t_n for column-stochastic X)
	for(uint r = 0; r < R; r++) {
		for(uint n = 0; n < N; n++) {
			auto con = lp.make_con(eT(0), inf);
			lp.set_con_coeff(con, t[col_stoch ? n : r], eT(1));
			for(uint m = 0; m < M; m++)
				lp.set_con_coeff(con, HH[m][n], -B(m,r));
		}
	}

	if(!lp.solve())
		throw std::runtime_error("factorize_project_lp: dual LP failed (" + lp.status + "), this shouldn't happen");

	H.set_size(M, N);
	for(uint m = 0; m < M; m++)
		for(uint n = 0; n < N; n++)
			H(m,n) = lp.solution(HH[m][n]);

	if(factorize_is_witness(H, A, B, col_stoch))
		return Chan<eT>();

	// For rat the dual optimum is the (positive) distance, so H is always a witness. For double, the check can fail
	// due to rounding errors when the distance is tiny, in which case we consider A factorizable (up to precision).
	if constexpr (std::is_same<eT, rat>::value)
		throw std::runtime_error("factorize_project_lp: dual solution is not a witness, this shouldn't happen");

	H.clear();
	return res;
}


// project_qp method: computes the projection C = BX of A to { BX } wrt the euclidean norm. If C = A (up to the
// solver's precision) X is returned. Otherwise an empty matrix is returned, and H is a witness showing that A != BX
// for all X (see factorize_is_witness).
//
// For the projection C*, H = A-C* is always a witness when C* != A (by weak duality <H,A> - max_X <H,BX> >= |H|^2),
// and the iterates A-C_k converge to it. So we check H = A-C_k periodically while solving, and stop as soon as it is
// a witness, without waiting for convergence.
//
template<typename eT = eT_def>
Chan<eT> factorize_project_qp(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& H, const bool col_stoch = false) {
	// A: M x N
	// B: M x R
	// X: R x N   unknowns
	//
	uint M = A.n_rows,
		 N = A.n_cols,
		 R = B.n_cols;

	H.clear();
	if(B.n_rows != M)
		return Chan<eT>();

	// The (squared) euclidean distance between C and A can be written as ("." is the dot product)
	//   |C-A|_2^2 = (C-A).(C-A) = C.C -2A.C + A.A
	// A.A is constant so we need to minimize C.C -2A.C
	//
	// We have 2 types of variables:
	// - entries of C (M x N), constrained by C = BX
	// - entries of X (R x N), constrained to be row (or column) stochastic
	//
	qp::QuadraticProgram<eT> qp;

	auto C = qp.make_vars(M, N);
	auto X = qp.make_vars(R, N);

	// constraints for C
	for(uint m = 0; m < M; m++) {
		for(uint n = 0; n < N; n++) {
			// We need to add the constraint: C_mn - sum_r B_mr X_rn = 0
			// C_mn has coeff 1, each X_rn has coeff - B_mr, and the bounds are 0
			//
			auto con = qp.make_con(0, 0);
			qp.set_con_coeff(con, C[m][n], 1);	// coeff of C_mn

			for(uint r = 0; r < R; r++)
				qp.set_con_coeff(con, X[r][n], -B(m,r));	// coeff of X_rn
		}
	}

	// constraints for X: 0 <= X_rn <= 1
	for(uint r = 0; r < R; r++) {
		for(uint n = 0; n < N; n++) {
			auto con = qp.make_con(0, 1);
			qp.set_con_coeff(con, X[r][n], 1);
		}
	}

	// rows (or columns) of X sum up to 1
	if(col_stoch) {
		for(uint n = 0; n < N; n++) {
			auto con = qp.make_con(1, 1);
			for(uint r = 0; r < R; r++)
				qp.set_con_coeff(con, X[r][n], 1);
		}
	} else {
		for(uint r = 0; r < R; r++) {
			auto con = qp.make_con(1, 1);
			for(uint n = 0; n < N; n++)
				qp.set_con_coeff(con, X[r][n], 1);
		}
	}

	// cost, quadratic part: C.C
	// coeff 2 for the C variables (because of the 1/2 in the QP), 0 for the X variables
	for(uint m = 0; m < M; m++)
		for(uint n = 0; n < N; n++)
			qp.set_obj_coeff(C[m][n], C[m][n], 2);

	// cost, linear part: -2A.C
	// coeffs -2A for the C variables, 0 for the X variables
	for(uint m = 0; m < M; m++)
		for(uint n = 0; n < N; n++)
			qp.set_obj_coeff(C[m][n], -2 * A(m,n));

	// precision of the solution (Kostas: not 100% sure here)
	// if no value is set by the user, use 1e-5 instead of OSQP's defaults
	if(qp::Defaults::osqp_eps_abs < 0.0)	qp.osqp_eps_abs = 1e-5;
	if(qp::Defaults::osqp_eps_rel < 0.0)	qp.osqp_eps_rel = 1e-5;

	// H = A - C for the current solution, true if it is a witness
	auto is_witness = [&](const Col<eT>& sol) -> bool {
		H.set_size(M, N);
		for(uint m = 0; m < M; m++)
			for(uint n = 0; n < N; n++)
				H(m,n) = A(m,n) - sol(C[m][n]);
		return factorize_is_witness(H, A, B, col_stoch);
	};
	qp.early_stop = is_witness;

	qp.solve();

	if(qp.status == qp::Status::STOPPED)
		return Chan<eT>();		// H already set by is_witness
	if(qp.status != qp::Status::OPTIMAL)
		throw std::runtime_error("factorize_project_qp: QP infeasible, this shouldn't happen");
	if(is_witness(qp.solution()))
		return Chan<eT>();

	H.clear();
	Chan<eT> res(R, N);
	for(uint r = 0; r < R; r++)
		for(uint n = 0; n < N; n++)
			res(r,n) = qp.solution(X[r][n]);
	return res;
}

template<typename eT = eT_def>
Chan<eT> _factorize(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>* H, const bool col_stoch, std::string method) {
	// deprecated names (1.5.0)
	if(method == "lp")		method = "solve_lp";
	if(method == "subgrad")	method = "project_subgrad";

	constexpr bool is_rat = std::is_same<eT, rat>::value;

	if(method != "auto" && method != "solve_lp" && method != "project_subgrad" && method != "project_lp" && method != "project_qp")
		throw std::runtime_error("invalid method: " + method);
	if(is_rat && (method == "project_subgrad" || method == "project_qp"))
		throw std::runtime_error(method + " method not available for rat");
	if(H && (method == "solve_lp" || method == "project_subgrad"))
		throw std::runtime_error(method + " method does not give a witness");

	Mat<eT> H_;
	if(method == "project_lp")
		return factorize_project_lp(A, B, H ? *H : H_, col_stoch);

	if constexpr (is_rat) {
		// auto or solve_lp
		return H ? factorize_project_lp(A, B, *H, col_stoch) : factorize_solve_lp(A, B, col_stoch);

	} else {
		if(method == "auto") {
			if(H) {
				method = "project_qp";		// faster than project_lp
			} else {
				// The subgradient method is usually fast, but can be very slow for some inputs (eg when A is close to,
				// but not, factorizable). So we try it with a budget of iterations, and if it doesn't decide we fall back
				// to project_qp, or solve_lp if B has more columns than rows (the QP is slow for such B).
				bool decided;
				Chan<eT> X = factorize_project_subgrad(A, B, col_stoch, eT(1e-4), 10000, &decided);
				if(decided)
					return X;
				method = B.n_cols > B.n_rows ? "solve_lp" : "project_qp";
			}
		}

		if(method == "solve_lp")
			return factorize_solve_lp(A, B, col_stoch);
		if(method == "project_subgrad")
			return factorize_project_subgrad(A, B, col_stoch);

		return factorize_project_qp(A, B, H ? *H : H_, col_stoch);
	}
}

// Returns a channel X such that A = B X, or an empty matrix if none exists. If col_stoch == true then the
// returned X is column-stochastic.
// method: "solve_lp", "project_subgrad" (not for rat), "project_lp", "project_qp" (not for rat), or "auto": for rat
// solve_lp, for double project_subgrad with a budget of iterations, falling back to project_qp (or solve_lp, see
// _factorize).
// See factorize_* above.
//
template<typename eT = eT_def>
inline
Chan<eT> factorize(const Chan<eT>& A, const Chan<eT>& B, const bool col_stoch = false, const std::string& method = "auto") {
	return _factorize<eT>(A, B, nullptr, col_stoch, method);
}

// Same as factorize(A, B), also giving a witness H when no X exists (see factorize_is_witness).
// Only for methods giving a witness: "project_lp", "project_qp", and "auto" (project_qp for double, project_lp for rat).
//
template<typename eT = eT_def>
inline
Chan<eT> factorize(const Chan<eT>& A, const Chan<eT>& B, Mat<eT>& H, const bool col_stoch = false, const std::string& method = "auto") {
	return _factorize<eT>(A, B, &H, col_stoch, method);
}


// Returns a channel X such that A = X B (method: see factorize)
//
template<typename eT = eT_def>
inline
Chan<eT> left_factorize(const Chan<eT>& A, const Chan<eT>& B, const bool col_stoch = false, const std::string& method = "auto") {
	// A = X B if A^t = B^t X^t, so we use factorize asking for an "inversly"-stochastic matrix
	//
	Chan<eT> X = factorize((Chan<eT>)A.t(), (Chan<eT>)B.t(), !col_stoch, method);
	arma::inplace_trans(X);
	return X;
}


// sum of column minima
//
template<typename eT = eT_def>
eT sum_column_min(const Chan<eT>& C) {
	return arma::accu(arma::min(C, 0));
}

// sample an input, then an output
template<typename eT = eT_def>
inline
std::pair<uint,uint> sample(const Chan<eT>& C, const Prob<eT>& pi) {
	uint x = probab::sample<eT>(pi);
	uint y = probab::sample<eT>(C.row(x));
	return { x, y };
}

// efficient batch sampling via the joint distribution
template<typename eT = eT_def>
inline
Mat<uint> sample(const Chan<eT>& C, const Prob<eT>& pi, uint n) {
	// build the joint	
	Mat<eT> J = C;
	J.each_col() %= pi.t();

	// sample from it (note that prbab::sample just needs an iterable object, and the sampled indexes are the positions in the iteration)
	Row<uint> sampled = probab::sample(J, n);

	Mat<uint> res(n, 2);
	for(uint i = 0; i < n; i++) {
		// J is iterated by probab::sample column-by-column. The index we get need to be converted to
		// indexes for i and j
		//
		uint joint_i = sampled(i);
		res(i, 0) = joint_i % J.n_rows;
		res(i, 1) = joint_i / J.n_rows;
	}
	return res;
}


} // namespace channel
