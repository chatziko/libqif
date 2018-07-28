
namespace channel {

template <typename eT> using ME = lp::MatrixEntry<eT>;

template<typename eT>
inline
Chan<eT>& identity(Chan<eT>& C) {
	if(!C.is_square()) throw std::runtime_error("not square");
	C.eye();
	return C;
}

template<typename eT>
inline
Chan<eT> identity(uint n) {
	Chan<eT> C(n, n);
	return identity(C);
}


template<typename eT>
inline
Chan<eT>& no_interference(Chan<eT>& C) {
	C.zeros();
	C.col(0).fill(eT(1));
	return C;
}

template<typename eT>
inline
Chan<eT> no_interference(uint n, uint cols = 1) {
	Chan<eT> C(n, cols);
	return no_interference(C);
}


template<typename eT>
inline
Chan<eT>& randu(Chan<eT>& C) {
	for(uint i = 0; i < C.n_rows; i++)
		C.row(i) = probab::randu<eT>(C.n_cols);

	return C;
}

template<typename eT>
inline
Chan<eT> randu(uint n) {
	Chan<eT> C(n, n);
	return randu(C);
}

template<typename eT>
inline
Chan<eT> randu(uint n, uint m) {
	Chan<eT> C(n, m);
	return randu(C);
}

template<typename eT>
inline
Chan<eT> deterministic(arma::ucolvec map, uint n_cols) {
	Chan<eT> C(map.n_rows, n_cols, arma::fill::zeros);
	for(uint i = 0; i < map.n_rows; i++)
		C(i, map(i)) = 1;
	return C;
}


template<typename eT>
inline
bool is_proper(const Chan<eT>& C, const eT& mrd = def_max_rel_diff<eT>()) {
	for(uint i = 0; i < C.n_rows; i++)
		if(!probab::is_proper<eT>(C.row(i), mrd))
			return false;

	return true;
}


template<typename eT>
inline
void check_proper(const Chan<eT>& C) {
	if(!is_proper(C))
		throw std::runtime_error("not a proper matrix");
}


template<typename eT>
inline
void check_prior_size(const Prob<eT>& pi, const Chan<eT>& C) {
	if(C.n_rows != pi.n_cols)
		throw std::runtime_error("invalid prior size");
}


template<typename eT>
inline bool equal(const Chan<eT>& A, const Chan<eT>& B, const eT& md = def_max_diff<eT>(), const eT& mrd = def_max_rel_diff<eT>()) {
	if(A.n_rows != B.n_rows || A.n_cols != B.n_cols)
		return false;

	for(uint i = 0; i < A.n_rows; i++)
		for(uint j = 0; j < A.n_cols; j++)
			if(!qif::equal(A.at(i, j), B.at(i, j), md, mrd))
				return false;

	return true;
}


// lexicographic order
template<typename eT>
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
template<typename eT>
inline
Prob<eT> posterior(const Chan<eT>& C, const Prob<eT>& pi, uint y) {
	return (arma::trans(C.col(y)) % pi) / arma::dot(C.col(y), pi);
}


// returns the hyper produced by C and pi
//
template<typename eT>
inline
Prob<eT> hyper(const Chan<eT>& C, const Prob<eT>& pi, Mat<eT>& inners) {
	Prob<eT> outer = pi * C;
	inners = C;

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
			outer(col) = 0;
		} else
			first = col;
	}

	// remove zero probability columns
	for(int i = outer.n_elem - 1; i >= 0; i--)	// inverse order to avoid changing the index
		if(qif::equal(outer(i), eT(0))) {
			inners.shed_col(i);
			outer.shed_col(i);
		}

	return outer;
}


// returns the prior estimate produced by C, given observed output distribution out
//
template<typename eT>
inline
uint bayesian_update(const Chan<eT>& C, const Prob<eT>& out, Prob<eT>& pi, eT max_diff = eT(1e-6), uint max_reps = 0) {
	eT almost_zero(1e-6);

	if(pi.is_empty())
		pi = probab::uniform<eT>(C.n_rows);

	if(C.n_rows != pi.n_cols || C.n_cols != out.n_cols)
		throw std::runtime_error("invalid sizes");

	for(uint count = 1; ; count++) {
		// out_cur[i] == 0 implies out[i] == 0. Since we want to have out[i]/out_cur[i] = 0
		// in such cases, we change change out_cur[i] to 1 to avoid NaNs.
		arma::Row<eT> out_cur = pi * C;
		out_cur.elem( find(out_cur < almost_zero) ).ones();

		Prob<eT> new_pi = pi % trans(C * trans(out / out_cur));

		pi -= new_pi;
		eT diff = qif::norm1(pi);
		pi = new_pi;

		if(diff <= max_diff || count == max_reps)
			return count;
	}
}


// Returns a channel X such that A = B X
//
template<typename eT>
inline
Chan<eT> factorize_lp(const Chan<eT>& A, const Chan<eT>& B) {
	// A: M x N
	// B: M x R
	// X: R x N   unknowns
	//
	uint M = A.n_rows,
		 N = A.n_cols,
		 R = B.n_cols,
		 n_vars = R * N,			// one var for each element of X
		 n_cons = M * N + R,		// one constraint for each element of A (A[i,j] = dot(B[i,:], X[: j]), plus one constraint for row of X (sum = 1)
		 n_cons_elems = (M+1)*N*R;	// M*N*R elements for the A = B X constrains and N*R elements for the sum=1 constraints

	if(B.n_rows != M)
		return Chan<eT>();

	lp::LinearProgram<eT> lp;
	lp.b.set_size(n_cons);
	lp.c = arma::zeros<Chan<eT>>(n_vars);			// we don't really care to optimize, so cost function = 0
	lp.sense.set_size(n_cons);
	lp.sense.fill('=');

	std::list<ME<eT>> entries;						// for batch entry into A

	// Build equations for A = B X
	// We have R x N variables, that will be unfolded in a vector.
	// The varialbe X[r,n] will have variable number rN+n.
	// For each element m,n of A we have an equation A[i,j] = dot(B[i,:], X[: j])
	//
	for(uint m = 0; m < M; m++) {
		for(uint n = 0; n < N; n++) {
			uint row = m*N+n;
			lp.b(row) = A(m, n);				// sum of row is A[m,n]

			for(uint r = 0; r < R; r++)
				// coeff B[m,r] for variable X[r,n]
				entries.push_back(ME<eT>(row, r*N+n, B(m, r)));
		}
	}

	// equalities for summing up to 1
	//
	for(uint r = 0; r < R; r++) {
		uint row = M*N+r;
		lp.b(row) = eT(1);						// sum of row = 1

		for(uint n = 0; n < N; n++)
			// coeff 1 for variable X[r,n]
			entries.push_back(ME<eT>(row, r*N+n, eT(1)));
	}

	n_cons_elems -= entries.size();
	assert(0 == n_cons_elems);				// added all constraint elements
	lp.fill_A(entries);

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> X(R, N);
	for(uint r = 0; r < R; r++)
		for(uint n = 0; n < N; n++)
			X(r, n) = lp.x(r*N+n);

	return X;
}

template<typename eT>
inline
Chan<eT>& project_to_simplex(Chan<eT>& C) {
	Prob<eT> temp(C.n_cols);
	for(uint i = 0; i < C.n_rows; i++) {
		temp = C.row(i);
		C.row(i) = probab::project_to_simplex(temp);
	}
	return C;
}


// factorize using a subgradient method.
// see: http://see.stanford.edu/materials/lsocoee364b/02-subgrad_method_notes.pdf
//
template<typename eT>
inline
Chan<eT> factorize_subgrad(const Chan<eT>& A, const Chan<eT>& B, const eT max_diff = 1e-4) {
	using arma::dot;

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
	std::ostream nullstream(0);				// temporarily disable
	ARMA_SET_CERR(nullstream);				// error messages
	arma::solve(X, B, A);
	ARMA_SET_CERR(std::cerr);

	if(!X.n_cols) return X;
	project_to_simplex(X);

	// G = max l2-norm of B's rows
	eT G(0);
	for(uint i = 0; i < M; i++)
		G = std::max(G, arma::norm(B.row(i), 2));

	const eT R = sqrt(2 * L);
		  	 //RG = R * G;
	const eT inf = infinity<eT>();

	Chan<eT> Z(M, N);
	Chan<eT> S = arma::zeros<Chan<eT>>(L, N);

	uint k;
	eT min(1), bound;
	eT sum1(- R * R), sum2(0);

	for(k = 1; true; k++) {
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

		project_to_simplex(X);

		// compute the lower bound given by the stopping criterion in section 3.4 of the lecture notes.
		// if the bound is positive then the optimal is also positive, so the channel cannot be factorized
		//
		sum1 += alpha * (2 * f - alpha * g_norm_sq);
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

// Returns a channel X such that A = B X
//
template<typename eT>
inline
Chan<eT> factorize(const Chan<eT>& A, const Chan<eT>& B) {
	return factorize_subgrad(A, B);
}

template<>
inline
rchan factorize(const rchan& A, const rchan& B) {
	return factorize_lp(A, B);
}


template<typename eT>
inline void share_memory(Mat<eT>& A, Mat<eT>& B) {
	Mat<eT> temp(B.memptr(), B.n_rows, B.n_cols, false, false);

	A.reset();
	A.swap(temp);
}

// sum of column minima
//
template<typename eT>
eT sum_column_min(const Chan<eT>& C) {
	// this doesn't work for rat, investigae
	// return arma::accu(arma::min(C, 0));

	eT res(0);
	for(uint y = 0; y < C.n_cols; y++) {
		eT min(1);
		for(uint x = 0; x < C.n_rows; x++)
			if(C(x,y) < min)
				min = C(x,y);
		res += min;
	}
	return res;
}


// draw an input, then an output
template<typename eT>
inline
arma::Row<uint> draw(const Prob<eT>& pi, const Chan<eT>& C) {
	uint x = probab::draw<eT>(pi);
	uint y = probab::draw<eT>(C.row(x));
	return { x, y };
}

// efficient batch sampling via the joint distribution
template<typename eT>
inline
Mat<uint> draw(const Prob<eT>& pi, const Chan<eT>& C, uint n) {
	// build the joint	
	Mat<eT> J = C;
	J.each_col() %= pi.t();

	// draw from it (note that prbab::draw just needs an iterable object, and the drawn indexes are the positions in the iteration)
	Row<uint> drawn = probab::draw(J, n);

	Mat<uint> res(n, 2);
	for(uint i = 0; i < n; i++) {
		// J is iterated by probab::draw column-by-column. The index we get need to be converted to
		// indexes for i and j
		//
		uint joint_i = drawn(i);
		res(i, 0) = joint_i % J.n_rows;
		res(i, 1) = joint_i / J.n_rows;
	}
	return res;
}


} // namespace channel
