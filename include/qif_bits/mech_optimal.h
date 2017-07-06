
namespace mechanism {

template <typename eT> using ME = lp::MatrixEntry<eT>;

// Returns the mechanism satisfying eps*d privacy and having the best utility wrt pi and loss
//
template<typename eT>
Chan<eT> optimal_utility(const Prob<eT>& pi, uint n_cols, Metric<eT, uint> d_priv, Metric<eT, uint> loss) {

	uint M = pi.n_cols,
		 N = n_cols,
		 n_adj = 0;
	eT inf = infinity<eT>();

	// find how many adjacent elements do we have
	for(uint x1 = 0; x1 < M; x1++)
	for(uint x2 = x1+1; x2 < M; x2++)
		if(d_priv.is_adjacent(x1, x2) && less_than(d_priv(x1, x2), inf))
			n_adj++;

	// C: M x N   unknowns
	// We have M x N variables, that will be unfolded in a vector.
	// The varialbe C[x,y] will have variable number xN+y.
	//
	uint n_vars = M * N,					// one var for each element of C
		 n_cons = 2*n_adj*N+M,				// one constraint for each C_xy, C_x'y for x adj x', plus M sum=1 constraints
		 n_cons_elems = 4*n_adj*N+M*N;		// 2 elems for each DP constraint + M*N elements for the sum=1 constraints

	lp::LinearProgram<eT> lp;
	lp.b.set_size(n_cons);
	lp.sense.set_size(n_cons);

	// cost function: minimize sum_xy pi_x C_xy loss(x,y)
	lp.maximize = false;
	lp.c = Col<eT>(n_vars);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			lp.c(x*N+y) = pi(x) * loss(x, y);

	std::list<ME<eT>> entries;
	uint cons_i = 0;

	// Build equations for C_xy <= exp(eps d_priv(x,x')) C_x'y
	//
	for(uint x1 = 0; x1 < M; x1++) {
	for(uint x2 = 0; x2 < M; x2++) {
		if(x1 == x2 || !d_priv.is_adjacent(x1, x2)) continue;		// constraints for non-adjacent inputs are redundant
		if(!less_than(d_priv(x1, x2), inf)) continue;				// inf distance, i.e. no constraint
		for(uint y = 0; y < N; y++) {

			lp.sense(cons_i) = '<';
			lp.b(cons_i) = eT(0);

			entries.push_back(ME<eT>(cons_i, x1*N+y, eT(1)));
			entries.push_back(ME<eT>(cons_i, x2*N+y, - std::exp(d_priv(x1, x2))));

			cons_i++;
		}
	}}

	// equalities for summing up to 1
	//
	for(uint x = 0; x < M; x++) {
		lp.b(cons_i) = eT(1);					// sum of row = 1
		lp.sense(cons_i) = '=';

		for(uint y = 0; y < N; y++)
			// coeff 1 for variable C[x,y]
			entries.push_back(ME<eT>(cons_i, x*N+y, eT(1)));

		cons_i++;
	}

	assert(cons_i == n_cons);					// added all constraints
	n_cons_elems *= 1;							// avoid unused warning
	assert(entries.size() == n_cons_elems);		// added all constraint elements

	lp.fill_A(entries);

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> C(M, N);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			C(x, y) = lp.x(x*N+y);

	return C;
}

template<typename eT>
Chan<eT> dist_optimal_utility(Prob<eT> pi, uint n_cols, Metric<eT, uint> d_priv, Metric<eT, uint> loss) {

	uint M = pi.n_cols,
		 N = n_cols;

	// insert all distances in a std::set to keep unique ones.
	// Then collect in dists vector and sort
	//
	std::set<eT> dists_set;
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			dists_set.insert(d_priv(x, y));

	uint D = dists_set.size();
	Col<eT> dists(D);
	uint i = 0;
	for(eT v : dists_set)
		dists(i++) = v;
	dists = arma::sort(dists);

	// DI : MxN matrix, DIxy is the index of d_priv(x,y) in dists
	Mat<uint> DI(M, N);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			DI(x, y) = arma::find(dists == d_priv(x, y), 1).eval()(0);

	// variable X[d,y] will have index d*N+y
	//
	uint n_vars = D*N,						// one var for each distinct distance and output
		 n_cons = 2*(D-1)*N+M,				// 2 constraints for each consecutive distances and output, plus M sum=1 constraints
		 n_cons_elems = 4*(D-1)*N+M*N;		// 2 elems for each DP constraint + M*N elements for the sum=1 constraints

	lp::LinearProgram<eT> lp;
	lp.b.set_size(n_cons);
	lp.sense.set_size(n_cons);

	// cost function: minimize sum_xy pi_x X[d(x,y),y] loss(x,y)
	lp.maximize = false;
	lp.c = Col<eT>(n_vars, arma::fill::zeros);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			lp.c(DI(x,y)*N+y) += pi(x) * loss(x, y);

	std::list<ME<eT>> entries;
	uint cons_i = 0;

	// Build equations for X_d_i <= exp(eps |d_i - d_i+1|) X_d_j
	//
	for(uint dist_i = 0; dist_i < D-1; dist_i++) {
		eT diff = dists(dist_i+1) - dists(dist_i);

		for(uint y = 0; y < N; y++) {
			lp.sense(cons_i) = '<';
			lp.b(cons_i) = eT(0);

			entries.push_back(ME<eT>(cons_i, dist_i*N+y, eT(1)));
			entries.push_back(ME<eT>(cons_i, (dist_i+1)*N+y, - std::exp(diff)));

			cons_i++;

			lp.sense(cons_i) = '<';
			lp.b(cons_i) = eT(0);

			entries.push_back(ME<eT>(cons_i, (dist_i+1)*N+y, eT(1)));
			entries.push_back(ME<eT>(cons_i, dist_i*N+y, - std::exp(diff)));

			cons_i++;
		}
	}

	// equalities for summing up to 1
	// Note: if the same distance appears multiple times in the same row, we're going to insert multiple 1's
	// for the same cell, which are summed due to the "true" param in the batch-insert below
	//
	for(uint x = 0; x < M; x++) {
		lp.b(cons_i) = eT(1);					// sum of row = 1
		lp.sense(cons_i) = '=';

		for(uint y = 0; y < N; y++)
			entries.push_back(ME<eT>(cons_i, DI(x,y)*N+y, eT(1)));

		cons_i++;
	}

	assert(cons_i == n_cons);					// added all constraints
	n_cons_elems *= 1;							// avoid unused warning
	assert(entries.size() == n_cons_elems);		// added all constraint elements

	lp.fill_A(entries);

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> C(M, N);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			C(x, y) = lp.x(DI(x,y)*N+y);

	return C;
}

// This is the first version, that had variables X_d instead of X_d,y
//
template<typename eT>
Chan<eT> dist_optimal_utility_strict(Prob<eT> pi, uint n_cols, Metric<eT, uint> d_priv, Metric<eT, uint> loss) {

	uint M = pi.n_cols,
		 N = n_cols;

	// insert all distances in a std::set to keep unique ones.
	// Then collect in dists vector and sort
	//
	std::set<eT> dists_set;
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			dists_set.insert(d_priv(x, y));

	uint D = dists_set.size();
	Col<eT> dists(D);
	uint i = 0;
	for(eT v : dists_set)
		dists(i++) = v;
	dists = arma::sort(dists);

	// DI : MxN matrix, DIxy is the index of d_priv(x,y) in dists
	Mat<uint> DI(M, N);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			DI(x, y) = arma::find(dists == d_priv(x, y), 1).eval()(0);

	// C: D variables
	//
	uint n_vars = D,						// one var for each distinct distance
		 n_cons = 2*(D-1)+M,				// 2 constraints for each consecutive distances, plus M sum=1 constraints
		 n_cons_elems = 4*(D-1)+M*N;		// 2 elems for each DP constraint + M*N elements for the sum=1 constraints

	lp::LinearProgram<eT> lp;
	lp.b.set_size(n_cons);
	lp.sense.set_size(n_cons);

	// cost function: minimize sum_xy pi_x var<C_xy> loss(x,y)
	lp.maximize = false;
	lp.c = Col<eT>(n_vars, arma::fill::zeros);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			lp.c(DI(x,y)) += pi(x) * loss(x, y);

	std::list<ME<eT>> entries;
	uint cons_i = 0;

	// Build equations for X_d_i <= exp(eps |d_i - d_i+1|) X_d_j
	//
	for(uint dist_i = 0; dist_i < D-1; dist_i++) {
		eT diff = dists(dist_i+1) - dists(dist_i);

		lp.sense(cons_i) = '<';
		lp.b(cons_i) = eT(0);

		entries.push_back(ME<eT>(cons_i, dist_i, eT(1)));
		entries.push_back(ME<eT>(cons_i, dist_i+1, - std::exp(diff)));

		cons_i++;

		lp.sense(cons_i) = '<';
		lp.b(cons_i) = eT(0);

		entries.push_back(ME<eT>(cons_i, dist_i+1, eT(1)));
		entries.push_back(ME<eT>(cons_i, dist_i, - std::exp(diff)));

		cons_i++;
	}

	// equalities for summing up to 1
	// Note: if the same distance appears multiple times in the same row, we're going to insert multiple 1's
	// for the same cell, which are summed due to the "true" param in the batch-insert below
	//
	for(uint x = 0; x < M; x++) {
		lp.b(cons_i) = eT(1);					// sum of row = 1
		lp.sense(cons_i) = '=';

		for(uint y = 0; y < N; y++)
			entries.push_back(ME<eT>(cons_i, DI(x,y), eT(1)));

		cons_i++;
	}

	assert(cons_i == n_cons);					// added all constraints
	n_cons_elems *= 1;							// avoid unused warning
	assert(entries.size() == n_cons_elems);		// added all constraint elements

	lp.fill_A(entries, true);					// add duplicate entries

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> C(M, N);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			C(x, y) = lp.x(DI(x,y));

	return C;
}

} // namespace mechanism
