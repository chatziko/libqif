
// Mechanism design for Bayes vulnerability

namespace mechanism::bayes_vuln {

// Returns the mechanism having the smallest expected loss (wrt pi, loss) given the V(pi,C) <= max_vuln constraint.
// Same as mechanism::g_vuln::min_loss_given_max_vuln for the identity gain function, but faster to construct the LP.
//
template<typename eT>
Chan<eT> min_loss_given_max_vuln(
	const Prob<eT>& pi,
	uint n_cols,
	eT max_vuln,
	Metric<eT, uint> loss,
	eT max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > max_loss
) {
	uint M = pi.n_cols,
		 N = n_cols;

	// C: M x N   unknowns
	// We have M x N variables
	lp::LinearProgram<eT> lp;

	uint zero_var = infinity<uint>();										// these vars are forced to 0
	std::vector<std::vector<uint>> vars(M, std::vector<uint>(N, zero_var));	// MxN vector of vars

	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			if(less_than_or_eq(loss(x, y), max_loss))
				vars[x][y] = lp.make_var(0, 1);

	// cost function: minimize sum_xy pi_x C_xy loss(x,y)
	lp.maximize = false;
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			if(vars[x][y] != zero_var)
				lp.set_obj_coeff(vars[x][y], pi(x) * loss(x, y));

	// Build equations for V(pi, C) = sum_y max_x pi_x C_x,y <= max_vuln
	//
	// For this we use auxiliary variables:
	//    vuln_y = max_x pi_x C_x,y <= max_vuln
	//
	auto vuln_y = lp.make_vars(n_cols);

	// For each variable to really represent the max_x ..., we need to set constraints:
	//    vuln_y >= pi_x C_x,y     for each y,x
	//
	for(uint y = 0; y < N; y++) {
		for(uint x = 0; x < M; x++) {
			if(vars[x][y] != zero_var) {
				auto con = lp.make_con(-infinity<eT>(), 0);
				lp.set_con_coeff(con, vuln_y[y], -1);
				lp.set_con_coeff(con, vars[x][y], pi(x));
			}
		}
	}

	// finally, the actual vulnerability constraint:
	//    sum_y vuln_y <= max_vuln
	//
	auto max_vuln_con = lp.make_con(-infinity<eT>(), max_vuln);
	for(uint y = 0; y < N; y++)
		lp.set_con_coeff(max_vuln_con, vuln_y[y], 1);

	// equalities for summing up to 1
	//
	for(uint x = 0; x < M; x++) {
		auto con = lp.make_con(1, 1);

		// coeff 1 for variable C[x,y]
		for(uint y = 0; y < N; y++)
			if(vars[x][y] != zero_var)
				lp.set_con_coeff(con, vars[x][y], 1);
	}

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> C(M, N);
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			C(x, y) = vars[x][y] == zero_var ? eT(0) : lp.solution(vars[x][y]);

	return C;
}

// Treating C's row-th row as variables, computes the row that minimizes
// the channel's posterior vulnerability. The row is updated inside C and the
// optimal vulnerability is returned
//
template<typename eT = eT_def>
inline
eT min_vuln_for_row(const Prob<eT>& pi, Chan<eT>& C, uint row) {
	uint N = C.n_cols;
	eT inf = infinity<eT>();

	// q: N variables
	lp::LinearProgram<eT> lp;
	auto vars = lp.make_vars(N, eT(0), eT(1));

	// We set z_y >= pi[x] C[x,y]   for each x
	// Since all rows but 'row' are fixed, this simply means
	//   z_y >= max_{x != r} pi[x] C[x,y]  and
	//   z_y >= pi[row] C[row,y]  and
	//
	C.row(row).fill(0);
	arma::Row<eT> maxes = arma::max(C.each_col() % pi.t(), 0);

	auto z = lp.make_vars(N, eT(0), eT(1));

	for(uint y = 0; y < N; y++) {
		auto con = lp.make_con(maxes(y), inf);
		lp.set_con_coeff(con, z[y], eT(1));

		con = lp.make_con(0, inf);
		lp.set_con_coeff(con, z[y], eT(1));
		lp.set_con_coeff(con, vars[y], -pi(row));
	}

	// cost function: minimize sum_y z_y
	//
	lp.maximize = false;
	for(uint y = 0; y < N; y++)
		lp.set_obj_coeff(z[y], eT(1));

	// equalities for summing up to 1
	//
	auto con = lp.make_con(1, 1);
	for(uint y = 0; y < N; y++)
		lp.set_con_coeff(con, vars[y], 1);

	// solve program
	//
	if(!lp.solve())
		throw std::runtime_error("min_l1_enclosing_ball: lp should be always solvable");

	// reconstruct q from solution
	//
	for(uint y = 0; y < N; y++)
		C(row, y) = lp.solution(vars[y]);

	return lp.objective();
}


} // namespace mechanism::bayes_vuln
