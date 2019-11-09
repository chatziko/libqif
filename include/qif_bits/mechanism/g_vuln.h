
namespace mechanism::g_vuln {

// Returns the mechanism having the smallest expected loss (wrt pi, loss) given the Vg(pi,C) <= max_vuln constraint
//
template<typename eT>
Chan<eT> min_loss_given_max_vuln(
	const Prob<eT>& pi,
	uint n_cols,
	uint n_guesses,
	eT max_vuln,
	Metric<eT, uint> gain,
	Metric<eT, uint> loss
) {
	uint M = pi.n_cols,
		 N = n_cols;

	// C: M x N   unknowns
	// We have M x N variables
	lp::LinearProgram<eT> lp;
	auto vars = lp.make_vars(M, N, 0, 1);

	// cost function: minimize sum_xy pi_x C_xy loss(x,y)
	lp.maximize = false;
	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			lp.set_obj_coeff(vars[x][y], pi(x) * loss(x, y));

	// Build equations for Vg(pi, C) = sum_y max_w sum_x pi_x C_x,y g(w,x) <= max_vuln
	//
	// For this we use auxiliary variables:
	//    vuln_y = max_w sum_x pi_x C_x,y g(w,x) <= max_vuln
	//
	auto vuln_y = lp.make_vars(n_cols);

	// For each variable to really represent the max_w ..., we need to set constraints:
	//    vuln_y >= sum_x pi_x C_x,y g(w,x)     for each y,w
	//
	for(uint y = 0; y < N; y++) {
		for(uint w = 0; w < n_guesses; w++) {
			auto con = lp.make_con(-infinity<eT>(), 0);
			lp.set_con_coeff(con, vuln_y[y], -1);

			for(uint x = 0; x < M; x++)
				lp.set_con_coeff(con, vars[x][y], pi(x) * gain(w,x));
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
			C(x, y) = lp.solution(vars[x][y]);

	return C;
}

} // namespace mechanism::g_vuln
