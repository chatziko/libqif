
// Mechanism design for Bayes vulnerability

namespace mechanism::bayes_vuln {

using measure::g_vuln::g_id;

// Returns the mechanism having the smallest E[loss] given the V(pi,C) <= max_vuln constraint.
// Just calls g_vuln::min_loss_given_max_vuln for the identity gain function.
//
template<typename eT>
Chan<eT> min_loss_given_max_vuln(
	const Prob<eT>& pi,
	uint n_cols,
	eT max_vuln,
	Metric<eT, uint> loss,
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	// Just use the correspondig g_vuln function with g_id. It is optimized to take advantage of a "sparce" g.
	//
	return g_vuln::min_loss_given_max_vuln(pi, n_cols, pi.n_elem, max_vuln, g_id<eT>, loss, hard_max_loss);
}

// Returns the mechanism having the smallest Bayes vulnerabiliy given the E[loss] <= max_loss constraint.
// Same as mechanism::g_vuln::min_vuln_given_max_loss for the identity gain function, but faster to construct the LP.
//
template<typename eT>
Chan<eT> min_vuln_given_max_loss(
	const Prob<eT>& pi,
	uint n_cols,
	eT max_loss,
	Metric<eT, uint> loss,
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	// Just use the correspondig g_vuln function with g_id. It is optimized to take advantage of a "sparce" g.
	//
	return g_vuln::min_vuln_given_max_loss(pi, n_cols, pi.n_elem, max_loss, g_id<eT>, loss, hard_max_loss);
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
		throw std::runtime_error("min_vuln_for_row: lp should be always solvable");

	// reconstruct q from solution
	//
	for(uint y = 0; y < N; y++)
		C(row, y) = lp.solution(vars[y]);

	return lp.objective();
}


} // namespace mechanism::bayes_vuln
