
// Mechanism design for g-vulnerability

namespace mechanism::g_vuln {

// Returns the mechanism having the smallest E[loss] given the Vg(pi,C) <= max_vuln constraint
//
template<typename eT>
Chan<eT> min_loss_given_max_vuln(
	const Prob<eT>& pi,
	uint n_cols,
	uint n_guesses,
	eT max_vuln,
	Metric<eT, uint> gain,
	Metric<eT, uint> loss,
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	uint M = pi.n_cols,
		 N = n_cols;

	// C: M x N   unknowns, except those with loss(x,y) > hard_max_loss
	// We have (at most) M x N variables
	lp::LinearProgram<eT> lp;

	std::vector< std::list<std::pair<uint,uint>> > vars(M);	// vars[x] is a list of <y, var>

	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			if(less_than_or_eq(loss(x, y), hard_max_loss))
				vars[x].push_back(std::pair(y, lp.make_var(0, 1)));

	// cost function: minimize sum_xy pi_x C_xy loss(x,y)
	lp.maximize = false;
	for(uint x = 0; x < M; x++)
		for(auto& [y, var] : vars[x])
			lp.set_obj_coeff(var, pi(x) * loss(x, y));

	// Build equations for Vg(pi, C) = sum_y max_w sum_x pi_x C_x,y g(w,x) <= max_vuln
	//
	// For this we use auxiliary variables:
	//    vuln_y = max_w sum_x pi_x C_x,y g(w,x) = p(y) Vg(delta^y)
	//
	auto vuln_y = lp.make_vars(n_cols, eT(0));		// p(y) Vg(delta^y) >= 0 so this is safe, and useful (see below)

	// the actual vulnerability constraint:
	//    sum_y vuln_y <= max_vuln
	//
	auto max_vuln_con = lp.make_con(-infinity<eT>(), max_vuln);
	for(uint y = 0; y < N; y++)
		lp.set_con_coeff(max_vuln_con, vuln_y[y], 1);

	// For each variable to really represent the max_w ..., we need to set constraints:
	//    vuln_y >= sum_x pi_x C_x,y g(w,x)     for each y,w
	//
	for(uint w = 0; w < n_guesses; w++) {
		// Store the constraints (y,w) for _this_ w. If for a specific combination (y,w) we have forall x:(C_xy = 0 OR pi_x g(w,x) = 0),
		// the RHS of the constraint is 0 so we can completely remove the constraint, since vuln_y >= 0 is set anyway! So we
		// store the constraints in a map, we might only have a few of them for this w.
		//
		std::unordered_map<uint,uint> cons;

		// IMPORTANT: optimize for the case where most pi_x*g(w,x) are zero! So we have w,x in the outer loops
		for(uint x = 0; x < M; x++) {
			eT g = pi(x) * gain(w, x); 
			if(equal(g, eT(0)))		// nothing to do if 0
				continue;

			for(auto& [y, var] : vars[x]) {
				auto& con = cons[y];
				if(!con) {									// 0 always means "new in the map", because "max_vuln_con" above has already taken the 0 constrain index
					con = lp.make_con(-infinity<eT>(), 0);
					lp.set_con_coeff(con, vuln_y[y], -1);
				}
				lp.set_con_coeff(con, var, g);
			}
		}
	}

	// equalities for summing up to 1
	//
	for(uint x = 0; x < M; x++) {
		auto con = lp.make_con(1, 1);

		// coeff 1 for variable C[x,y]
		for(auto& [y, var] : vars[x]) {
			(void)y; // avoid unused warning
			lp.set_con_coeff(con, var, 1);
		}
	}

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> C(M, N, arma::fill::zeros);
	for(uint x = 0; x < M; x++)
		for(auto& [y, var] : vars[x])
			C(x, y) = lp.solution(var);

	return C;
}

// Returns the mechanism having the smallest Vg[pi,C] given the E[loss] <= max_loss constraint
//
template<typename eT>
Chan<eT> min_vuln_given_max_loss(
	const Prob<eT>& pi,
	uint n_cols,
	uint n_guesses,
	eT max_loss,
	Metric<eT, uint> gain,
	Metric<eT, uint> loss,
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	uint M = pi.n_cols,
		 N = n_cols;

	// C: M x N   unknowns, except those with loss(x,y) > hard_max_loss
	// We have (at most) M x N variables
	lp::LinearProgram<eT> lp;

	std::vector< std::list<std::pair<uint,uint>> > vars(M);	// vars[x] is a list of <y, var>

	for(uint x = 0; x < M; x++)
		for(uint y = 0; y < N; y++)
			if(less_than_or_eq(loss(x, y), hard_max_loss))
				vars[x].push_back(std::pair(y, lp.make_var(0, 1)));

	// loss constraint: sum_xy pi_x C_xy loss(x,y) <= max_loss
	auto max_loss_con = lp.make_con(-infinity<eT>(), max_loss);
	for(uint x = 0; x < M; x++)
		for(auto& [y, var] : vars[x])
			lp.set_con_coeff(max_loss_con, var, pi(x) * loss(x, y));

	// Objective function:
	//    minimize sum_y vuln_y	
	// using auxiliary variables:
	//    vuln_y  =  max_w sum_x pi_x C_x,y g(w,x)  =  p(y) Vg(delta^y)
	//
	auto vuln_y = lp.make_vars(n_cols, eT(0));		// p(y) Vg(delta^y) >= 0 so this is safe, and useful (see below)

	lp.maximize = false;
	for(uint y = 0; y < N; y++)
		lp.set_obj_coeff(vuln_y[y], 1);

	// For each aux variable to really represent the max_w ..., we need to set constraints:
	//    vuln_y >= sum_x pi_x C_x,y g(w,x)     for each y,w
	//
	for(uint w = 0; w < n_guesses; w++) {
		// Store the constraints (y,w) for _this_ w. If for a specific combination (y,w) we have forall x:(C_xy = 0 OR pi_x g(w,x) = 0),
		// the RHS of the constraint is 0 so we can completely remove the constraint, since vuln_y >= 0 is set anyway! So we
		// store the constraints in a map, we might only have a few of them for this w.
		//
		std::unordered_map<uint,uint> cons;

		// IMPORTANT: optimize for the case where most pi_x*g(w,x) are zero! So we have w,x in the outer loops
		for(uint x = 0; x < M; x++) {
			eT g = pi(x) * gain(w, x); 
			if(equal(g, eT(0)))		// nothing to do if 0
				continue;

			for(auto& [y, var] : vars[x]) {
				auto& con = cons[y];
				if(!con) {									// 0 always means "new in the map", because "max_loss_con" above has already taken the 0 constrain index
					con = lp.make_con(-infinity<eT>(), 0);
					lp.set_con_coeff(con, vuln_y[y], -1);
				}
				lp.set_con_coeff(con, var, g);
			}
		}
	}

	// equalities for summing up to 1
	//
	for(uint x = 0; x < M; x++) {
		auto con = lp.make_con(1, 1);

		// coeff 1 for variable C[x,y]
		for(auto& [y, var] : vars[x]) {
			(void)y; // avoid unused warning
			lp.set_con_coeff(con, var, 1);
		}
	}

	// solve program
	//
	if(!lp.solve())
		return Chan<eT>();

	// reconstrict channel from solution
	//
	Chan<eT> C(M, N, arma::fill::zeros);
	for(uint x = 0; x < M; x++)
		for(auto& [y, var] : vars[x])
			C(x, y) = lp.solution(var);

	return C;
}

} // namespace mechanism::g_vuln
