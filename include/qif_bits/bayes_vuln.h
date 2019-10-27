#define _USE_MATH_DEFINES

#include <cmath>

namespace bayes_vuln {

template<typename eT>
eT prior(const Prob<eT>& pi) {
	return arma::max(pi);
}

//sum y max x pi(x) C[x,y]
//
template<typename eT>
eT posterior(const Prob<eT>& pi, const Chan<eT>& C) {
	channel::check_prior_size(pi, C);

	if(probab::is_uniform(pi)) {
		// common case that can be optimized
		return arma::accu(arma::max(C, 0)) / (int)pi.n_cols;

	} else {
		eT s = eT(0);
		for(uint y = 0; y < C.n_cols; y++)
			s += arma::max(trans(pi) % C.col(y));
		return s;
	}
}

template<typename eT>
eT add_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return posterior(pi, C) - prior(pi);
}

template<typename eT>
eT mult_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return posterior(pi, C) / prior(pi);
}

template<typename eT>
eT min_entropy_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return real_ops<eT>::log2(mult_leakage(pi, C));
}

// sum of column maxima
//
template<typename eT>
eT mult_capacity(const Chan<eT>& C) {
	return arma::accu(arma::max(C, 0));
}

template<typename eT>
arma::ucolvec strategy(const Prob<eT>& pi, const Chan<eT>& C) {
	channel::check_prior_size(pi, C);

	arma::ucolvec strategy(pi.n_elem);
	for(uint y = 0; y < C.n_cols; y++)
		(trans(pi) % C.col(y)).max( strategy.at(y) );

	return strategy;
}

// upper bound to cap_b(n), from the recurrence formula and the bound for cap_2(n)
// see Geoffrey's POST paper
//
template<typename eT>
eT cap(uint b, uint n) {
	if(b == 1)
		return 1;

	eT cap1 = 1;														// cap_1(n) = 1
	eT cap2 = sqrt(M_PI * n / 2) + 2.0/3 + sqrt(M_PI / (2*n)) / 12;		// bound for cap_2(n)

	for(uint i = 3; i <= b; i++) {
		eT temp = cap2;
		cap2 += cap1 * n/(i-2);
		cap1 = temp;
	}
	return cap2;
}

template<typename eT>
eT mult_capacity_bound_cap(const Chan<eT>& C, uint n) {
	return cap<eT>(C.n_cols, n);
}

// Treating C's row-th row as variables, computes the row that minimizes
// the channels posterior vulnerability. The row is updated inside C and the
// optimal vulnerability is returned
//
template<typename eT = eT_def>
inline
eT optimal_row(const Prob<eT>& pi, Chan<eT>& C, uint row) {
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


} // namespace bayes_vuln
