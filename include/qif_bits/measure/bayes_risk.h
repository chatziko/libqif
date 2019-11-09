
namespace measure::bayes_risk {

template<typename eT>
eT prior(const Prob<eT>& pi) {
	return 1 - bayes_vuln::prior(pi);
}

template<typename eT>
eT posterior(const Prob<eT>& pi, const Chan<eT>& C) {
	return 1 - bayes_vuln::posterior(pi, C);
}

template<typename eT>
eT add_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return bayes_vuln::add_leakage(pi, C);
}

template<typename eT>
eT mult_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	eT pr = prior(pi);
	if(equal(pr, eT(0))) return eT(1);
	
	eT post = posterior(pi, C);
	if(equal(pr, eT(0))) return infinity<eT>();

	return pr / post;
}

// try all 2-secret sub-uniform distributions
// (still not 100% sure whether this is true)
//
template<typename eT>
eT mult_capacity(const Chan<eT>& C, uint& res_x1, uint& res_x2) {
	auto tv = metric::total_variation<eT, Prob<eT>>();
	eT tv_max(0);

	for(uint x1 = 0; x1 < C.n_rows; x1++) {
		for(uint x2 = x1+1; x2 < C.n_rows; x2++) {
			eT tv_cur = tv(C.row(x1), C.row(x2));

			if(less_than(tv_max, tv_cur)) {
				tv_max = tv_cur;
				res_x1 = x1;
				res_x2 = x2;

				if(equal(tv_max, eT(1)))
					return infinity<eT>();
			}
		}
	}
	return eT(1) / (eT(1) - tv_max);
}

template<typename eT>
eT mult_capacity(const Chan<eT>& C) {
	uint x1, x2;
	return mult_capacity<eT>(C, x1, x2);
}

// Computes by mult_capacity by embeding C's rows in Rinf^d where d = 2^m (m is n_cols).
// The embedding phi(x) in R^d of a vector x in R^m has one coordinate for every bitstring b in {0,1}^m.
// The value of that coordinate is cdot(x, coeff) where coeff_i = (-1)^{b_i}.
//
// In Linf, the diamater is simply given by max_i (max_x x_i - min_x x_i). So for each coordinate b, we need to
// compute { phi(x)_b }_x, then compute the max - min, and finally keep the maximum of those.
//
template<typename eT>
eT mult_capacity_via_linf(const Chan<eT>& C) {
	eT diam(0), one(1);

	arma::Col<eT> coeff(C.n_cols);	// we only keep the coeff vector (b is implicit)
	coeff.fill(one);				// start with b = 00..0, coeff = (1,...,1)

	// we have a loop for every value of b. The loop finishes when b becomes 11..1, in which case
	// the update loop (which stops at the first 0) will finish all steps and i will become C.n_cols
	uint i = 0;
	while(i < C.n_cols) {
		// with a single multiplication we get cdot(coeff, row) for all rows.
		// we need to compute max-min, and update diam if the new value is bigger
		auto trans = C * coeff;
		diam = max(diam, arma::max(trans) - arma::min(trans));

		// compute the coeff for the next value of b. For b we would need to change 0s to 1s and vice versa, until we find the first 0 (which became 1).
		// In coeff, this translates changing 1s to -1s and vice versa, until we find the first -1 (which became 1).
		for(i = 0; i < C.n_cols && (coeff(i) *= -one) > eT(0); i++)
			;
	}

	return one / (one - diam/eT(2));
}

// bound mult_capacity given by the max tv distance from row
// For the lower bound to be valid, row needs to be a convex combination of the rows of C
//
template<typename eT>
std::pair<eT,eT> mult_capacity_bound1(const Chan<eT>& C, const Prob<eT>& row) {
	auto tv = metric::total_variation<eT, Prob<eT>>();
	eT tv_max(0);

	for(uint x = 0; x < C.n_rows; x++) {
		eT tv_cur = tv(C.row(x), row);

		if(less_than(tv_max, tv_cur)) {
			tv_max = tv_cur;
		}
	}

	eT one(1);
	eT d(2 * tv_max);
	eT lower(one / (one - d/2));
	eT upper(less_than(d, one) ? one / (one - d) : infinity<eT>());
	return std::pair<eT,eT>(lower, upper);
}

// bound1 with the middle row
template<typename eT>
std::pair<eT,eT> mult_capacity_bound2(const Chan<eT>& C) {
	Prob<eT> row = C.row(C.n_rows/2);
	return mult_capacity_bound1(C, row);
}

// bound1 with the average of all rows
template<typename eT>
std::pair<eT,eT> mult_capacity_bound3(const Chan<eT>& C) {
	Prob<eT> row = arma::mean(C, 0);
	return mult_capacity_bound1(C, row);
}

// bound1 with a uniform row (lower bound might not be valid)
template<typename eT>
std::pair<eT,eT> mult_capacity_bound4(const Chan<eT>& C) {
	return mult_capacity_bound1(C, probab::uniform<eT>(C.n_cols));
}

template<typename eT>
std::pair<eT,eT> mult_capacity_bound5(const Chan<eT>& C) {
	auto euclid = metric::euclidean<eT, Prob<eT>>();
	eT euclid_max(0);

	for(uint x1 = 0; x1 < C.n_rows; x1++) {
		for(uint x2 = x1+1; x2 < C.n_rows; x2++) {
			eT euclid_cur = euclid(C.row(x1), C.row(x2));

			if(less_than(euclid_max, euclid_cur)) {
				euclid_max = euclid_cur;
			}
		}
	}

	eT one(1);
	eT lbound(euclid_max / 2);
	eT ubound(euclid_max * std::sqrt(C.n_cols) / 2);
	return std::pair<eT,eT>(
		one / (one - lbound),
		less_than(ubound, one) ? one / (one - ubound) : infinity<eT>()
	);
}

template<typename eT>
arma::ucolvec strategy(const Prob<eT>& pi, const Chan<eT>& C) {
	return bayes_vuln::strategy(pi, C);
}

template<typename eT>
eT posterior_bound_via_risk_mult_cap(const Prob<eT>& pi, const Chan<eT>& C) {
	return prior<eT>(pi) / mult_capacity<eT>(C);
}

template<typename eT>
eT posterior_bound_via_vuln_mult_cap(const Prob<eT>& pi, const Chan<eT>& C) {
	return max(1 - bayes_vuln::prior<eT>(pi) * bayes_vuln::mult_capacity<eT>(C), eT(0));
}

template<typename eT = eT_def>
inline
eT optimal_row(const Prob<eT>& pi, Chan<eT>& C, uint row) {
	return 1 - bayes_vuln::optimal_row(pi, C, row);
}

} // namespace bayes_risk
