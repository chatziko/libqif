
namespace bayes_risk {

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
eT mult_capacity(const Chan<eT>& C) {
	Chan<eT> Csmall(2, C.n_cols);

	eT max(0);
	for(uint x1 = 0; x1 < C.n_rows; x1++) {
		Csmall.row(0) = C.row(x1);

		for(uint x2 = x1+1; x2 < C.n_rows; x2++) {
			Csmall.row(1) = C.row(x2);

			eT sum_of_maxima = arma::accu(arma::max(Csmall, 0));
			if(equal(sum_of_maxima, eT(2))) return infinity<eT>();

			eT l = eT(1) / (eT(2) - sum_of_maxima);
			if(less_than(max, l))
				max = l;
		}
	}
	return max;
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


} // namespace bayes_risk
