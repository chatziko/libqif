
namespace l {

template<typename eT>
eT uncertainty(const Mat<eT>& L, const Prob<eT>& pi) {
	return -g::vulnerability<eT>(-L, pi);
}

template<typename eT>
eT post_uncertainty(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return -g::post_vulnerability<eT>(-L, pi, C);
}

template<typename eT>
eT add_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return post_uncertainty(L, pi, C) - uncertainty(L, pi);
}

template<typename eT>
eT mult_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return uncertainty(L, pi) / post_uncertainty(L, pi, C);
}

template<typename eT>
eT mulg_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return real_ops<eT>::log2(mult_leakage(L, pi, C));
}

template<typename eT>
arma::ucolvec strategy(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return g::strategy<eT>(-L, pi, C);
}

/*
OLD CODE, we now use the g:: methods with -L

// same as cond_vulnerability but with min. G is assumed to be a loss function
//
template<typename eT>
eT GLeakage<eT>::bayes_risk(const Prob<eT>& pi) {
	check_prior(pi);

	eT s = eT(0);
	for(uint y = 0; y < C.n_cols; y++)
		s += arma::min(G * (trans(pi) % C.col(y)));
	return s;
}

template<typename eT>
arma::ucolvec GLeakage<eT>::bayes_strategy(const Prob<eT>& pi) const {
	check_prior(pi);

	arma::ucolvec strategy(pi.n_elem);
	for(uint y = 0; y < C.n_cols; y++)
		(G * (trans(pi) % C.col(y))).min( strategy.at(y) );

	return strategy;
}
*/

} // namespace l
