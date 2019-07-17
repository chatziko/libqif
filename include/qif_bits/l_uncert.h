
namespace l_uncert {

template<typename eT>
eT prior(const Mat<eT>& L, const Prob<eT>& pi) {
	return -g_vuln::prior<eT>(-L, pi);
}

template<typename eT>
eT posterior(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return -g_vuln::posterior<eT>(-L, pi, C);
}

template<typename eT>
eT add_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return prior(L, pi) - posterior(L, pi, C);
}

template<typename eT>
eT mult_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return prior(L, pi) / posterior(L, pi, C);
}

template<typename eT>
eT mulg_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return real_ops<eT>::log2(mult_leakage(L, pi, C));
}

template<typename eT>
arma::ucolvec strategy(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return g_vuln::strategy<eT>(-L, pi, C);
}


// same methods, with Metric loss instead of Mat L.
// Guesses are assumed to be the inputs
//
template<typename eT>
inline
Mat<eT> metric_to_mat(Metric<eT, uint> l, uint n) {
	Mat<eT> L(n, n);
	for(uint i = 0; i < n; i++)
		for(uint j = 0; j < n; j++)
			L(i, j) = l(i, j);
	return L;
}

template<typename eT>
eT prior(const Metric<eT, uint>& l, const Prob<eT>& pi) {
	return prior(metric_to_mat(l, pi.n_cols), pi);
}

template<typename eT>
eT posterior(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return posterior(metric_to_mat(l, pi.n_cols), pi, C);
}

template<typename eT>
eT add_leakage(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return add_leakage(metric_to_mat(l, pi.n_cols), pi, C);
}

template<typename eT>
eT mult_leakage(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return mult_leakage(metric_to_mat(l, pi.n_cols), pi, C);
}

template<typename eT>
eT mulg_leakage(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return mulg_leakage(metric_to_mat(l, pi.n_cols), pi, C);
}

template<typename eT>
arma::ucolvec strategy(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return strategy(metric_to_mat(l, pi.n_cols), pi, C);
}

template<typename eT>
eT add_capacity(const Prob<eT>& pi, const Chan<eT>& C, bool one_spanning_g = false) {
	return g_vuln::add_capacity(pi, C, one_spanning_g);
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
