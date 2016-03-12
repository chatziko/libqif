
namespace g {

//		void * compare_over_prior(chan& other_channel);
//		void * compare_over_gain(chan& other_channel,Prob<eT>& prior);

template<typename eT>
inline
void check_g_size(const Mat<eT>& G, const Prob<eT>& pi) {
	if(G.n_cols != pi.n_cols)
		throw std::runtime_error("invalid prior size");
}

// max_w sum_x pi[x] G[w, x]
//
template<typename eT>
eT vulnerability(const Mat<eT>& G, const Prob<eT>& pi) {
	check_g_size(G, pi);

	return arma::max(G * trans(pi));
}

// sum_y max_w sum_x pi[x] C[x, y] G[w, x]
//
template<typename eT>
eT post_vulnerability(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	check_g_size(G, pi);
	channel::check_prior_size(pi, C);

	eT s = eT(0);
	for(uint y = 0; y < C.n_cols; y++)
		s += arma::max(G * (trans(pi) % C.col(y)));
	return s;
}

template<typename eT>
eT add_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return vulnerability(G, pi) - post_vulnerability(G, pi, C);
}

template<typename eT>
eT mult_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return post_vulnerability(G, pi, C) / vulnerability(G, pi);
}

template<typename eT>
eT mulg_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return real_ops<eT>::log2(mult_leakage(G, pi, C));
}

template<typename eT>
arma::ucolvec strategy(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	check_g_size(G, pi);
	channel::check_prior_size(pi, C);

	arma::ucolvec strategy(pi.n_elem);
	for(uint y = 0; y < C.n_cols; y++)
		(G * (trans(pi) % C.col(y))).max( strategy.at(y) );

	return strategy;
}

} // namespace g
