
namespace bayes {

template<typename eT>
eT vulnerability(const Prob<eT>& pi) {
	return arma::max(pi);
}

//sum y max x pi(x) C[x,y]
//
template<typename eT>
eT post_vulnerability(const Prob<eT>& pi, const Chan<eT>& C) {
	check_prior_size(pi, C);

	eT s = eT(0);
	for(uint y = 0; y < C.n_cols; y++)
		s += arma::max(trans(pi) % C.col(y));
	return s;
}

template<typename eT>
eT add_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return vulnerability(pi) - post_vulnerability(pi, C);
}

template<typename eT>
eT mult_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return post_vulnerability(pi, C) / vulnerability(pi);
}

template<typename eT>
eT mlog_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
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
	check_prior_size(pi, C);

	arma::ucolvec strategy(pi.n_elem);
	for(uint y = 0; y < C.n_cols; y++)
		(trans(pi) % C.col(y)).max( strategy.at(y) );

	return strategy;
}


} // namespace bayes
