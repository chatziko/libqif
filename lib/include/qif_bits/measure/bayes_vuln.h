namespace measure::bayes_vuln {

template<typename eT>
eT prior(const Prob<eT>& pi) {
	return arma::max(pi);
}

//sum y max x pi(x) C[x,y]
//
template<typename eT>
eT posterior(const Prob<eT>& pi, const Chan<eT>& C) {
	channel::check_prior_size(pi, C);
		
	// Use the joint formulation: V[pi, C] = sum_y max_w J_{w,y}
	//
	if(probab::is_uniform(pi))		// common case that can be optimized
		return arma::accu(arma::max(C)) / (int)pi.n_cols;
	else
		return arma::accu(arma::max(C.each_col() % pi));
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

} // namespace bayes_vuln
