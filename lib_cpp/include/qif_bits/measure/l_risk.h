
namespace measure::l_risk {

template<typename eT>
eT prior(const Mat<eT>& L, const Prob<eT>& pi) {
	// Note: negating L works fine, although eT(-1) * L is not a "proper" loss function. Faster than using loss_to_gain
	//
	return eT(-1) * g_vuln::prior<eT>(eT(-1) * L, pi);
}

template<typename eT>
eT prior(const Metric<eT, uint>& l, const Prob<eT>& pi) {
	return prior(metric::to_distance_matrix(l, pi.n_cols), pi);
}

template<typename eT>
eT posterior(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return eT(-1) * g_vuln::posterior<eT>(eT(-1) * L, pi, C);
}

template<typename eT>
eT posterior(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return posterior(metric::to_distance_matrix(l, pi.n_cols), pi, C);
}

template<typename eT>
eT add_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return prior(L, pi) - posterior(L, pi, C);
}

template<typename eT>
eT add_leakage(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return add_leakage(metric::to_distance_matrix(l, pi.n_cols), pi, C);
}

template<typename eT>
eT mult_leakage(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return prior(L, pi) / posterior(L, pi, C);
}

template<typename eT>
eT mult_leakage(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return mult_leakage(metric::to_distance_matrix(l, pi.n_cols), pi, C);
}

template<typename eT>
arma::ucolvec strategy(const Mat<eT>& L, const Prob<eT>& pi, const Chan<eT>& C) {
	return g_vuln::strategy<eT>(eT(-1) * L, pi, C);
}

template<typename eT>
arma::ucolvec strategy(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return strategy(metric::to_distance_matrix(l, pi.n_cols), pi, C);
}

template<typename eT>
eT add_capacity(const Prob<eT>& pi, const Chan<eT>& C, bool one_spanning_g = false) {
	return g_vuln::add_capacity(pi, C, one_spanning_g);
}

// Convert a loss function to a gain function by complementing wrt its max value
//
template<typename eT>
Metric<eT,uint> loss_to_gain(
	uint n_secrets,
	uint n_guesses,
	Metric<eT, uint> loss
) {
	eT ceiling(0);

	for(uint w = 0; w < n_guesses; w++)
		for(uint x = 0; x < n_secrets; x++)
			if(eT l = loss(w,x); less_than(ceiling, l))
				ceiling = l;

	return [=](uint w, uint x) -> auto {
		return ceiling - loss(w, x);
	};
}

template<typename eT>
const Metric<eT,uint> l_zero_one = metric::discrete<eT,uint>();

template<typename eT>
Mat<eT> L_zero_one(uint n) {
	return eT(1) - g_vuln::G_id<eT>(n);
}

} // namespace measure::l_risk
