
namespace measure::l_risk {

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
arma::ucolvec strategy(const Metric<eT, uint>& l, const Prob<eT>& pi, const Chan<eT>& C) {
	return strategy(metric_to_mat(l, pi.n_cols), pi, C);
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

} // namespace measure::l_risk
