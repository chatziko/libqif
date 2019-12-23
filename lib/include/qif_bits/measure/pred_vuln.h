
namespace measure::pred_vuln {

// Two-block gain function, given a predicate P (boolean vector expressing whether each secret satisfies P)
//
template<typename eT = eT_def>
Mat<eT> G_pred(const arma::urowvec& P) {
	Mat<eT> G(2, P.n_elem, arma::fill::zeros);
	auto PP = arma::conv_to<Row<eT>>::from(P);
	G.row(0) = PP;
	G.row(1) = eT(1) - PP;
	return G;
}

template<typename eT>
eT prior(const arma::urowvec& P, const Prob<eT>& pi) {
	return g_vuln::prior(G_pred<eT>(P), pi);
}

template<typename eT>
eT posterior(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	return g_vuln::posterior(G_pred<eT>(P), pi, C);
}

template<typename eT>
eT add_leakage(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	return g_vuln::add_leakage(G_pred<eT>(P), pi, C);
}

template<typename eT>
eT mult_leakage(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	return g_vuln::mult_leakage(G_pred<eT>(P), pi, C);
}

// Mult capacity is given by the max l1 distance between the rows P and not P
//
template<typename eT>
std::pair<eT,Prob<eT>> mult_capacity(const arma::urowvec& P, const Chan<eT>& C, std::string method = "direct") {

	auto [diam, x1, x2] = metric::optimize::l1_max_distance(
		(Chan<eT>)C.rows(arma::find(P == 1)),
		(Chan<eT>)C.rows(arma::find(P == 0)),
		method
	);

	Prob<eT> pi(C.n_rows, arma::fill::zeros);
	pi(x1) = pi(x2) = eT(1)/2;

	return { eT(1) + diam/2, pi };
}

// Transform (pi,C) to "binary" (pibin, Cbin), modeling a system with secrets "P" and "not x",
// where "not x" acts as the average of all secrets different than x.
//
template<typename eT = eT_def>
std::pair<Prob<eT>,Chan<eT>> binary_channel(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	// we can simply use g_to_bayes
	auto [rho, R, a, b] = g_vuln::g_to_bayes(G_pred<eT>(P), pi);

	// for two block gain functions a,b should always be 1 and 0 respectively.
	(void)a; (void)b;
	assert(equal(a, eT(1)));
	assert(equal(b, eT(0)));

	return { rho, R*C };
}

} // namespace g_vuln
