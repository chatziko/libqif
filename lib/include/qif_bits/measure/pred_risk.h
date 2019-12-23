
namespace measure::pred_risk {

// Two-block loss function, given a predicate P (boolean vector expressing whether each secret satisfies P)
//
template<typename eT = eT_def>
Mat<eT> L_pred(const arma::urowvec& P) {
	return pred_vuln::G_pred<eT>(P);		// gain function is identical to loss, the two rows simply swap meanings
}

template<typename eT>
eT prior(const arma::urowvec& P, const Prob<eT>& pi) {
	return l_risk::prior(L_pred<eT>(P), pi);
}

template<typename eT>
eT posterior(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	return l_risk::posterior(L_pred<eT>(P), pi, C);
}

template<typename eT>
eT add_leakage(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	return l_risk::add_leakage(L_pred<eT>(P), pi, C);
}

template<typename eT>
eT mult_leakage(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	return l_risk::mult_leakage(L_pred<eT>(P), pi, C);
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

	eT cap = equal(diam, eT(2))
		? infinity<eT>()
		: eT(1) / (eT(1) - diam/2);

	Prob<eT> pi(C.n_rows, arma::fill::zeros);
	pi(x1) = pi(x2) = eT(1)/2;

	return { cap, pi };
}

// Transform (pi,C) to "binary" (pibin, Cbin), modeling a system with secrets "P" and "not x",
// where "not x" acts as the average of all secrets different than x.
//
template<typename eT = eT_def>
std::pair<Prob<eT>,Chan<eT>> binary_channel(const arma::urowvec& P, const Prob<eT>& pi, const Chan<eT>& C) {
	return pred_vuln::binary_channel(P, pi, C);
}

} // namespace l_risk
