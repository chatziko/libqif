
// Mechanism design for Bayes risk

namespace mechanism::bayes_risk {

// Returns the mechanism having the smallest E[loss] given the R(pi,C) >= min_risk constraint.
//
template<typename eT>
Chan<eT> min_loss_given_min_risk(
	const Prob<eT>& pi,
	uint n_cols,
	eT min_risk,
	Metric<eT, uint> loss,
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	return bayes_vuln::min_loss_given_max_vuln(pi, n_cols, 1-min_risk, loss, hard_max_loss);
}

// Returns the mechanism having the largest Bayes risk given the E[loss] <= max_loss constraint.
//
template<typename eT>
Chan<eT> max_risk_given_max_loss(
	const Prob<eT>& pi,
	uint n_cols,
	eT max_loss,
	Metric<eT, uint> loss,
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	return bayes_vuln::min_vuln_given_max_loss(pi, n_cols, max_loss, loss, hard_max_loss);
}

} // namespace mechanism::bayes_risk
