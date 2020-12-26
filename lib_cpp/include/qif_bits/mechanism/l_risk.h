
// Mechanism design for l-risk

namespace mechanism::l_risk {

// Returns the mechanism having the smallest E[u_loss] given the Rl(pi,C) >= min_risk constraint.
//
template<typename eT>
Chan<eT> min_loss_given_min_risk(
	const Prob<eT>& pi,
	uint n_cols,
	uint n_guesses,
	eT min_risk,
	Metric<eT, uint> adv_loss,			// adversary loss
	Metric<eT, uint> loss,				// utility loss
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	auto adv_gain = measure::l_risk::loss_to_gain(pi.n_elem, n_guesses, adv_loss);
	eT max_vuln = adv_gain(0,0) + adv_loss(0,0) - min_risk;		// gain(0,0) + loss(0,0) is the ceiling

	return g_vuln::min_loss_given_max_vuln(pi, n_cols, pi.n_elem, max_vuln, adv_gain, loss, hard_max_loss);
}

// Returns the mechanism having the largest Rl[pi,C] given the E[u_loss] <= max_loss constraint.
//
template<typename eT>
Chan<eT> max_risk_given_max_loss(
	const Prob<eT>& pi,
	uint n_cols,
	uint n_guesses,
	eT max_loss,
	Metric<eT, uint> adv_loss,			// adversary loss
	Metric<eT, uint> loss,				// utility loss
	eT hard_max_loss = infinity<eT>()	// C[x,y] is forced to 0 if loss(x,y) > hard_max_loss
) {
	auto adv_gain = measure::l_risk::loss_to_gain(pi.n_elem, n_guesses, adv_loss);

	return g_vuln::min_vuln_given_max_loss(pi, n_cols, pi.n_elem, max_loss, adv_gain, loss, hard_max_loss);
}

} // namespace mechanism::l_risk
