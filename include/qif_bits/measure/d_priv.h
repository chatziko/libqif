
namespace measure::d_priv {

template<typename eT>
bool is_private(const Chan<eT>& C, Metric<eT, uint> d) {

	return metric::is_lipschitz<eT,uint,Prob<eT>>(
		[&](uint x) -> Prob<eT> { return C.row(x); },
		d,
		metric::mult_total_variation<eT, Prob<eT>>(),
		range<uint>(0, C.n_rows)
	);
}

template<typename eT>
eT smallest_epsilon(const Chan<eT>& C, Metric<eT, uint> d) {

	return metric::lipschitz_constant<eT,uint,Prob<eT>>(
		[&](uint x) -> Prob<eT> { return C.row(x); },
		d,
		metric::mult_total_variation<eT, Prob<eT>>(),
		range<uint>(0, C.n_rows)
	);
}

template<typename eT>
eT d_vulnerability(Metric<eT, uint> d, const Prob<eT>& pi) {

	eT res(0);
	for(uint i = 0; i < pi.n_cols; i++) {
		for(uint j = i+1; j < pi.n_cols; j++) {
			// chainable elements are redundant to check
			if(d.chainable(i, j)) continue;

			eT ratio = std::abs(std::log(pi(i)) - std::log(pi(j))) / d(i, j);
			if(less_than(res, ratio))
				res = ratio;
		}
	}
	return res;
}

} // namespace measure::d_priv
