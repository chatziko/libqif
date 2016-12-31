
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

	arma::ucolvec strategy(C.n_cols);
	for(uint y = 0; y < C.n_cols; y++)
		(G * (trans(pi) % C.col(y))).max( strategy.at(y) );

	return strategy;
}

template<typename eT>
eT add_capacity(const Prob<eT>& pi, const Chan<eT>& C, bool one_spanning = true) {
	channel::check_prior_size(pi, C);

	// we essentially compute the Kantorovich distance between the hypers
	// [pi] and [pi, C], where the underlying metric between the inners is pdist
	// (CSF'14, Theorem 17). Since [pi] is a point hyper, computing the transportation
	// problem is easy, since the full 1 mass of pi needs to be transferred to all
	// the posteriors. Each posterior sigma^y needs to receive mass exactly delta(y),
	// for a total cost of sum_y delta(y) * pdist(pi, sigma^y)
	//
	auto pdist = one_spanning
		? metric::total_variation<eT, Prob<eT>>()
		: metric::bounded_entropy_distance<eT, Prob<eT>>();
	Prob<eT> outer = pi * C;
	eT res(0);

	for(uint y = 0; y < outer.n_elem; y++)
		if(!equal(outer(y), eT(0)))		// ignore 0 prob outputs
			res += outer(y) * pdist(pi, channel::posterior(C, pi, y));

	return res;
}

// TODO: remove _kant method
template<typename eT>
eT add_capacity_kant(const Prob<eT>& pi, const Chan<eT>& C, bool one_spanning = true) {
	channel::check_prior_size(pi, C);

	// we need construct two hypers: [pi,C] and [pi]
	// we need them to have the same support, so we first construct
	// [pi, C] and extend the inners matrix (its support) with pi.
	// outer is then extended with a 0 and pointhyper is a dirac on that extra column.
	// (note: pi might exist in inners already but it's not really a problem)
	//
	Mat<eT> inners;
	Prob<eT> outer = channel::hyper(C, pi, inners);

	inners.insert_cols(inners.n_cols, pi.t());

	outer.resize(inners.n_cols);
	outer(outer.n_elem-1) = 0;

	Prob<eT> pointhyper = probab::dirac<eT>(inners.n_cols, inners.n_cols-1);

	// kantorovich needs a metric on the elements of the distribution (just uints)
	// This is pdist_inners, which applies pdist on the columns of inners
	//
	auto pdist = one_spanning
		? metric::total_variation<eT, Prob<eT>>()
		: metric::bounded_entropy_distance<eT, Prob<eT>>();
	auto pdist_inners = [&pdist, &inners](const uint& a, const uint& b) -> eT {
		return pdist(trans(inners.col(a)), trans(inners.col(b)));
	};
	auto kant = metric::kantorovich<eT, Prob<eT>>(pdist_inners);

	return kant(pointhyper, outer);		// CSF'14, Theorem 17
}

} // namespace g
