
namespace measure::g_vuln {

//		void * compare_over_prior(chan& other_channel);
//		void * compare_over_gain(chan& other_channel,Prob<eT>& prior);

template<typename eT>
inline
void check_g_size(const Mat<eT>& G, const Prob<eT>& pi) {
	if(G.n_cols != pi.n_cols)
		throw std::runtime_error("invalid prior size");
}

template<typename eT>
inline
void check_g_size(const Mat<eT>& G1, const Mat<eT>& G2) {
	if(G1.n_cols != G2.n_cols)
		throw std::runtime_error("invalid G size");
}

// max_w sum_x pi[x] G[w, x]
//
template<typename eT>
eT prior(const Mat<eT>& G, const Prob<eT>& pi) {
	check_g_size(G, pi);

	return arma::max(G * trans(pi));
}

// metric-loss version, guesses are assumed to be the inputs
template<typename eT>
eT prior(const Metric<eT, uint>& g, const Prob<eT>& pi) {
	return prior(metric::to_distance_matrix(g, pi.n_cols), pi);
}

// sum_y max_w sum_x pi[x] C[x, y] G[w, x]
//
template<typename eT>
eT posterior(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	check_g_size(G, pi);
	channel::check_prior_size(pi, C);

	eT s = eT(0);
	for(uint y = 0; y < C.n_cols; y++)
		s += arma::max(G * (trans(pi) % C.col(y)));
	return s;
}

template<typename eT>
eT posterior(const Metric<eT, uint>& g, const Prob<eT>& pi, const Chan<eT>& C) {
	return posterior(metric::to_distance_matrix(g, pi.n_cols), pi, C);
}

template<typename eT>
eT add_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return posterior(G, pi, C) - prior(G, pi);
}

template<typename eT>
eT add_leakage(const Metric<eT, uint>& g, const Prob<eT>& pi, const Chan<eT>& C) {
	return add_leakage(metric::to_distance_matrix(g, pi.n_cols), pi, C);
}

template<typename eT>
eT mult_leakage(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return posterior(G, pi, C) / prior(G, pi);
}

template<typename eT>
eT mult_leakage(const Metric<eT, uint>& g, const Prob<eT>& pi, const Chan<eT>& C) {
	return mult_leakage(metric::to_distance_matrix(g, pi.n_cols), pi, C);
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
arma::ucolvec strategy(const Metric<eT, uint>& g, const Prob<eT>& pi, const Chan<eT>& C) {
	return strategy(metric::to_distance_matrix(g, pi.n_cols), pi, C);
}


// additive capacity for fixed pi and g ranging over 1-spanning Vg's (larger class, default) or
// 1-spanning g's (if one_spanning_g == true)
//
template<typename eT>
eT add_capacity(const Prob<eT>& pi, const Chan<eT>& C, bool one_spanning_g = false) {
	channel::check_prior_size(pi, C);

	eT res(0);

	if(one_spanning_g) {
		// 1-spanning g's (smaller class). We compute the Kantorovich distance between the
		// hypers [pi] and [pi, C], where the underlying metric between the inners is tv
		// (CSF'14, Theorem 17). Since [pi] is a point hyper, computing the transportation
		// problem is easy, since the full 1 mass of pi needs to be transferred to all
		// the posteriors. Each posterior sigma^y needs to receive mass exactly delta(y),
		// for a total cost of sum_y delta(y) * tv(pi, sigma^y)
		//
		auto tv = metric::total_variation<eT, Prob<eT>>();
		// metric::convex_separation_quasi<eT, Prob<eT>>();
		Prob<eT> outer = pi * C;

		for(uint y = 0; y < outer.n_elem; y++)
			if(!equal(outer(y), eT(0)))		// ignore 0 prob outputs
				res += outer(y) * tv(pi, channel::posterior(C, pi, y));

		
	} else {
		// For the larger class of 1-spanning Vg's, the capacity only depends on the support of pi and is
		// equal to 1 - the sum of column minima (including only rows in the support of pi).
		// The same result can also be obtained via the Kantorovich above, replacing tv with convex_separation_quasi.
		// 
		res = 1;
		for(uint y = 0; y < C.n_cols; y++) {
			eT min(1);
			for(uint x = 0; x < C.n_rows; x++)
				if(!equal(pi(x), eT(0)) && C(x,y) < min)
					min = C(x,y);
			res -= min;
		}
	}

	return res;
}


// mult leakage bound (even for negative g) coming from the miracle theorem, adjusted so that the minimum gain is exactly 0
template<typename eT>
eT mult_leakage_bound1(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	// NOTE: wrt the notation of the book  (Thm 4.7), arma::min(G, 0) is _minus_ the kappa vector we need to add to G to make it non-negative.
	// So the z here and the lambda of the book are related by lambda = z / (z-1)
	//
	eT z = arma::cdot(pi, arma::min(G, 0)) / prior<eT>(G, pi); 
	// std::cout << "z:" << z << "\n";
	return bayes_vuln::mult_capacity(C) * (1-z) + z;
}

template<typename eT>
eT posterior_bound1(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return prior<eT>(G, pi) * mult_leakage_bound1(G, pi, C);
}

// mult leakage bound (even for negative g) coming from the additive theorem
template<typename eT>
eT mult_leakage_bound2(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return (G.max() - G.min()) * (eT(1) - channel::sum_column_min<eT>(C)) /  prior<eT>(G, pi) + eT(1);
}

// add leakage bound (even for g not bounded by 1) coming from the additive miracle theorem, adjusted so that the max gain is exactly 1
template<typename eT>
eT add_leakage_bound1(const Mat<eT>& G, const Chan<eT>& C) {
	return G.max() * (eT(1) - channel::sum_column_min<eT>(C));
}

template<typename eT>
eT posterior_bound2(const Mat<eT>& G, const Prob<eT>& pi, const Chan<eT>& C) {
	return prior<eT>(G, pi) + add_leakage_bound1(G, C);
}



////////////////// Gain function manipulation //////////////////////////


template<typename eT>
const Metric<eT,uint> g_id = [](uint w, uint x) -> eT {
	return x == w ? eT(1) : eT(0);
};

template<typename eT>
Mat<eT> G_id(uint n) {
	return arma::eye<arma::Mat<eT>>(n, n);
}

// Adding g1+g2 produces a g such that Vg = Vg1 + Vg2
//
template<typename eT>
Mat<eT> g_add(const Mat<eT>& G1, const Mat<eT>& G2) {
	check_g_size(G1, G2);

	Mat<eT> G(G1.n_rows * G2.n_rows, G1.n_cols);
	G.fill(eT(0));

	for(uint i = 0; i < G1.n_rows; i++)
		for(uint j = 0; j < G2.n_rows; j++)
			G.row(i * G2.n_rows + j) = G1.row(i) + G2.row(j);

	return G;
}

// g such that Vg(pi) = Vg'(pi,C)
//
// 20191216: Is this really independend of pi? It seems to work only for uniform!
//
template<typename eT>
Mat<eT> g_from_posterior(const Mat<eT>& G, const Chan<eT>& C) {
	if(G.n_cols != C.n_rows)
		throw std::runtime_error("invalid G size");

	Mat<eT> Gres(1, G.n_cols);
	Gres.fill(eT(0));

	for(uint y = 0; y < C.n_cols; y++) {
		Mat<eT> Gtemp = G;
		Gtemp.each_row() %= arma::trans(C.col(y));
		Gres = g_add(Gres, Gtemp);
	}

	return Gres;
}

} // namespace g_vuln
