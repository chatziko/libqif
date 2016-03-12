
namespace probab {

template<typename eT>
inline
Prob<eT>& uniform(Prob<eT>& pi) {
	// cast to uint cause rat is confused when dividing by const uint
	pi.fill( eT(1) / uint(pi.n_cols) );
	return pi;
}

template<typename eT>
inline Prob<eT> uniform(uint n) {
	Prob<eT> pi(n);
	return uniform(pi);
}


template<typename eT>
inline
Prob<eT>& dirac(Prob<eT>& pi, uint i = 0) {
	pi.zeros();
	pi.at(i) = eT(1);
	return pi;
}

template<typename eT>
inline
Prob<eT> dirac(uint n, uint i = 0) {
	Prob<eT> pi(n);
	return dirac(pi, i);
}


// Generate a distribution of n elements uniformly (among all the elements of the n-1 simplex)
//
// Algorithm: Generate n-1 numbers uniformly in [0,1], add 0 and 1 to the list,
// sort, and take the difference between consecutive elements.
// This algorithm has advantages over the simple "normalize n uniform elements" tehcnique:
//
// 1. The normalizing technique is not uniform! See the url below.
//
// 2. The algorithm involves no divisions, and the resulting sum is much closer to exactly 1.0
//    Using this algorithm with float the Kantorovich tests over random dists always pass, while with
//    the normalizing algorithm there are instabilities.
//
// 3. We avoid to sum all elements, which creates a big denominator under rats
//
// Note: if the n logn complexity is a problem, there's a linear algorithm involving logs in the following url:
// http://stats.stackexchange.com/questions/14059/generate-uniformly-distributed-weights-that-sum-to-unity
//
template<typename eT>
inline
Prob<eT>& randu(Prob<eT>& pi) {
	pi.randu();
	pi(pi.n_cols-1) = eT(1);		// add 1 to the list. We don't really need to add 0
	pi = arma::sort(pi);

	for(uint i = pi.n_cols-1; i > 0; i--)
		pi(i) -= pi(i-1);

	// naif normalize algorithm
	//pi.randu();
	//pi /= arma::accu(pi);

	return pi;
}

template<typename eT>
inline
Prob<eT> randu(uint n) {
	Prob<eT> pi(n);
	return randu(pi);
}


template<typename eT>
inline
bool is_proper(const Prob<eT>& pi, const eT& mrd = def_max_rel_diff<eT>()) {
	eT sum = 0;
	for(uint j = 0; j < pi.n_cols; j++) {
		// elements should be non-negative
		const eT& elem = pi.at(j);
		if(less_than(elem, eT(0)))
			return false;

		sum += elem;
	}

	// sum should be 1
	if(!equal(sum, eT(1), def_max_diff<eT>(), mrd))
		return false;

	return true;
}

template<typename eT>
inline
void check_proper(const Prob<eT>& pi) {
	if(!is_proper(pi))
		throw std::runtime_error("not a proper dist");
}


// Computes in place the (euclidean) projection of vector x \in R^n onto the
// n-simplex. That is, finds the point of the simplex that is closer  to x.
// Uses the algorithm of:
//    http://www.springerlink.com/content/q1636371674m36p1/
//
template<typename eT>
inline
Prob<eT>& project_to_simplex(Prob<eT>& x) {
	uint n = x.n_cols;
	Row<char> done(n);
	done.fill(0);

	while(true) {
		eT t = (arma::accu(x)-1)/n;
		uint n1 = 0;

		for(uint i = 0; i < x.n_cols; i++) {
			if(done(i)) continue;

			x(i) -= t;
			if(eT(0) > x(i)) {
				x(i) = eT(0);
				done(i) = 1;
				n1++;
			}
		}

		if(n1 == 0) break;		// no negative elements
		n -= n1;
	}

	return x;
}

} // namespace probab
