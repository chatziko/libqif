
namespace probab {

template<typename eT = eT_def>
inline
void uniform(Prob<eT>& pi) {
	// cast to uint cause rat is confused when dividing by const uint
	pi.fill( eT(1) / uint(pi.n_cols) );
}

template<typename eT = eT_def>
inline Prob<eT> uniform(uint n) {
	Prob<eT> pi(n);
	uniform(pi);
	return pi;		// separate return to allow move semantics!
}


template<typename eT = eT_def>
inline
void point(Prob<eT>& pi, uint i = 0) {
	pi.zeros();
	pi.at(i) = eT(1);
}

template<typename eT = eT_def>
inline
Prob<eT> point(uint n, uint i = 0) {
	Prob<eT> pi(n);
	point(pi, i);
	return pi;		// separate return to allow move semantics!
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
template<typename eT = eT_def>
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

template<typename eT = eT_def>
inline
Prob<eT> randu(uint n) {
	Prob<eT> pi(n);
	randu(pi);
	return pi;		// separate return to allow move semantics!
}

template<typename eT = eT_def>
inline void normalize(Prob<eT>& pi) {
	pi /= arma::accu(pi);
}

template<typename eT = eT_def>
Prob<eT> normalize(const Prob<eT>& pi) {
	Prob<eT> res = pi;
	normalize(res);
	return res;
}


// sample from pi. Allow pi to be anything iterable
//
template<typename eT = eT_def, typename T = Prob<eT>>
inline
uint sample(const T& pi) {
	eT p = rng::randu<eT>();

	eT accu(0);
	uint cur = 0;
	for(auto pi_it = pi.begin(); pi_it != pi.end(); ++pi_it) {
		if(!less_than(accu += *pi_it, p))
			return cur;
		cur++;
	}

	return pi.n_cols - 1;
}

// sample multiple samples efficiently (with a single iteration over pi)
//
template<typename eT = eT_def, typename T = Prob<eT>>
inline
Row<uint> sample(const T& pi, uint n) {
	// we need n numbers uniformly sampled in [0,1]. Sort them and keep the indexes of the sorted list in 'order'
	Row<eT> ps(n);
	ps.randu();
	arma::uvec order = arma::sort_index(ps);

	eT accu(0);
	int cur = -1;
	auto pi_it = pi.begin();	// single iteration over pi for all samples
	Row<uint> res(n);

	// sample n elements. 
	for(uint i = 0; i < n; i++) {
		// we need to visit elements in sorted order of p
		eT p = ps(order(i));

		while(less_than(accu, p) && pi_it != pi.end()) {
			accu += *pi_it;
			pi_it++;
			cur++;
		}

		// place the result in the same position in res as p was in ps.
		res(order(i)) = cur;
	}

	return res;
}

template<typename eT = eT_def>
inline
bool is_uniform(const Prob<eT>& pi, const eT& mrd = def_mrd<eT>) {
	eT v = eT(1) / (int)pi.n_cols;
	for(uint j = 0; j < pi.n_cols; j++)
		if(!qif::equal(v, pi(j), eT(0), mrd))
			return false;

	return true;
}

template<typename eT = eT_def>
inline
bool is_proper(const Prob<eT>& pi, const eT& mrd = def_mrd<eT>) {
	eT sum(0);
	for(uint j = 0; j < pi.n_cols; j++) {
		// elements should be non-negative
		const eT& elem = pi.at(j);
		if(less_than(elem, eT(0)))
			return false;

		sum += elem;
	}

	// sum should be 1
	if(!qif::equal(sum, eT(1), def_md<eT>, mrd))
		return false;

	return true;
}

template<typename eT = eT_def>
inline
void assert_proper(const Prob<eT>& pi) {
	if(!is_proper(pi))
		throw std::runtime_error("not a proper dist");
}

template<typename eT = eT_def>
inline bool equal(const Prob<eT>& A, const Prob<eT>& B, const eT& md = def_md<eT>, const eT& mrd = def_mrd<eT>) {
	if(A.n_cols != B.n_cols)
		return false;

	for(uint i = 0; i < A.n_cols; i++)
		if(!qif::equal(A.at(i), B.at(i), md, mrd))
			return false;

	return true;
}

// Takes a probability pi on a grid of given width.
// Returns a grid representation of that probability.
//
template<typename eT = eT_def>
Mat<eT> to_grid(const Prob<eT>& pi, uint width) {
	return arma::trans(arma::reshape(pi, pi.n_elem / width, width));
}

// Vectorizes a probability on a grid.
// Element (i,j) is mapped to i+j*width.
//
template<typename eT = eT_def>
Prob<eT> from_grid(const Mat<eT>& grid) {
	return arma::vectorise(grid, 1);
}

} // namespace probab
