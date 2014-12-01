#ifndef _QIF_aux_h_
#define _QIF_aux_h_

#include <cmath>  /* for std::abs(double) */

const double epsilon = 1e-6;

template<typename eT>
inline bool equal(const eT& x, const eT& y) {
	// default comparison using ==
	return x == y;
}

template<>
inline bool equal(const double& x, const double& y) {
	// comparison for double, see Knuth section 4.2.2 pages 217-218
	return std::abs(x - y) <= epsilon * std::abs(x);
}

template<>
inline bool equal(const float& x, const float& y) {
	// comparison for float, see Knuth section 4.2.2 pages 217-218
	return std::abs(x - y) <= epsilon * std::abs(x);
}

template<typename eT>
inline bool less_than(const eT& x, const eT& y) {
	// in the perl library a tolerance of epsilon is allowed, but I don't
	// know why. We leave it simple for now
	return x < y;
}

template<typename eT>
inline bool less_than_or_eq(const eT& x, const eT& y) {
	return x < y || equal(x, y);
}

// like abs(x-y), but avoids negative values, cause eT might not support them!
template<typename eT>
inline eT abs_diff(const eT& x, const eT& y) {
	return x > y ? x - y : y - x;
}


// errors close to 1 are translated by log2 to errors close to 0, which are harder to test
namespace qif {
	inline float log2(float a) {
		return equal(a, 1.0f) ? 0 : std::log2(a);
	}
	inline double log2(double a) {
		return equal(a, 1.0) ? 0 : std::log2(a);
	}
	inline double log2(int a) {
		return equal(a, 1) ? 0 : std::log2(a);
	}
	inline double log2(uint a) {
		return std::log2(a);
	}
}



#endif
