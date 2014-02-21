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

#endif
