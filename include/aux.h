#ifndef _QIF_aux_h_
#define _QIF_aux_h_

#include <cmath>  /* for std::abs(double) */

const double epsilon = 1e-6;

inline bool equal(double x, double y) {
	// see Knuth section 4.2.2 pages 217-218
	return std::abs(x - y) <= epsilon * std::abs(x);
}

inline bool less_than(double x, double y) {
	// in the perl library a tolerance of epsilon is allowed, but I don't
	// know why. We leave it simple for now
	return x < y;
}

inline bool less_than_or_eq(double x, double y) {
	return x < y || equal(x, y);
}

#endif
