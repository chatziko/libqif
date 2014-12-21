#ifndef _QIF_aux_h_
#define _QIF_aux_h_

#include <cmath>  /* for std::abs(double) */
#include <limits>
#include "types.h"


const double  inf = std::numeric_limits<double>::infinity();
const float  finf = std::numeric_limits<float>::infinity();

template<typename eT> inline eT     def_max_diff() { return eT(0); }
template<>            inline double def_max_diff() { return 1e-7; }
template<>            inline float  def_max_diff() { return 1e-7; }

template<typename eT> inline eT     def_max_rel_diff() { return eT(0); }
template<>            inline double def_max_rel_diff() { return 100 * std::numeric_limits<double>::epsilon(); }
template<>            inline float  def_max_rel_diff() { return  10 * std::numeric_limits<float >::epsilon(); }


template<typename eT>
inline bool equal(const eT& x, const eT& y, const eT& = def_max_diff<eT>(), const eT& = def_max_rel_diff<eT>()) {
	// default comparison using ==
	return x == y;
}

template<>
inline bool equal(const rat& x, const rat& y, const rat&, const rat&) {
	// for some weird reason, == doesn't work on rat
	return cmp(x, y) == 0;
}

// mixed absolute/relative error comparison. We use absolute for comparisons with 0.0, and relative for
// positive numbers. We also consider inf == inf.
// see: http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
//
template<>
inline bool equal(const double& x, const double& y, const double& max_diff, const double& max_rel_diff) {
	if(x == y) return true;

	double diff = std::abs(x - y),
		   ax = std::abs(x),
		   ay = std::abs(y),
		   largest = (ay > ax ? ay : ax);

	return	largest == inf			? false :
			x == 0.0 || y == 0.0	? (diff <= max_diff) :
									  (diff <= largest * max_rel_diff);
}

template<>
inline bool equal(const float& x, const float& y, const float& max_diff, const float& max_rel_diff) {
	if(x == y) return true;

	float diff = std::abs(x - y),
		  ax = std::abs(x),
		  ay = std::abs(y),
		  largest = (ay > ax ? ay : ax);

	return	largest == finf			? false :
			x == 0.0 || y == 0.0	? (diff <= max_diff) :
									  (diff <= largest * max_rel_diff);
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

template<typename eT>
inline eT abs(const eT& x) {
	return x < eT(0) ? -x : x;
}

// like abs(x-y), but avoids negative values, cause eT might not support them!
template<typename eT>
inline eT abs_diff(const eT& x, const eT& y) {
	return x > y ? x - y : y - x;
}

template<typename eT>
inline eT infinity() {
	return std::numeric_limits<eT>::has_infinity
		? std::numeric_limits<eT>::infinity()
		: std::numeric_limits<eT>::max();
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

	// helper template class, for operations that are defined on double/float
	// but compile and return runtime error on other types
	//
	template<typename eT>
	struct real_ops {
		inline static eT log2(eT x) { throw std::runtime_error("not supported on this datatype"); }
	};
	template<>
	struct real_ops<double> {
		inline static double log2(double x) { return qif::log2(x); }
	};
	template<>
	struct real_ops<float> {
		inline static float log2(float x) { return qif::log2(x); }
	};
}



#endif
