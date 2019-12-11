
// for debugging define values
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "="  VALUE(var)
//#pragma message(VAR_NAME_VALUE(ARMA_USE_CXX11))

// armadillo 8.0 changed set_stream2 to set_cerr_stream
#ifdef ARMA_COUT_STREAM
#define ARMA_SET_CERR(stream) arma::set_cerr_stream(stream)
#else
#define ARMA_SET_CERR(stream) arma::set_stream_err2(stream)
#endif

template<typename eT> const eT     def_md			= eT(0);
template<>            const double def_md<double>	= 1e-7;
template<>            const float  def_md<float>	= 1e-7;

template<typename eT> const eT     def_mrd			= eT(0);
template<>            const double def_mrd<double>	= 100 * std::numeric_limits<double>::epsilon();
template<>            const float  def_mrd<float>	=  10 * std::numeric_limits<float >::epsilon();


template<typename eT>
inline bool equal(const eT& x, const eT& y, const eT& md = def_md<eT>, const eT& mrd = def_mrd<eT>) {
	if constexpr (std::is_same<eT, double>::value || std::is_same<eT, float>::value) {
		// mixed absolute/relative error comparison. We use absolute for comparisons with 0.0, and relative for
		// positive numbers. Absolute error can be forced by passing max_rel_diff == 0.0. We also consider inf == inf.
		// see: http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
		//
		if(x == y) return true;

		eT diff = std::abs(x - y),
			ax = std::abs(x),
			ay = std::abs(y),
			largest = (ay > ax ? ay : ax);

		return largest == std::numeric_limits<eT>::infinity()
			? false
			: diff <= (x == 0.0 || y == 0.0 || mrd == 0.0 ? md : mrd);

	} else {
		// default comparison using ==
		return x == y;
	}
}


template<typename eT>
inline bool less_than_or_eq(const eT& x, const eT& y, const eT& md = def_md<eT>, const eT& mrd = def_mrd<eT>) {
	return x < y || equal(x, y, md, mrd);
}

template<typename eT>
inline bool less_than(const eT& x, const eT& y, const eT& md = def_md<eT>, const eT& mrd = def_mrd<eT>) {
	// strictly less_than
	return !less_than_or_eq(y, x, md, mrd);
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
inline eT max(const eT& a, const eT& b) {
	return less_than<eT>(a, b) ? b : a;
}

template<typename eT>
inline eT infinity() {
	if constexpr (std::is_same<eT, rat>::value) {
		return rat(std::numeric_limits<long>::max(), 1);			// is this sufficiently large?

	} else if constexpr (std::numeric_limits<eT>::has_infinity) {
		return std::numeric_limits<eT>::infinity();

	} else {
		return std::numeric_limits<eT>::max();
	}
}


// errors close to 1 are translated by log2 to errors close to 0, which are harder to test
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
#if defined(ARMA_64BIT_WORD)
// arma::uword is different than uint
inline double log2(arma::uword a) {
	return std::log2(a);
}
#endif

// helper template class, for operations that are defined on double/float
// but compile and return runtime error on other types
//
template<typename eT>
struct real_ops {
	inline static eT log2(eT x) {
		if constexpr (std::is_same<eT, float>::value || std::is_same<eT, double>::value) {
			return qif::log2(x);
		} else {
			throw std::runtime_error("not supported on this datatype");
		}
	}
};

// precise sum using Kahan algorithm (https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
//
template<typename eT>
class LargeSum {
	private:
	eT c = 0;
	eT val = 0;

	public:
	inline eT add(eT n) {
		eT y = n - c; 		// So far, so good: c is zero.
		eT t = val + y;		// Alas, val is big, y small, so low-order digits of y are lost.
		c = (t - val) - y;	// (t - val) recovers the high-order part of y; subtracting y recovers -(low part of y)
		return val = t;		// Algebraically, c should always be zero. Beware eagerly optimising compilers!
	}

	inline eT value() {
		return val;
	}
};

// iterative mean (http://www.heikohoffmann.de/htmlthesis/node134.html)
//
template<typename eT>
class LargeAvg {
	private:
	LargeSum<eT> sum;
	uint samples = 1;

	public:
	inline eT add(eT n) {
		eT diff = (n - sum.value()) / samples++;
		return sum.add(diff);
	}

	inline eT value() {
		return sum.value();
	}
};

// L-1 norm, with faster armadillo implementation for double/float
template<typename T>
inline typename T::elem_type norm1(const T& vec) {
	if constexpr (std::is_same<T, rprob>::value) {
		return arma::accu(arma::abs(vec));
	} else {
		return arma::norm(vec, 1);
	}
}

// exp with rat support, by converting to double
template<typename eT>
inline eT exp(eT x) {
	if constexpr (std::is_same<eT, rat>::value) {
		return rat(std::exp((double)x));
	} else {
		return std::exp(x);
	}
}

// pi with rat support, by converting to double
template<typename eT>
inline eT pi() {
	return eT(std::atan(1)*4);
}

// convert float/double/rat to double
template<typename eT>
inline double to_double(eT x) {
	return (double)x;
}