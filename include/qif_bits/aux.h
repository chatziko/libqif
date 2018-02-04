
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
// positive numbers. Absolute error can be forced by passing max_rel_diff == 0.0. We also consider inf == inf.
// see: http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
//
template<>
inline bool equal(const double& x, const double& y, const double& md, const double& mrd) {
	if(x == y) return true;

	double diff = std::abs(x - y),
		   ax = std::abs(x),
		   ay = std::abs(y),
		   largest = (ay > ax ? ay : ax);

	return largest == inf
		? false
		: diff <= (x == 0.0 || y == 0.0 || mrd == 0.0 ? md : mrd);
}

template<>
inline bool equal(const float& x, const float& y, const float& md, const float& mrd) {
	if(x == y) return true;

	float diff = std::abs(x - y),
		  ax = std::abs(x),
		  ay = std::abs(y),
		  largest = (ay > ax ? ay : ax);

	return largest == finf
		? false
		: diff <= (x == 0.0 || y == 0.0 || mrd == 0.0 ? md : mrd);
}


template<typename eT>
inline bool less_than_or_eq(const eT& x, const eT& y, const eT& md = def_max_diff<eT>(), const eT& mrd = def_max_rel_diff<eT>()) {
	return x < y || equal(x, y, md, mrd);
}

template<typename eT>
inline bool less_than(const eT& x, const eT& y, const eT& md = def_max_diff<eT>(), const eT& mrd = def_max_rel_diff<eT>()) {
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
inline eT infinity() {
	return std::numeric_limits<eT>::has_infinity
		? std::numeric_limits<eT>::infinity()
		: std::numeric_limits<eT>::max();
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
	inline static eT log2(eT) { throw std::runtime_error("not supported on this datatype"); }
};
template<>
struct real_ops<double> {
	inline static double log2(double x) { return qif::log2(x); }
};
template<>
struct real_ops<float> {
	inline static float log2(float x) { return qif::log2(x); }
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
	return arma::norm(vec, 1);
}

template<>
inline rat norm1(const Row<rat>& vec) {
	rat sum(0);
	for(auto e : vec)
		sum += abs(e);
	return sum;
}

// exp with rat support, by converting to double
template<typename eT>
inline eT exp(eT x) {
	return std::exp(x);
}

template<>
inline rat exp(rat x) {
	return rat(std::exp(x.get_d()));
}
