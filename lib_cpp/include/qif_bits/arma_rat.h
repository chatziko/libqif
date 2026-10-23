
namespace arma {

using qif::rat;

// register rational<1> as a real type.
template<>
struct arma_real_only<rat> {
	typedef rat result;
};

// use the generic (non-BLAS) implementations of dot/norm for rat. Armadillo selects them via
// arma_fp16_real_only (fp16 is another real type not supported by BLAS), so we register rat there.
//
template<>
struct arma_fp16_real_only<rat> {
	typedef rat result;
};

// for abs
//
template<>
arma_inline rat eop_aux::arma_abs<rat>(const rat x) {
	return mppp::abs(x);
}

// eop_aux::sqrt is only defined for real/cx/integral types. Mat's eOp constructors instantiate
// eop_sqrt (to optimize pow(X, 0.5)) for all types, even though it is only used for real/cx ones.
// So we provide sqrt for rat (via double, like armadillo does for integral types).
//
template<> template<>
arma_inline rat eop_core<eop_sqrt>::process<rat>(const rat val, const rat) {
	return rat(std::sqrt(double(val)));
}

// most_neg/most_pos are the initial values of max/min (and other reductions). They use std::numeric_limits, which is
// not specialized for rat, so they would be 0: max() of a negative matrix (or min() of a positive one) would return 0.
// We use -/+ qif's rat infinity instead.
//
namespace priv {
template<> inline rat most_neg<rat>() { return -qif::infinity<rat>(); }
template<> inline rat most_pos<rat>() { return  qif::infinity<rat>(); }
}

// ------------------------------------------------

// register as supported
template<>
struct is_supported_elem_type<rat> {
	static const bool value = true;
};

// armadillo's memory::acquire/release uses malloc/free for speed, which ignores the constructor
// We override it to use C++'s new/delete.
//
template<>
inline
arma_malloc
rat*
memory::acquire<rat>(const uword n_elem) {
	if(n_elem == 0) return nullptr;

	arma_conform_check(
		( size_t(n_elem) > (std::numeric_limits<size_t>::max() / sizeof(rat)) ),
		"arma::memory::acquire(): requested size is too large"
	);
	rat* out_memptr = new(std::nothrow) rat[n_elem];
	arma_check_bad_alloc( (out_memptr == nullptr), "arma::memory::acquire(): out of memory" );
	return out_memptr;
}

template<>
arma_inline
void
memory::release<rat>(rat* mem) {
	delete [] mem;
}

template<>
arma_inline
void
memory::release<const rat>(const rat* mem) {
	delete [] mem;
}

// armadillo uses memcpy to copy memory, we can't do that so we manually copy
//
template<>
arma_hot
arma_inline
void
arrayops::copy<rat>(rat* dest, const rat* src, const uword n_elem) {
	for(uint32_t i = 0; i < n_elem; i++)
		dest[i] = src[i];
}

// conversion from strings 
//
template<>
inline
bool
diskio::convert_token(rat& val, const std::string& token) {
	val = token;
	return true;
}

// This is called from constructors taking fill:: arguments (even if fill::randn is never used).
// And _only_ when ARMA_USE_EXTERN_CXX11_RNG == false (which is the case in MSVC)
template<>
inline
void
arma_rng_cxx03::randn_dual_val(rat&, rat&) {
	throw "not implemented";
}

// The generic Datum<eT> computes log_min/log_max via std::log, which is not available for rat.
// So we specialize the class and declare the math constants (the physical ones are skipped).
//
template<>
struct Datum<rat>
  {
  static const rat pi;           //!< ratio of any circle's circumference to its diameter
  static const rat tau;          //!< ratio of any circle's circumference to its radius (replacement of 2*pi)
  static const rat e;            //!< base of the natural logarithm
  static const rat euler;        //!< Euler's constant, aka Euler-Mascheroni constant
  static const rat gratio;       //!< golden ratio
  static const rat sqrt2;        //!< square root of 2
  static const rat sqrt2pi;      //!< square root of 2*pi
  static const rat log_sqrt2pi;  //!< log of square root of 2*pi
  static const rat eps;          //!< the difference between 1 and the least value greater than 1 that is representable
  static const rat log_min;      //!< log of the minimum representable value
  static const rat log_max;      //!< log of the maximum representable value
  static const rat nan;          //!< "not a number"
  static const rat inf;          //!< positive infinity
  static const rat pos_inf;      //!< positive infinity
  static const rat neg_inf;      //!< negative infinity
  };

inline const rat Datum<rat>::pi          = rat(Datum<double>::pi);
inline const rat Datum<rat>::tau         = rat(Datum<double>::tau);
inline const rat Datum<rat>::e           = rat(Datum<double>::e);
inline const rat Datum<rat>::euler       = rat(Datum<double>::euler);
inline const rat Datum<rat>::gratio      = rat(Datum<double>::gratio);
inline const rat Datum<rat>::sqrt2       = rat(Datum<double>::sqrt2);
inline const rat Datum<rat>::sqrt2pi     = rat(Datum<double>::sqrt2pi);
inline const rat Datum<rat>::log_sqrt2pi = rat(Datum<double>::log_sqrt2pi);
inline const rat Datum<rat>::eps         = std::numeric_limits<rat>::epsilon();
inline const rat Datum<rat>::log_min     = rat(Datum<double>::log_min);
inline const rat Datum<rat>::log_max     = rat(Datum<double>::log_max);
inline const rat Datum<rat>::nan         = priv::Datum_helper::nan<rat>();
inline const rat Datum<rat>::inf         = priv::Datum_helper::pos_inf<rat>();
inline const rat Datum<rat>::pos_inf     = priv::Datum_helper::pos_inf<rat>();
inline const rat Datum<rat>::neg_inf     = priv::Datum_helper::neg_inf<rat>();

} // namespace arma
