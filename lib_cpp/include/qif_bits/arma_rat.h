
namespace arma {

using qif::rat;

// register rational<1> as a real type.
template<>
struct arma_real_only<rat> {
	typedef rat result;
};

// use direct_dot_arma (generic dot product implementation) for direct_dot<rat>
//
template<>
arma_hot inline rat op_dot::direct_dot<rat>(const uword n_elem, const rat* const A, const rat* const B) {
	return op_dot::direct_dot_arma<rat>(n_elem, A, B);
}

// for abs
//
template<>
arma_inline rat eop_aux::arma_abs<rat>(const rat x) {
	return mppp::abs(x);
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
	arma_debug_check(
		( size_t(n_elem) > (std::numeric_limits<size_t>::max() / sizeof(rat)) ),
		"arma::memory::acquire(): requested size is too large"
	);
	return ( new(std::nothrow) rat[n_elem] );
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
arma_rng_cxx98::randn_dual_val(rat&, rat&) {
	throw "not implemented";
}

// Datum<eT> declares various constexpr constants, but a constexpr cannot have type rat.
// So we specialize the class and declare the same constants as inline const.
//
template<>
class Datum<rat>
  {
  public:
  static const rat pi;       //!< ratio of any circle's circumference to its diamrater
  static const rat e;        //!< base of the natural logarithm
  static const rat euler;    //!< Euler's constant, aka Euler-Mascheroni constant
  static const rat gratio;   //!< golden ratio
  static const rat sqrt2;    //!< square root of 2
  static const rat sqrt2pi;  //!< square root of 2*pi
  static const rat eps;      //!< the difference bratween 1 and the least value greater than 1 that is representable
  static const rat log_min;  //!< log of the minimum representable value
  static const rat log_max;  //!< log of the maximum representable value
  static const rat nan;      //!< "not a number"
  static const rat inf;      //!< infinity 
  };
  
inline const rat Datum<rat>::pi        = rat(Datum<double>::pi);
inline const rat Datum<rat>::e         = rat(Datum<double>::e);
inline const rat Datum<rat>::euler     = rat(Datum<double>::euler);
inline const rat Datum<rat>::gratio    = rat(Datum<double>::gratio);
inline const rat Datum<rat>::sqrt2     = rat(Datum<double>::sqrt2);
inline const rat Datum<rat>::sqrt2pi   = rat(Datum<double>::sqrt2pi);
inline const rat Datum<rat>::eps       = std::numeric_limits<rat>::epsilon();
inline const rat Datum<rat>::log_min   = rat(Datum<double>::log_min);
inline const rat Datum<rat>::log_max   = rat(Datum<double>::log_max);
inline const rat Datum<rat>::nan       = priv::Datum_helper::nan<rat>();
inline const rat Datum<rat>::inf       = priv::Datum_helper::inf<rat>();

} // namespace arma