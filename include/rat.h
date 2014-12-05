#ifndef _QIF_rat_h_
#define _QIF_rat_h_

#include <gmpxx.h>


typedef mpq_class rat;


// armadillo internal voodoo to make it play well with rat
//
namespace arma {
	// make armadillo support mpq_class as element type
	template<>
	struct is_supported_elem_type<mpq_class> {
		static const bool value = true;
	};

	// armadillo's memory::acquire/release uses malloc/free for speed, which overrides the constructor
	// and causes segfaults when used with rat. We override it to use C++'s new/delete
	//
	template<>
	inline
	arma_malloc
	rat*
	memory::acquire<rat>(const uword n_elem) {
//		std::cout << "ACQUIRE\n";
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
//		std::cout << "RELEASE\n";
		delete [] mem;
	}

	template<>
	arma_inline
	void
	memory::release<const rat>(const rat* mem) {
//		std::cout << "CONST RELEASE\n";
		delete [] mem;
	}

	// armadillo uses memcpy to copy memory, can't do that with rat
	//
	template<>
	arma_hot
	arma_inline
	void
	arrayops::copy<rat>(rat* dest, const rat* src, const uword n_elem) {
		for(uint i = 0; i < n_elem; i++)
			dest[i] = src[i];
	}
}

#endif
