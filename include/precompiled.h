#ifndef _QIF_precompiled_h_
#define _QIF_precompiled_h_

// this header gets precompiled, include here all heavy stuff

#include <functional>	// std::function
#include <armadillo>
#include <vector>
#include <cassert>
#include <set>
#include <list>
#include <map>
#include <unordered_map>
#include <iterator>		// needed by
#include <type_traits>	// range.hpp

#include <glpk.h>
#include <gmpxx.h>		// for rats

extern "C" {
	#include <gsl/gsl_sf.h>				// gsl_sf_lambert_Wm1
	#include <gsl/gsl_monte_miser.h>
}

#ifdef QIF_USE_ORTOOLS
#include <ortools/linear_solver/linear_solver.h>
#include <ortools/linear_solver/linear_solver.pb.h>
#endif

#endif
