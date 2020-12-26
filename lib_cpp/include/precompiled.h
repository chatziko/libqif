#ifndef _QIF_precompiled_h_
#define _QIF_precompiled_h_

// this header gets precompiled, include here all heavy stuff

#include <functional>	// std::function
#include <armadillo>
#include <vector>
#include <cassert>
#define _USE_MATH_DEFINES	// for PI
#include <cmath>
#include <set>
#include <list>
#include <map>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <iterator>		// needed by
#include <type_traits>	// range.hpp

#include <mp++/mp++.hpp>	// for rats

extern "C" {
	#include <gsl/gsl_sf.h>				// gsl_sf_lambert_Wm1
	#include <gsl/gsl_monte_miser.h>
}

#ifdef QIF_USE_ORTOOLS
#include <ortools/linear_solver/linear_solver.h>
#include <ortools/linear_solver/linear_solver.pb.h>
#endif

#ifdef QIF_USE_GLPK
#include <glpk.h>
#endif QIF_USE_GLPK

#endif
