
extern "C" {
	#include <gsl/gsl_sf.h>				// gsl_sf_lambert_Wm1
	#include <gsl/gsl_monte_miser.h>
}

#include "qif"


namespace qif {
namespace mechanism {

using std::max;
using std::min;

//double max_error = 0;
//double sum_error = 0;

double _coeff, _epsilon;					// params of palnar_laplace_pdf, need to be global to pass function pointer to gsl_monte_miser_integrate


double planar_laplace_pdf(double *x, size_t, void *) {
	return _coeff * exp( - _epsilon * sqrt(x[0] * x[0] + x[1] * x[1]) );
}

double inverse_cumulative_gamma(double epsilon, double p) {
	double x = (p-1) / M_El;
	return - (gsl_sf_lambert_Wm1(x) + 1) / epsilon;
}

// restrict val within the bound
int bound(int val, int min_val, int max_val) {
	return min(max(val, min_val), max_val);
}

gsl_rng *r;
gsl_monte_function L = { &planar_laplace_pdf, 2, 0 };

double integrate_laplace(double epsilon, const arma::vec& a, const arma::vec& b, int calls) {
	if(!r) {
		gsl_rng_env_setup();
		r = gsl_rng_alloc(gsl_rng_default);
	}

	_epsilon = epsilon;
	_coeff = (_epsilon * _epsilon) / (2 * M_PIl);

	double res, err;

	gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
	gsl_monte_miser_integrate(&L, a.colptr(0), b.colptr(0), 2, calls, r, s, &res, &err);
	gsl_monte_miser_free(s);

	// store errors
	//if(err > max_error) max_error = err;
	//sum_error += err;

	return res;
}

}} // namespace qif::mechanism

