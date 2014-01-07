//
// Fast implementation of the plannar laplace and tight constraint
// mechanisms on a grid
//
//
#include <iostream>
#include <algorithm>		// max

//#define ARMA_NO_DEBUG
#include <armadillo>

extern "C" {
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_sf.h>				// gsl_sf_lambert_Wm1
	#include <gsl/gsl_monte.h>
	#include <gsl/gsl_monte_miser.h>
}

using namespace std;
using namespace arma;


int bound(int val, int min_val, int max_val);
double euclidean(double* x1, double* x2);

double integrate_laplace(double epsilon, const vec& a, const vec& b, int calls);
double inverse_cumulative_gamma(double epsilon, double p);
mat grid_integration(int width, int height, double step, double epsilon);

mat plannar_laplace_grid(int width, int height, double step, double epsilon);
mat tight_grid(int width, int height, double step, double epsilon);

double smallest_epsilon(const mat& m, int width, double step, bool only_first_row = false);
double min_capacity(const mat& m);
double vulnerability(const mat& m, const vec& prior);


double max_error = 0;
double sum_error = 0;
int debug = 0;

int main (int argc, char** argv) {
	if(argc != 2) {
		cerr << "usage geo <size>\n";
		exit(-1);
	}

	int size = strtol(argv[1], 0, 10);
	double step = 1;

	vec uniform = vec(size*size);
	uniform.fill(1.0 / (size*size));

	cout << "size: " << size << ", step: " << step << "\n\n";
	
	for(double epsilon = 0.4; epsilon <= 1.3; epsilon += 0.02) {
		mat m;
		try {
			m = tight_grid(size, size, step, epsilon);
			double v_tight = vulnerability(m, uniform);
			m.reset();

			m = plannar_laplace_grid(size, size, step, epsilon);
			double v_lap = vulnerability(m, uniform);
			m.reset();

			cout << epsilon << ", " << v_tight << ", " << v_lap << "\n";

		} catch(...) {
			cout << epsilon << " not in range\n";
		}
	}


	/*
	if(argc != 5) {
		cout << "syntax: laplace <width> <height> <step> <epsilon>\n";
		return -1;
	}

	int width = strtol(argv[1], 0, 10);
	int height = strtol(argv[2], 0, 10);
	double step = strtod(argv[3], 0);
	double epsilon = strtod(argv[4], 0);

	mat m = plannar_laplace_grid(width, height, step, epsilon);
  	m.save(cout, arma_ascii);

	if(debug) {
		cerr << "done\n";
		cerr << "smallest epsilon: " << smallest_epsilon(m, width, step, true) << "\n";
		cerr << "max error: " << max_error << "\n";
		cerr << "sum error: " << sum_error << "\n";
	}
	*/

	return 0;
}

double min_capacity(const mat& m) {
	return log2(accu(max(m)));
}

double vulnerability(const mat& m, const vec& prior) {
	double s = 0;
	for(uint j = 0; j < m.n_cols; j++)
		s += max( m.col(j) % prior );
	return s;
}

double smallest_epsilon(const mat& m, int width, double step, bool only_first_row) {
	double max = 0;

	for(uint i1 = 0; i1 < m.n_rows; i1++) {
		double p1[2] = {
			(i1 % width) * step,
			(i1 / width) * step
		};

		for(uint i2 = i1+1; i2 < m.n_rows; i2++) {
			double p2[2] = {
				(i2 % width) * step,
				(i2 / width) * step
			};
			double dist = euclidean(p1, p2);

			for(uint j = 0; j < m.n_cols; j++) {
				double small = m(i1,j) < m(i2,j) ? m(i1,j) : m(i2,j);
				double large = m(i1,j) < m(i2,j) ? m(i2,j) : m(i1,j);

				// ignore very small probabilities, errors are more likely there
				if(small < 1e-5) continue;

				double eps = log(large/small)/dist;
				if(eps > max)
					max = eps;
			}
		}

		if(only_first_row)		// less reliable but fast
			break;
	}
	return max;
}

double euclidean(double* x1, double* x2) {
	double a = x1[0] - x2[0];
	double b = x1[1] - x2[1];
	return sqrt(a*a + b*b);
}

// restrict val within the bound
int bound(int val, int min_val, int max_val) {
	return min(max(val, min_val), max_val);
}

mat plannar_laplace_grid(int width, int height, double step, double epsilon) {

	int size = width * height;
	mat m = zeros<mat>(size, size);

	int cx = width - 1;
	int cy = height - 1;

	mat big_matrix = grid_integration(2*width-1, 2*height-1, step, epsilon);

	// we loop over the small grid, the current element is the chanel input
	//
	for(int input_x = 0; input_x < width; input_x++) {
		for(int input_y = 0; input_y < height; input_y++) {
			int input = input_y * width + input_x;

			// we put this element in the center of the big grid
			// we loop over the big grid, and add this probability to
			// the closest point in the small grid
			//
			for(int x = -cx; x <= cx; x++) {
				int output_x = bound(input_x + x, 0, width-1);

				for(int y = -cy; y <= cy; y++) {
					int output_y = bound(input_y + y, 0, height-1);
					int output = output_y * width + output_x;

					m(input,output) += big_matrix(cy+y, cx+x);
				}
			}
		}
	}

	return m;
}


mat grid_integration(int width, int height, double step, double epsilon) {
	if(width % 2 == 0 || height % 2 == 0) throw "width/height must by odd";
	if(width == 1 || height == 1) throw "width/height must be > 1";

	// integration is in theory more accurate when the area of integration is small.
	// In our case, we can simply scale down distances by adjusting epsilon
	// So we use a_step instead of the requested one, while properly adjusting
	// epsilon.
	// (Kostas: in practice I didn't see any difference, but still it seems a good idea to use a small step)
	//
	double a_step = 1;
	epsilon *= step / a_step;		// we _multiply_ epsilon by the same factor used to _divide_ step. i.e. when a_step < step then epsilon is increased

	mat m = zeros<mat>(height, width);	// rows = height, cols = width,  point (x,y) has index (cy+y, cx+x)

	int cx = (width-1)/2;									// corner x index
	int cy = (height-1)/2;									// corner y index
	const double grid_bound = (max(cx, cy) + 0.5) * a_step;	// boundary of the grid

	// we set the limit of integration to the distance in which 0.9999 of the
	// probability is included (and at least as big as the grid boundary)
	const double far_away = max(inverse_cumulative_gamma(epsilon, 0.9999), grid_bound);

	if(debug) {
		cerr << "integration epsilon: " << epsilon << "\n";
		cerr << "far_away: " << far_away << "\n";
		cerr << "grid_bound: " << grid_bound << "\n";
	}

	int c = 0;
	for(int i = 0; i <= cx; i++) {
		for(int j = i; j <= cy; j++) {
			// get integration rectangle
			vec a = { 
				(i - 0.5) * a_step,
				(j - 0.5) * a_step
			};
			vec b = {
				i == cx ? far_away : a(0) + a_step,
				j == cy ? far_away : a(1) + a_step,
			};

			// integrate. we use more calls when we go outside the grid boundary
			bool out_of_bound = b(0) > (grid_bound+0.001) || b(1) > (grid_bound+0.001);
			int calls = 50000 * (out_of_bound ? 10 : 1);

			double prob = integrate_laplace(epsilon, a, b, calls);
			c++;

			// due to symmetry, the value is copied in up to 8 cells.
			// a point (x,y) has index (cy+y, cx+x)
			// i.e. the (-cx,-cy) corner has index (0,0)
			m(cy+j,cx+i) =
			m(cy+j,cx-i) =
			m(cy-j,cx+i) =
			m(cy-j,cx-i) = prob;

			if(j <= cx && i <= cy)		// if width != height these might fall outside the grid
				m(cy+i,cx+j) =
				m(cy+i,cx-j) =
				m(cy-i,cx+j) =
				m(cy-i,cx-j) = prob;
		}
	}
	if(debug)
		cerr << "integrations " << c << "\n";

	return m;
}


double _coeff, _epsilon;

double plannar_laplace_pdf (double *x, size_t, void *) {
	return _coeff * exp( - _epsilon * sqrt(x[0] * x[0] + x[1] * x[1]) );
}

gsl_rng *r;
gsl_monte_function L = { &plannar_laplace_pdf, 2, 0 };

double integrate_laplace(double epsilon, const vec& a, const vec& b, int calls = 50000) {
	if(!r) {
		gsl_rng_env_setup();
		r = gsl_rng_alloc(gsl_rng_default);
	}

	_epsilon = epsilon;
	_coeff = (_epsilon * _epsilon) / (2 * M_PI);

	double res, err;

	gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
	gsl_monte_miser_integrate (&L, a.colptr(0), b.colptr(0), 2, calls, r, s, &res, &err);
	gsl_monte_miser_free (s);

	if(err > max_error) max_error = err;
	sum_error += err;

	return res;
}

double inverse_cumulative_gamma(double epsilon, double p) {
	double x = (p-1) / M_E;
	return - (gsl_sf_lambert_Wm1(x) + 1) / epsilon;
}




mat tight_grid(int width, int height, double step, double epsilon) {

	int size = width * height;

	mat phi = zeros<mat>(size, size);
	double e = exp(-epsilon);

	if(debug) cerr << "building phi\n";
	for(int i = 0; i < size; i++) {
		double x[2] = {
			(i % width) * step,
			(i / width) * step
		};
		phi(i,i) = 1;

		for(int j = i+1; j < size; j++) {
			double y[2] = {
				(j % width) * step,
				(j / width) * step
			};

			phi(i, j) = phi(j, i) = pow(e, euclidean(x, y));
		}
	}

	if(debug) cerr << "inverting\n";
	mat phi_inv = phi.i();

	mat uniform = zeros<mat>(1, size);
	uniform.fill(1.0/size);

	if(debug) cerr << "creating diag\n";
	mat diag = uniform * phi_inv;
	phi_inv.clear();		// free

	for(int i = 0; i < size; i++)
		if(diag(0,i) < -1e-4) {
			//cout << "negative diagonal i: " << i << ", diag(0,i): " << diag(0,i) << "\n";
			throw "negative";
		}

	// fill m
	//
	if(debug) cerr << "filling m\n";

	mat m = zeros<mat>(size, size);
	for(uint i = 0; i < m.n_rows; i++) {
		for(uint j = 0; j < m.n_rows; j++) {
			m(i,j) = size * diag(0,j) * phi(i,j);
		}
	}

	return m;
}


