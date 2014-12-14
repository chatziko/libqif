#ifndef _QIF_Metric_h_
#define _QIF_Metric_h_

#include <type_traits>
#include <exception>
#include "types.h"
#include "aux.h"
#include "LinearProgram.h"


// R: result type, what we measure distances in
// T: metric space type, what elements we measure the distance of

template<typename R, typename T> using Metric = std::function<R(const T&, const T&)>;

class metrics {
	public:

	template<typename R, typename T, EnableIf<std::is_arithmetic<T>>...>
	static Metric<R, T>
	euclidean() {
		return [](const T& a, const T& b) -> R {
			return abs_diff(a, b);
		};
	}

	template<typename R, typename T>
	static Metric<R, T>
	discrete(R val = R(1)) {
		return [val](const T& a, const T& b) -> R {
			return equal(a, b) ? R(0) : val;
		};
	}

	template<typename R, typename T>
	static Metric<R, T>
	scale(Metric<R, T> d, R coeff) {
		return [d, coeff](const T& a, const T& b) -> R {
			return coeff * d(a, b);
		};
	}

	template<typename R, typename T, EnableIf<is_Point<T>>...>
	static Metric<R, T>
	euclidean() {
		return [](const T& a, const T& b) -> R {
			auto v1 = a.x - b.x,
				 v2 = a.y - b.y;
			return std::sqrt(v1*v1 + v2*v2);
		};
	}

	template<typename R, typename T, EnableIf<is_Point<T>>...>
	static Metric<R, T>
	manhattan() {
		return [](const T& a, const T& b) -> R {
			return abs_diff(a.x, b.x) + abs_diff(a.y, b.y);
		};
	}

	// transform a metric on points to a metric on indexes on a grid of the given width.
	// The cell of index 0 is (0,0) (bottom left), and the cell of index i
	// is (i%width, i/width).
	// The step is fixed to 1, for a different step just scale the resulting metric
	//
	// Note: PT is the point type. The metric space type (denoted by T above) is fixed to uint
	//
	template<typename R, typename PT, EnableIf<is_Point<PT>>...>
	static Metric<R, uint>
	grid(uint width, Metric<R, PT> d = metrics::euclidean<R, PT>()) {
		return [width, d](const uint& a, const uint& b) -> R {
			return d( PT(a%width, a/width), PT(b%width, b/width) );
		};
	}


	//////////////////////// METRICS ON PROBABILITY DISTRIBUTIONS ////////////////////////

	template<typename R, typename T, EnableIf<is_Prob<T>>...>
	static Metric<R, T>
	total_variation() {
		static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

		return [](const T& a, const T& b) -> R {
			if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

			return arma::accu(arma::abs(a - b)) / 2;
		};
	}

	template<typename R, typename T, EnableIf<is_Prob<T>>...>
	static Metric<R, T>
	mult_total_variation() {
		static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

		return [](const T& a, const T& b) -> R {
			if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

			R res = R(0);
			for(uint i = 0; i < a.n_cols; i++) {
				R ai = a(i),
				  bi = b(i);
				bool ai_is_zero = equal(ai, R(0)),
					 bi_is_zero = equal(bi, R(0));

				if(ai_is_zero && bi_is_zero)			// both are the same, no diff
					continue;
				else if(ai_is_zero || bi_is_zero)		// only one is zero, diff is infty
					return infinity<R>();
				else {
					R diff = std::abs(std::log(ai) - std::log(bi));
					if(less_than(res, diff))
						res = diff;
				}
			}
			return res;
		};
	}

	template<typename R, typename T, EnableIf<is_Prob<T>>...>
	static Metric<R, T>
	bounded_entropy_distance() {
		static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

		return [](const T& a, const T& b) -> R {
			if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

			R res = R(0);
			for(uint i = 0; i < a.n_cols; i++) {
				R m = std::max(a.at(i), b.at(i));
				if(equal(m, R(0)))
					continue;

				R d = abs(a.at(i) - b.at(i)) / m;
				if(d > res)
					res = d;
			}
			return res;
		};
	}

	// kantorovich through linear programming
	// WISHLIST: faster method: jorlin.scripts.mit.edu/docs/publications/26-faster strongly polynomial.pdf
	//
	template<typename R, typename T, EnableIf<is_Prob<T>>...>
	static Metric<R, T>
	kantorovich(Metric<R, uint> d) {
		static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

		return [d](const T& a, const T& b) -> R {
			if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

			// linear program for the transportation problem
			// vars:     x_ij >= 0             forall i,j
			// minimize  sum_ij x_ij d(i,j)
			// s.t.      sum_j x_ij = a[i]     forall i
			//           sum_i x_ij = b[j]     forall j
			//
			LinearProgram<R> lp;
			uint n      = a.n_cols,
				 n_vars = n * n,
				 n_cons = 2 * n;

			lp.maximize = false;
			lp.A = arma::zeros<Mat<R>>(n_cons, n_vars);
			lp.b = arma::trans(arma::join_rows(a, b));
			lp.c.set_size(n_vars);
			lp.sense.set_size(n_cons);
			lp.sense.fill('=');

			for(uint i = 0; i < n; i++) {
				for(uint j = 0; j < n; j++) {
					uint xij_id = i * n + j;	// index of the variable x_ij in A,c

					lp.c(xij_id) = d(i, j);		// coefficient

					lp.A(i,     xij_id) = R(1);	// var x_ij participates in the i-th constraint of the first family
					lp.A(j + n, xij_id) = R(1);	// and the j-th constraint of the second
				}
			}

			if(!lp.solve())
				throw std::runtime_error(std::string("Kantorovich program failed: ") +
										 (lp.status == LinearProgram<R>::status_t::infeasible ? "infeasible" : "unbounded"));
			return lp.optimum();
		};
	}

	// multiplicative kantorovich through linear programming
	//
	template<typename R, typename T, EnableIf<is_Prob<T>>...>
	static Metric<R, T>
	mult_kantorovich(Metric<R, uint> d) {
		// we need to solve 2 linear programs, with the role of a, b exchanged. one_side solves one of them.
		//
		auto one_side = [d](const T& a, const T& b) -> R {
			static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

			uint n = a.n_cols;
			if(n != b.n_cols) throw std::runtime_error("size mismatch");

			// linear program for the transportation problem
			// vars:     x_ij, r_i, z >= 0    forall i,j
			// minimize  z
			// s.t.      sum_j x_ij             - r_i           = a[i]  forall i
			//           sum_i x_ij exp(d(i,j)) - r_j - b[j] z <= 0     forall j
			//
			LinearProgram<R> lp;
			uint n_vars = n * n + n + 1,
				 n_cons = 2 * n;

			lp.maximize = false;
			lp.A = arma::zeros<Mat<R>>(n_cons, n_vars);
			lp.b = arma::join_cols(a.t(), arma::zeros<Col<R>>(n));
			lp.c = arma::zeros<Col<R>>(n_vars);
			lp.c(n_vars - 1) = R(1);
			lp.sense.set_size(n_cons);
			lp.sense.rows(0, n-1     ).fill('=');
			lp.sense.rows(n, n_cons-1).fill('<');

			lp.A.col(n_vars-1).rows(n, n_cons-1) = - b.t();

			for(uint i = 0; i < n; i++) {
				uint ri_id = n * n + i;			// index of the variable r_i in A,c
				lp.A(i,     ri_id) = R(-1);		// var r_i participates in the i-th constraint of the first family
				lp.A(i + n, ri_id) = R(-1);		// and the i-th constraint of the second

				for(uint j = 0; j < n; j++) {
					uint xij_id = i * n + j;	// index of the variable x_ij in A,c

					// when m(i,j) == inf, it forces the x_ij variable to be 0. We capture the same effect
					// by leaving the A coefficients of x_ij to be 0. So all (in)equalities behave as if x_ij was fixed to 0.
					R dist = d(i, j);
					if(equal(dist, infinity<R>())) continue;

					lp.A(i,     xij_id) = R(1);					// var x_ij participates in the i-th constraint of the first family
					lp.A(j + n, xij_id) = R(std::exp(dist));	// and the j-th constraint of the second
				}
			}

			// an infeasible problem means that there's no finite z to satisfy it, so the distance is infinite
			return lp.solve()
				? std::log(lp.optimum())
				: infinity<R>();
		};

		return [one_side](const T& a, const T& b) -> R {
			// solve both programs and keep the max
			return std::max(one_side(a, b), one_side(b, a));
		};
	}
};


#endif
