
namespace metric {

// NOTE ON chainable()
// chainable(a, b) should true only if there is a "tight chain" between a, b, i.e. there
// is a chain a = e_1, ..., e_n = b such that d(a, b) = sum_i d(e_i, e_i+1).
// Returning false is safe if we cannot guarantee this chain, and is the default if chainable() is not defined.


// Euclidean on arithmetic T except uint
template<typename R, typename T, EnableIf<std::is_arithmetic<T>> = _, DisableIf<std::is_same<T, uint>> = _>
Metric<R, T>
euclidean() {
	return [](const T& a, const T& b) -> R {
		return abs_diff(a, b);
	};
}

// on uint's (discrete space, particularly useful for measuring distances
// between channel inputs), all non-consecutive elements are chainable
//
template<typename R, typename T, EnableIf<std::is_same<T, uint>> = _>
Metric<R, T>
euclidean() {
	Metric<R, T> d = [](const T& a, const T& b) -> R {
		return abs_diff(a, b);
	};
	d.chainable = [](const T& a, const T& b) -> bool {
		return abs_diff(a, b) != 1;
	};
	return d;
}

template<typename R, typename T>
Metric<R, T>
discrete() {
	return [](const T& a, const T& b) -> R {
		return R(equal(a, b) ? 0 : 1);
	};
}

template<typename R, typename T>
Metric<R, T>
scale(Metric<R, T> d, R coeff) {
	Metric<R, T> d2 = [d, coeff](const T& a, const T& b) -> R {
		// separate treatment of 0 allows to scale by infinity and still get d(x,x) == 0
		R r = d(a, b);
		if(r != R(0)) r *= coeff;
		return r;
	};
	d2.chainable = d.chainable;
	return d2;
}

// min of two metrics (technically not a metric)
//
template<typename R, typename T>
Metric<R, T>
min(Metric<R, T> d1, Metric<R, T> d2) {
	Metric<R, T> d = [d1, d2](const T& a, const T& b) -> R {
		R r1 = d1(a, b);
		R r2 = d2(a, b);
		return less_than(r1, r2) ? r1 : r2;
	};
	return d;
}

// max of two metrics (always a metric)
//
template<typename R, typename T>
Metric<R, T>
max(Metric<R, T> d1, Metric<R, T> d2) {
	Metric<R, T> d = [d1, d2](const T& a, const T& b) -> R {
		R r1 = d1(a, b);
		R r2 = d2(a, b);
		return less_than(r1, r2) ? r2 : r1;
	};
	return d;
}

// min( d(a,b), thres ) (always a metric)
//
template<typename R, typename T>
Metric<R, T>
threshold(Metric<R, T> d, R thres) {
	Metric<R, T> d2 = [d, thres](const T& a, const T& b) -> R {
		R r = d(a, b);
		return less_than(r, thres) ? r : thres;
	};
	d2.chainable = [d, thres](const T& a, const T& b) -> bool {
		// a,b are chainable in d2 if they are chainable in d and below the threshold
		return d.chainable(a, b) && less_than(d(a, b), thres);
	};
	return d2;
}

// 0 if below threshold t, 1 otherwise (technically not a metric)
// Kostas: is this used somewhere?
//
template<typename R, typename T>
Metric<R, T>
threshold_bin(Metric<R, T> d, R thres) {
	Metric<R, T> d2 = [d, thres](const T& a, const T& b) -> R {
		return less_than(d(a, b), thres) ? 0 : 1;
	};
	d2.chainable = [d, thres](const T& a, const T& b) -> bool {
		// a,b are chainable in d2 if they are chainable in d and below the threshold
		return d.chainable(a, b) && less_than(d(a, b), thres);
	};
	return d2;
}

// inf if above threshold t, same as d otherwise (technically not a metric)
//
template<typename R, typename T>
Metric<R, T>
threshold_inf(Metric<R, T> d, R thres) {
	Metric<R, T> d2 = [d, thres](const T& a, const T& b) -> R {
		R res = d(a, b);
		return less_than_or_eq(res, thres) ? res : inf;
	};
	d2.chainable = [d, thres](const T& a, const T& b) -> bool {
		// a,b are chainable in d2 if they are chainable in d and below the threshold
		return d.chainable(a, b) && less_than_or_eq(d(a, b), thres);
	};
	return d2;
}

// Euclidean distance on non-uint Points
//
template<typename R, typename T, EnableIf<is_Point<T>> = _, DisableIf<std::is_same<T, Point<uint>>> = _>
Metric<R, T>
euclidean() {
	return [](const T& a, const T& b) -> R {
		auto v1 = a.x - b.x,
			 v2 = a.y - b.y;
		return std::sqrt(v1*v1 + v2*v2);
	};
}

// Euclidean distance on discrete (uint) cartesian Points. The only chainable points
// are on the same line or diagonal, and at index difference more than one
//
template<typename R, typename T, EnableIf<std::is_same<T, Point<uint>>> = _>
Metric<R, T	>
euclidean() {
	Metric<R, T> d = [](const T& a, const T& b) -> R {
		uint v1 = abs_diff(a.x, b.x),
			 v2 = abs_diff(a.y, b.y);
		return std::sqrt(v1*v1 + v2*v2);
	};
	d.chainable = [](const T& a, const T& b) -> bool {
		uint v1 = abs_diff(a.x, b.x),
			 v2 = abs_diff(a.y, b.y);
		return (v1 == 0 || v2 == 0 || v1 == v2) && (v1 > 1 || v2 > 1);
	};
	return d;
}

template<typename R, typename T, EnableIf<is_Point<T>> = _, DisableIf<std::is_same<T, Point<uint>>> = _>
Metric<R, T>
manhattan() {
	return [](const T& a, const T& b) -> R {
		return abs_diff(a.x, b.x) + abs_diff(a.y, b.y);
	};
}

// Mahattan distance on discrete (uint) cartesian Points. All points are chainable except
// those whose index differs by at most 1
//
template<typename R, typename T, EnableIf<std::is_same<T, Point<uint>>> = _>
Metric<R, T>
manhattan() {
	Metric<R, T> d = [](const T& a, const T& b) -> R {
		return abs_diff(a.x, b.x) + abs_diff(a.y, b.y);
	};
	d.chainable = [](const T& a, const T& b) -> bool {
		return abs_diff(a.x, b.x) > 1 || abs_diff(a.y, b.y) > 1;
	};
	return d;
}

template<typename R>
Metric<R, uint>
from_distance_matrix(Mat<R>& M) {
	return [&M](const uint& a, const uint& b) -> R {
		return M(a, b);
	};
}

template<typename R>
Mat<R>
to_distance_matrix(Metric<R, uint> d, uint size) {

	Mat<R> Dist(size, size);

	for(uint i = 0; i < size; i++)
		for(uint j = 0; j <= i; j++)
			Dist(i, j) = Dist(j, i) = d(i, j);

	return Dist;
}

// transform a metric on Point<uint> to a metric on indexes (uint) on a grid of the given width.
// The cell of index 0 is (0,0) (bottom left), and the cell of index i
// is (i%width, i/width).
// The step is fixed to 1, for a different step just scale the resulting metric
//
// Note: The metric space type (denoted by T above) is fixed to uint
//
template<typename R>
Metric<R, uint>
grid(uint width, Metric<R, Point<uint>> d = metric::euclidean<R, Point<uint>>()) {
	Metric<R, uint> d2 = [width, d](const uint& a, const uint& b) -> R {
		return d( Point<uint>(a%width, a/width), Point<uint>(b%width, b/width) );
	};
	d2.chainable = [width, d](const uint& a, const uint& b) -> bool {
		return d.chainable( Point<uint>(a%width, a/width), Point<uint>(b%width, b/width) );
	};
	return d2;
}


//////////////////////// METRICS ON PROBABILITY DISTRIBUTIONS ////////////////////////

template<typename R, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
total_variation() {
	static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

	return [](const T& a, const T& b) -> R {
		if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

		return arma::accu(arma::abs(a - b)) / 2;
	};
}

template<typename R, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
mult_total_variation() {
	static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

	return [](const T& a, const T& b) -> R {
		if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

		R res = R(0);
		for(uint i = 0; i < a.n_cols; i++) {
			R log_a = std::log(a(i)),
			  log_b = std::log(b(i));
			bool inf_a = equal(log_a, -infinity<R>()),
				 inf_b = equal(log_b, -infinity<R>());

			if(inf_a && inf_b)					// both are zero, no diff
				continue;
			else if(inf_a || inf_b)				// only one is zero, diff is infty
				return infinity<R>();
			else {
				R diff = std::abs(log_a - log_b);
				if(less_than(res, diff))
					res = diff;
			}
		}
		return res;
	};
}

template<typename R, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
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
//
template<typename R, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
kantorovich_lp(Metric<R, uint> d) {
	static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

	return [d](const T& a, const T& b) -> R {
		if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

		// linear program for the transportation problem
		// vars:     x_ij >= 0             forall i,j
		// minimize  sum_ij x_ij d(i,j)
		// s.t.      sum_j x_ij = a[i]     forall i
		//           sum_i x_ij = b[j]     forall j
		//
		lp::LinearProgram<R> lp;
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
									 (lp.status == lp::status_t::infeasible ? "infeasible" : "unbounded"));
		return lp.optimum();
	};
}

// kantorovich using the FastEMD algorithm from:
// http://ofirpele.droppages.com//ICCV2009.pdf
// https://dl.dropboxusercontent.com/s/i5g3a8tqsm2hcpl/FastEMD-3.1.zip?dl=0
//
// It assumes that d is a metric!
// Only T = double is supported. (emd_hat_fb_metric works with: int, long int, long long int, double)
//
// TODO: the double implementation is via ints, see: emd_hat_impl<double,FLOW_TYPE> in emd_hat_impl.hpp.
//       we could do the same to support float, rat
// 
// (Note: an alternative algorithm: jorlin.scripts.mit.edu/docs/publications/26-faster strongly polynomial.pdf)
//
template<typename R, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
kantorovich_fastemd(Metric<R, uint> d) {
	static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

	return [d](const T& a, const T& b) -> R {
		if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

		typedef std::vector<R> vec;

		// create distance matrix
		auto Dist = to_distance_matrix<R>(d, a.n_cols);

		std::vector<vec> Dist_v( a.n_cols );
		for(uint i = 0; i < a.n_cols; i++)
			Dist_v[i] = arma::conv_to<vec>::from(Dist.row(i));
		Dist.clear();

		auto a_v = arma::conv_to<vec>::from(a);
		auto b_v = arma::conv_to<vec>::from(b);

		return fastemd::emd_hat_gd_metric<R>()(a_v, b_v, Dist_v);
	};
}

//  kantorovich. use FastEMD for doubles, LP for all others
//
template<typename R, typename T, EnableIf<is_Prob<T>> = _>
inline
Metric<R, T>
kantorovich(Metric<R, uint> d) {
	return kantorovich_lp<R, T>(d);
}

template<>
inline
Metric<double, prob>
kantorovich(Metric<double, uint> d) {
	return kantorovich_fastemd<double, prob>(d);
}

// multiplicative kantorovich through linear programming
//
template<typename R, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
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
		lp::LinearProgram<R> lp;
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

} // namespace metric

template<typename R, typename T>
Metric<R, T>
operator*(R coeff, const Metric<R, T>& d) {
    return metric::scale(d, coeff);
}

