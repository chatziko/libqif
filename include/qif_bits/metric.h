
namespace metric {

// use eT_def as default type for R (return type of metrics)
typedef eT_def R_def;

// NOTE ON chainable()
// chainable(a, b) should true only if there is a "tight chain" between a, b, i.e. there
// is a chain a = e_1, ..., e_n = b such that d(a, b) = sum_i d(e_i, e_i+1).
// Returning false is safe if we cannot guarantee this chain, and is the default if chainable() is not defined.


// Euclidean on arithmetic T except uint
template<typename R = R_def, typename T, EnableIf<std::is_arithmetic<T>> = _, DisableIf<std::is_same<T, uint>> = _>
Metric<R, T>
euclidean() {
	return [](const T& a, const T& b) -> R {
		return abs_diff(a, b);
	};
}

// on uint's (discrete space, particularly useful for measuring distances
// between channel inputs), all non-consecutive elements are chainable
//
template<typename R = R_def, typename T, EnableIf<std::is_same<T, uint>> = _>
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

template<typename R = R_def, typename T>
Metric<R, T>
discrete() {
	return [](const T& a, const T& b) -> R {
		return R(equal(a, b) ? 0 : 1);
	};
}

template<typename R = R_def, typename T>
Metric<R, T>
mult_reals() {
	return [](const T& a, const T& b) -> R {
		return abs_diff(std::log(a), std::log(b));
	};
}

template<typename R = R_def, typename T>
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
template<typename R = R_def, typename T>
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
template<typename R = R_def, typename T>
Metric<R, T>
max(Metric<R, T> d1, Metric<R, T> d2) {
	Metric<R, T> d = [d1, d2](const T& a, const T& b) -> R {
		R r1 = d1(a, b);
		R r2 = d2(a, b);
		return less_than(r1, r2) ? r2 : r1;
	};
	return d;
}

// swap a and b
//
template<typename R = R_def, typename T>
Metric<R, T>
mirror(Metric<R, T> d) {
	Metric<R, T> d2 = [d](const T& a, const T& b) -> R {
		return d(b, a);
	};
	d2.chainable = [d](const T& a, const T& b) -> bool {
		return d.chainable(b, a);
	};
	return d2;
}

// min( d(a,b), thres ) (always a metric)
//
template<typename R = R_def, typename T>
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
template<typename R = R_def, typename T>
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
template<typename R = R_def, typename T>
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
template<typename R = R_def, typename T, EnableIf<is_Point<T>> = _, DisableIf<std::is_same<T, Point<uint>>> = _>
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
template<typename R = R_def, typename T, EnableIf<std::is_same<T, Point<uint>>> = _>
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

template<typename R = R_def, typename T, EnableIf<is_Point<T>> = _, DisableIf<std::is_same<T, Point<uint>>> = _>
Metric<R, T>
manhattan() {
	return [](const T& a, const T& b) -> R {
		return abs_diff(a.x, b.x) + abs_diff(a.y, b.y);
	};
}

// Mahattan distance on discrete (uint) cartesian Points. All points are chainable except
// those whose index differs by at most 1
//
template<typename R = R_def, typename T, EnableIf<std::is_same<T, Point<uint>>> = _>
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

template<typename R = R_def>
Metric<R, uint>
from_distance_matrix(Mat<R>& M) {
	return [&M](const uint& a, const uint& b) -> R {
		return M(a, b);
	};
}

template<typename R = R_def>
Mat<R>
to_distance_matrix(Metric<R, uint> d, uint size) {

	Mat<R> Dist(size, size);

	for(uint i = 0; i < size; i++)
		for(uint j = 0; j <= i; j++)
			Dist(i, j) = Dist(j, i) = d(i, j);

	return Dist;
}


// Compose a metric on T2, and function f:T1->T2, into a metric on T1
//
template<typename R = R_def, typename T1, typename T2>
Metric<R, T1>
compose(Metric<R,T2> d, std::function<T2(T1)> f) {
	Metric<R, T1> d2 = [=](const T1& a, const T1& b) -> R {
		return d( f(a), f(b) );
	};
	d2.chainable = [=](const T1& a, const T1& b) -> bool {
		return d.chainable( f(a), f(b) );		// maybe not always correct?
	};
	return d2;
}
// same but f takes a reference
template<typename R = R_def, typename T1, typename T2>
Metric<R, T1>
compose(Metric<R,T2> d, std::function<T2(const T1&)> f) {
	Metric<R, T1> d2 = [=](const T1& a, const T1& b) -> R {
		return d( f(a), f(b) );
	};
	d2.chainable = [=](const T1& a, const T1& b) -> bool{
		return d.chainable( f(a), f(b) );		// maybe not always correct?
	};
	return d2;
}
// same with two functions f1,f2
template<typename R = R_def, typename T1, typename T2>
Metric<R, T1>
compose(Metric<R,T2> d, std::function<T2(T1)> f1, std::function<T2(T1)> f2) {
	Metric<R, T1> d2 = [=](const T1& a, const T1& b) -> R {
		return d( f1(a), f2(b) );
	};
	d2.chainable = [=](const T1& a, const T1& b) -> bool {
		return d.chainable( f1(a), f2(b) );		// maybe not always correct?
	};
	return d2;
}
// same but f takes a reference
template<typename R = R_def, typename T1, typename T2>
Metric<R, T1>
compose(Metric<R,T2> d, std::function<T2(const T1&)> f1, std::function<T2(const T1&)> f2) {
	Metric<R, T1> d2 = [=](const T1& a, const T1& b) -> R {
		return d( f1(a), f2(b) );
	};
	d2.chainable = [=](const T1& a, const T1& b) -> bool{
		return d.chainable( f1(a), f2(b) );		// maybe not always correct?
	};
	return d2;
}


// OBSOLETE. Use compose(euclidean(), cell_to_point(...))
// 
// transform a metric on Point<uint> to a metric on indexes (uint) on a grid of the given width.
// The cell of index 0 is (0,0) (bottom left), and the cell of index i
// is (i%width, i/width).
// The step is fixed to 1, for a different step just scale the resulting metric
//
// Note: The metric space type (denoted by T above) is fixed to uint
//
template<typename R = R_def>
Metric<R, uint>
grid(uint width, Metric<R, Point<uint>> d = metric::euclidean<R, Point<uint>>()) {
	return compose<R, uint, Point<uint>>(d, geo::cell_to_point<uint>(width));
}


//////////////////////// METRICS ON PROBABILITY DISTRIBUTIONS ////////////////////////

template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
total_variation() {
	static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

	return [](const T& a, const T& b) -> R {
		if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

		return arma::accu(arma::abs(a - b)) / 2;
	};
}

template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
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

// the quasi metric that gives additive pi-capacity for 1-spanning Vg's
//
template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
convex_separation_quasi() {
	static_assert(std::is_same<R, typename T::elem_type>::value, "result and prob element type should be the same");

	return [](const T& a, const T& b) -> R {
		if(a.n_cols != b.n_cols) throw std::runtime_error("size mismatch");

		R res = R(0);
		for(uint i = 0; i < a.n_cols; i++) {
			if(equal(a.at(i), R(0)))
				continue;

			R d = 1 - b.at(i) / a.at(i);
			if(d > res)
				res = d;
		}
		return res;
	};
}

// the old symmetric variant, less useful
//
template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
convex_separation() {
	auto q = convex_separation_quasi<R, T>();
	return max(q, mirror(q));
}

// kantorovich through linear programming
//
template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
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
		lp.maximize = false;
		uint n = a.n_cols;
		auto vars = lp.make_vars(n, n, 0, infinity<R>());

		// minimize  sum_ij x_ij d(i,j)
		for(uint i = 0; i < n; i++)
		for(uint j = 0; j < n; j++)
			lp.set_obj_coeff(vars[i][j], d(i, j));

		// s.t.      sum_j x_ij = a[i]     forall i
		for(uint i = 0; i < n; i++) {
			auto con = lp.make_con(a[i], a[i]);
			for(uint j = 0; j < n; j++)
				lp.set_con_coeff(con, vars[i][j], 1);
		}

		//           sum_i x_ij = b[j]     forall j
		for(uint j = 0; j < n; j++) {
			auto con = lp.make_con(b[j], b[j]);
			for(uint i = 0; i < n; i++)
				lp.set_con_coeff(con, vars[i][j], 1);
		}

		// for floats, set the default solver to GLPK, CLP failed tests (numerical instability?)
		if(std::is_same<R,float>::value && lp.solver == lp::Solver::AUTO)
			lp.solver = lp::Solver::GLPK;

		if(!lp.solve())
			throw std::runtime_error(std::string("Kantorovich program failed: ") +
									 (lp.status == lp::Status::INFEASIBLE ? "infeasible" : "unbounded"));
		return lp.objective();
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
template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
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
template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
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
template<typename R = R_def, typename T, EnableIf<is_Prob<T>> = _>
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

		auto vars_x = lp.make_vars(n, n, R(0), infinity<R>());
		auto vars_r = lp.make_vars(n, R(0), infinity<R>());
		auto var_z = lp.make_var(0, infinity<R>());

		lp.maximize = false;
		lp.set_obj_coeff(var_z, 1);

		// s.t.      sum_j x_ij             - r_i           = a[i]  forall i
		for(uint i = 0; i < n; i++) {
			auto con = lp.make_con(a[i], a[i]);

			lp.set_con_coeff(con, vars_r[i], -1);

			for(uint j = 0; j < n; j++) {
				// when m(i,j) == inf, it forces the x_ij variable to be 0. We capture the same effect
				// by leaving the coefficients of x_ij to be 0. So all (in)equalities behave as if x_ij was fixed to 0.
				R dist = d(i, j);
				if(!equal(dist, infinity<R>()))
					lp.set_con_coeff(con, vars_x[i][j], 1);
			}
		}

		//           sum_i x_ij exp(d(i,j)) - r_j - b[j] z <= 0     forall j
		for(uint j = 0; j < n; j++) {
			auto con = lp.make_con(-infinity<R>(), 0);

			lp.set_con_coeff(con, vars_r[j], -1);
			lp.set_con_coeff(con, var_z, - b[j]);

			for(uint i = 0; i < n; i++) {
				// ignore inf again, see comment above
				R dist = d(i, j);
				if(!equal(dist, infinity<R>()))
					lp.set_con_coeff(con, vars_x[i][j], std::exp(dist));
			}
		}

		// an infeasible problem means that there's no finite z to satisfy it, so the distance is infinite
		return lp.solve()
			? std::log(lp.objective())
			: infinity<R>();
	};

	return [one_side](const T& a, const T& b) -> R {
		// solve both programs and keep the max
		return std::max(one_side(a, b), one_side(b, a));
	};
}

template<typename R = R_def, typename A, typename B, typename D>
bool
is_lipschitz(std::function<B(A)> f, Metric<R,A> da, Metric<R,B> db, const D& domain) {

	for(auto a1 = domain.begin(); a1 != domain.end(); a1++) {
		auto a2 = a1;
		for(++a2; a2 != domain.end(); a2++) {
			// chainable elements are redundant to check
			if(da.chainable(*a1, *a2)) continue;

			if(!less_than_or_eq(db(f(*a1), f(*a2)), da(*a1, *a2)))
				return false;
		}
	}
	return true;
}

template<typename R = R_def, typename A, typename B, typename D>
R
lipschitz_constant(std::function<B(A)> f, Metric<R,A> da, Metric<R,B> db, const D& domain) {

	R res(0);
	for(auto a1 = domain.begin(); a1 != domain.end(); a1++) {
		auto a2 = a1;
		for(++a2; a2 != domain.end(); a2++) {
			// chainable elements are redundant to check
			if(da.chainable(*a1, *a2)) continue;

			R ratio = db(f(*a1), f(*a2)) / da(*a1, *a2);
			if(less_than(res, ratio))
				res = ratio;
		}
	}
	return res;
}

} // namespace metric

template<typename R = metric::R_def, typename T>
Metric<R, T>
operator*(R coeff, const Metric<R, T>& d) {
    return metric::scale(d, coeff);
}

