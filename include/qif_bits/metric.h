namespace metric {

// NOTE ON ADJACENCY
// is_adjacent(a, b) should return false only if there is a path a = e_1, ..., e_n = b
// such that d(a, b) = sum_i d(e_i, e_i+1).
// Returning true is safe if we cannot guarantee this path, and is the default if is_adjacent is not defined.


// Euclidean on arithmetic T except uint
template<typename R, typename T, EnableIf<std::is_arithmetic<T>> = _, DisableIf<std::is_same<T, uint>> = _>
Metric<R, T>
euclidean() {
	return [](const T& a, const T& b) -> R {
		return abs_diff(a, b);
	};
}

// on uint's (discrete space, particularly useful for measuring distances
// between channel inputs), only consecutive elements are adjacent
//
template<typename R, typename T, EnableIf<std::is_same<T, uint>> = _>
Metric<R, T>
euclidean() {
	Metric<R, T> d = [](const T& a, const T& b) -> R {
		return abs_diff(a, b);
	};
	d.is_adjacent = [](const T& a, const T& b) -> bool {
		return abs_diff(a, b) == 1;
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
	d2.is_adjacent = d.is_adjacent;
	return d2;
}

template<typename R, typename T, EnableIf<is_Point<T>> = _, DisableIf<std::is_same<T, Point<uint>>> = _>
Metric<R, T>
euclidean() {
	return [](const T& a, const T& b) -> R {
		auto v1 = a.x - b.x,
			 v2 = a.y - b.y;
		return std::sqrt(v1*v1 + v2*v2);
	};
}

// Euclidean distance on discrete (uint) points. The only _non_-adjacent points
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
	d.is_adjacent = [](const T& a, const T& b) -> bool {
		uint v1 = abs_diff(a.x, b.x),
			 v2 = abs_diff(a.y, b.y);
		return !(v1 == 0 || v2 == 0 || v1 == v2) || (v1 <= 1 && v2 <= 1);
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

// Mahattan distance on discrete (uint) points. The only _adjacent_ points
// are those whose index differs by at most 1
//
template<typename R, typename T, EnableIf<std::is_same<T, Point<uint>>> = _>
Metric<R, T>
manhattan() {
	Metric<R, T> d = [](const T& a, const T& b) -> R {
		return abs_diff(a.x, b.x) + abs_diff(a.y, b.y);
	};
	d.is_adjacent = [](const T& a, const T& b) -> bool {
		return abs_diff(a.x, b.x) <= 1 && abs_diff(a.y, b.y) <= 1;
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
	d2.is_adjacent = [width, d](const uint& a, const uint& b) -> bool {
		return d.is_adjacent( Point<uint>(a%width, a/width), Point<uint>(b%width, b/width) );
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
// WISHLIST: faster method: jorlin.scripts.mit.edu/docs/publications/26-faster strongly polynomial.pdf
//
template<typename R, typename T, EnableIf<is_Prob<T>> = _>
Metric<R, T>
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

} // namespace metric

template<typename R, typename T>
Metric<R, T>
operator*(R coeff, const Metric<R, T>& d) {
    return metric::scale(d, coeff);
}

