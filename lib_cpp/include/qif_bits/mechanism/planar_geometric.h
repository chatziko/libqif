namespace mechanism::geo_ind {

using std::cerr;
using std::max;
using ::qif::geo::GridWalk;


double inverse_cumulative_gamma(double epsilon, double p);
int bound(int val, int min_val, int max_val);


// -- For sampling -----------------------------------------------

// sample from a planar geometric centered at (0,0) (if origin is different, just add it to the result)
//
template<typename eT = eT_def>
eT
_planar_geometric_coeff(eT cell_size, eT eps) {
	LargeSum<eT> sum;
	Point<eT> zero(0,0);
	auto d = eps * metric::euclidean<double, Point<eT>>();

	eT coeff_cur(0);
	for(Point<eT> p : GridWalk<eT>(cell_size)) {
		eT prob = exp<eT>(-d(p, zero));
		sum.add(prob);

		// stop if adding a new value does not change the result at all
		eT coeff_new(1 / sum.value());
		if(equal(coeff_new, coeff_cur))
			break;
		coeff_cur = coeff_new;
	}
	return coeff_cur;
}

template<typename eT = eT_def>
Point<eT>
planar_geometric_sample(eT cell_size, eT eps) {
	auto euclid = metric::euclidean<double, Point<eT>>();
	double coeff = _planar_geometric_coeff<eT>(cell_size, eps);

	eT p = rng::randu<eT>();
	LargeSum<eT> accu;
	uint cnt = 0;
	Point<eT> zero(0,0);

	for(Point<eT> z : GridWalk<eT>(cell_size)) {
		eT prob = coeff * exp<eT>(-eps * euclid(z, zero));
		accu.add(prob);

		if(accu.value() > p)
			return z;

		if(cnt++ == 1e8)
			throw std::runtime_error("infinite loop?");
	}
	return Point<eT>(0,0);		// unreachable, just avoid the warning
}

// efficient batch sampling
template<typename eT = eT_def>
std::vector<Point<eT>>
planar_geometric_sample(eT cell_size, eT eps, uint n) {
	auto euclid = metric::euclidean<double, Point<eT>>();
	double coeff = _planar_geometric_coeff<eT>(cell_size, eps);

	// we need n numbers uniformly sampled in [0,1]. Sort them and keep the indexes of the sorted list in 'order'
	Row<eT> ps(n);
	ps.randu();
	arma::uvec order = arma::sort_index(ps);

	// single grid walk for all samples
	GridWalk<eT> gw(cell_size);
	auto gw_it = gw.begin();

	LargeSum<eT> accu;
	accu.add(coeff);		// gw points at the first element (x itself), so accu should contain its probability (= coeff)!
	uint cnt = 0;
	Point<eT> zero(0,0);
	std::vector<Point<eT>> res(n);

	// sample n elements. 
	for(uint i = 0; i < n; i++) {
		// we need to visit elements in sorted order of p
		eT p = ps(order(i));

		while(less_than(accu.value(), p)) {
			++gw_it;	// first increment, to get the probability of the new point
			eT prob = coeff * exp<eT>(-eps * euclid(*gw_it, zero));
			accu.add(prob);
			if(cnt++ == 1e8)
				throw std::runtime_error("infinite loop?");
		}

		// place the result in the same position in res as p was in ps.
		res[order(i)] = *gw_it;
	}

	return res;
}



// -- For computing the channel matrix -----------------------------------------------

template<typename eT>
Mat<eT>
grid_summation(uint width, uint height, eT step, eT epsilon) {
	if(width % 2 == 0 || height % 2 == 0) throw std::runtime_error("width/height must by odd");
	if(width == 1 || height == 1) throw std::runtime_error("width/height must be > 1");

	// same trick changing step and epsilon. Not sure if it's really useful
	//
	double a_step = 1;
	epsilon *= step / a_step;								// we _multiply_ epsilon by the same factor used to _divide_ step. i.e. when a_step < step then epsilon is increased

	Mat<eT> m = arma::zeros<Mat<eT>>(height, width);		// rows = height, cols = width,  point (x,y) has index (cy+y, cx+x)

	int cx = (width-1)/2;									// corner x index
	int cy = (height-1)/2;									// corner y index
	const double grid_bound = (max(cx, cy) + 0.5) * a_step;	// boundary of the grid

	// we set the limit of integration to the distance in which 1-1e-6 of the
	// probability is included (and at least as big as the grid boundary)
	const int far_away = max(inverse_cumulative_gamma(epsilon, 1-1e-6), grid_bound) / a_step;

	if(debug) {
		cerr << "integration epsilon: " << epsilon << "\n";
		cerr << "far_away: " << far_away << "\n";
		cerr << "grid_bound: " << grid_bound << "\n";
	}

	int c = 0;
	for(int i = 0; i <= cx; i++) {
		for(int j = i; j <= cy; j++) {
			// summation from (i,j) = (end_i, end_j)
			int end_i = i == cx ? far_away : i;
			int end_j = j == cy ? far_away : j;

			// summation
			eT prob(0);

			for(int si = i; si <= end_i; si++)
				for(int sj = j; sj <= end_j; sj++)
					prob += std::exp( - epsilon * a_step * std::sqrt(si*si + sj*sj) );

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
	m /= arma::accu(m);

	if(debug)
		cerr << "summations " << c << "\n";

	return m;
}


template<typename eT>
Chan<eT>
planar_geometric_grid(uint width, uint height, eT step, eT epsilon) {

	int size = width * height;

	Chan<eT> C = arma::zeros<Mat<eT>>(size, size);

	int cx = width - 1;
	int cy = height - 1;

	Mat<eT> big_matrix = grid_summation(2*width-1, 2*height-1, step, epsilon);

	// we loop over the small grid, the current element is the chanel input
	//
	for(uint input_x = 0; input_x < width; input_x++) {
		for(uint input_y = 0; input_y < height; input_y++) {
			uint input = input_y * width + input_x;

			// we put this element in the center of the big grid
			// we loop over the big grid, and add this probability to
			// the closest point in the small grid
			//
			for(int x = -cx; x <= cx; x++) {
				uint output_x = bound(input_x + x, 0, width-1);

				for(int y = -cy; y <= cy; y++) {
					uint output_y = bound(input_y + y, 0, height-1);
					uint output = output_y * width + output_x;

					C(input, output) += big_matrix(cy+y, cx+x);
				}
			}
		}
	}

	return C;
}

} // namespace mechanism::geo_ind

