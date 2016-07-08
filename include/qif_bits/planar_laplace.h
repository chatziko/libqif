namespace mechanism {

using std::cerr;
using std::max;

const int integration_calls = 50000;		// how many times to call the function for each integration
const bool debug = false;

double planar_laplace_pdf(double *x, size_t, void *);
double inverse_cumulative_gamma(double epsilon, double p);
int bound(int val, int min_val, int max_val);
double integrate_laplace(double epsilon, const arma::vec& a, const arma::vec& b, int calls = integration_calls);


template<typename eT>
Mat<eT>
grid_integration(uint width, uint height, eT step, eT epsilon) {
	if(width % 2 == 0 || height % 2 == 0) throw std::runtime_error("width/height must by odd");
	if(width == 1 || height == 1) throw std::runtime_error("width/height must be > 1");

	// integration is in theory more accurate when the area of integration is small.
	// In our case, we can simply scale down distances by adjusting epsilon
	// So we use a_step instead of the requested one, while properly adjusting
	// epsilon.
	// (Kostas: in practice I didn't see any difference, but still it seems a good idea to use a small step)
	//
	double a_step = 1;
	epsilon *= step / a_step;								// we _multiply_ epsilon by the same factor used to _divide_ step. i.e. when a_step < step then epsilon is increased

	Mat<eT> m = arma::zeros<Mat<eT>>(height, width);		// rows = height, cols = width,  point (x,y) has index (cy+y, cx+x)

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
			arma::vec a = {
				(i - 0.5) * a_step,
				(j - 0.5) * a_step
			};
			arma::vec b = {
				i == cx ? far_away : a(0) + a_step,
				j == cy ? far_away : a(1) + a_step,
			};

			// integrate. we use more calls when we go outside the grid boundary
			bool out_of_bound = !less_than_or_eq(b(0), grid_bound) || !less_than_or_eq(b(1), grid_bound);
			int calls = integration_calls * (out_of_bound ? 10 : 1);

			eT prob = integrate_laplace(epsilon, a, b, calls);
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

template<typename eT>
Chan<eT>
planar_laplace_grid(uint width, uint height, eT step, eT epsilon) {

	int size = width * height;

	Chan<eT> C = arma::zeros<Mat<eT>>(size, size);

	int cx = width - 1;
	int cy = height - 1;

	Mat<eT> big_matrix = grid_integration(2*width-1, 2*height-1, step, epsilon);

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

} // namespace mechanism

