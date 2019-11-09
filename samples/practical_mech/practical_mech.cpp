
#include <string>
#include <iostream>     // std::cout, std::fixed
#include <map>

using std::cout;
using std::string;

#include "qif"
using namespace qif;
using namespace qif::measure;

typedef std::map<string, prob> probdict;



// scales a geo mechanism
// C is a mechanism working on a grid of width_c
// The number of cells on each axis is split in fact, so the
// total number of cells is multiplied by fact^2
//
chan scale_geo_mechanism(const chan& C, uint width_c, uint fact) {
	// The cell of index 0 is (0,0) (bottom left)
	// index i => cell (i%width, i/width)
	// cell (x,y) => index x*width + y
	//
	// fine   coord x_f => coarse coord x_f/fact
	// coarse coord x_c => fine   coord x_c * fact + fact/2
	//
	assert(C.n_rows == C.n_cols);

	uint height_c = C.n_rows / width_c;

	uint width_f  = width_c  * fact;
	uint height_f = height_c * fact;

	chan D(width_f * height_f, width_f * height_f);
	D.fill(0);

	// loop over input i_f=(x_f,y_f) and output ii_c(xx_c,xx_c)
	for(uint x_f = 0; x_f < width_f; x_f++) {
		for(uint y_f = 0; y_f < height_f; y_f++) {

			uint i_f = x_f        * width_f + y_f;
			uint i_c = x_f/fact * width_c + y_f/fact;

			for(uint xx_c = 0; xx_c < width_c; xx_c++) {
				for(uint yy_c = 0; yy_c < height_c; yy_c++) {

					uint ii_c = xx_c               * width_c + yy_c;
					uint ii_f = (xx_c*fact+fact/2) * width_f + (yy_c*fact+fact/2);

					//std::cout <<  " i_f: " <<  i_f <<  " i_c: " <<  i_c << " ii_f: " << ii_f << " ii_f: " << ii_c << "\n";

					D(i_f,ii_f) = C(i_c,ii_c);
				}
			}
		}
	}

	return D;
}


void print_mech(string name, const chan& C, const prob& pi_global, mat& priors, const mat& Loss) {
	if(C.n_cols == 0) {
		cout << "# " << name << " not available\n\n\n";
		return;
	}

	cout << name << " no remap";
	arma::vec errors(priors.n_rows);
	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);
//		cout << ", " << utility::expected_distance(Loss, pi, C);
		errors(i) = utility::expected_distance(Loss, pi, C);
	}
	cout << ": " << arma::mean(errors) << ", " << arma::median(errors);
	cout << "\n";

	chan R = channel::deterministic<double>(l_risk::strategy(Loss, pi_global, C), Loss.n_rows);
	chan CR = C * R;

	cout << name << " with remap";
	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);
//		cout << ", " << utility::expected_distance(Loss, pi, CR);
		errors(i) = utility::expected_distance(Loss, pi, CR);
	}
	cout << ": " << arma::mean(errors) << ", " << arma::median(errors);
	cout << "\n";
}

void print_mech2(string name, const chan& C, const probdict& globals, mat& priors, const mat& Loss) {
	if(C.n_cols == 0) {
		cout << "# " << name << " not available\n\n\n";
		return;
	}

	for(auto g : globals) {
		cout << "# " << name << "-" << g.first << "\n";

		chan M = g.second.is_empty()
			? C
			: C * channel::deterministic<double>(l_risk::strategy(Loss, g.second, C), Loss.n_rows);

		for(uint i = 0; i < priors.n_rows; i++) {
			prob pi = priors.row(i);
			cout << utility::expected_distance(Loss, pi, M) << "\n";
		}
		cout << "\n\n";
	}
}

void compare_noisy_remaps() {
	double
		cell_size = 0.2;
	uint
		factor = 10,
		width_c = 6,
		height_c = 14,
		width = width_c * factor,
		height = height_c * factor,
		n_inputs = width * height;
//		n_inputs_c = width_c * height_c;
//		n_outputs_c = n_inputs_c;

	double alphas[] = { 1.4, 1.7, 2.0, 2.3 };

//	string path = "/home/vagabond/Downloads/datasets/practical-mech-cache";
	string path = "/home/ubuntu/dataset";

	mat priors;
	priors.load(path + "/SF_fine_users");

	probdict globals;
	globals["noremap"];									// empty global to include unremaped mechanism
	globals["exact"].load(path + "/SF_fine_global");
	globals["noisy5"].load(path + "/SF_fine_global.noisy5");
	globals["noisy50"].load(path + "/SF_fine_global.noisy50");

//	if(!(priors.n_cols == n_inputs && globals["exact"].n_cols == n_inputs))
//		throw new std::runtime_error("invalid size");
//	assert(pi_global_c.n_cols == n_inputs_c);

	auto geo_d   = cell_size          * metric::euclidean<double, Point<uint>>();
	auto geo_d_c = cell_size * factor * metric::euclidean<double, Point<uint>>();

	auto dx   = metric::grid(width,   geo_d  );	// distance on grid from distance on points
	auto dx_c = metric::grid(width_c, geo_d_c);	// distance on grid from distance on points
	auto loss   = dx;
	auto loss_c = dx_c;
	mat Loss = l_risk::metric_to_mat(loss, n_inputs);

	for(double alpha : alphas) {
		double eps = std::log(alpha) / 0.1;

//		std::cout << "\n\n----------------------------\neps: ln("
//			<< alpha << ")/0.1 = " << eps << "\n";

//		chan opt_c = mechanism::d_priv::min_loss_given_d<double>(pi_global_c, n_outputs_c, eps * dx_c, loss_c);
//		chan opt = scale_geo_mechanism(opt_c, width_c, factor);
//		print_mech("optimal", opt, pi_global, priors, Loss);
//		opt.clear();
//		opt_c.clear();

//		chan tight = mechanism::d_priv::tight_constraints(n_inputs, eps * dx);
//		print_mech2(std::to_string(alpha)+"-tight", tight, pi_global, pi_global_n, priors, Loss);
//		tight.clear();

		chan laplace = mechanism::geo::planar_laplace_grid<double>(width, height, cell_size, eps);
		print_mech2(std::to_string(alpha)+"-laplace", laplace, globals, priors, Loss);
		laplace.clear();

		chan geom = mechanism::geo::planar_geometric_grid<double>(width, height, cell_size, eps);
		print_mech2(std::to_string(alpha)+"-geometric", geom, globals, priors, Loss);
		geom.clear();
	}
}




int coarse() {
	double
		cell_size = 2;
	uint
		width = 6,
		height = 14,
		n_inputs = width * height,
		n_outputs = n_inputs;

	double alphas[] = { 1.4, 1.7, 2, 2.3 };

	prob pi_global;
	mat priors;
	pi_global.load("/home/vagabond/Downloads/datasets/practical-mech-cache/coarse-global-prior");
	priors.load("/home/vagabond/Downloads/datasets/practical-mech-cache/coarse-user-priors");

	assert(pi_global.n_cols == n_inputs);
	assert(priors.n_cols == n_inputs);

	std::cout << "size: " << width << "x" << height << "\n";

	auto geo_d   = cell_size          * metric::euclidean<double, Point<uint>>();

	auto dx   = metric::grid(width,   geo_d  );	// distance on grid from distance on points
	auto loss   = dx;
	mat Loss = l_risk::metric_to_mat(loss, n_inputs);

	for(double alpha : alphas) {
		double eps = std::log(alpha) / 0.1;

		std::cout << "\n\n----------------------------\neps: ln("
			<< alpha << ")/0.1 = " << eps << "\n";

		chan opt = mechanism::d_priv::min_loss_given_d<double>(pi_global, n_outputs, eps * dx, loss);
		if(opt.is_empty())
			continue;
		std::cout << alpha << ", c[0,0]: " << opt(0,0) << "\n";
		continue;
		print_mech("optimal", opt, pi_global, priors, Loss);
		opt.clear();

		chan tight = mechanism::d_priv::tight_constraints(n_inputs, eps * dx);
		print_mech("tight", tight, pi_global, priors, Loss);
		tight.clear();

		chan laplace = mechanism::geo::planar_laplace_grid<double>(width, height, cell_size, eps);
		cout << "laplace smallest eps " << measure::d_priv::smallest_epsilon(laplace, dx) << "\n";
		print_mech("laplace", laplace, pi_global, priors, Loss);
		std::cout << laplace;
		laplace.clear();

		chan geom = mechanism::geo::planar_geometric_grid<double>(width, height, cell_size, eps);
		print_mech("geometric", geom, pi_global, priors, Loss);
		cout << "geom smallest eps " << measure::d_priv::smallest_epsilon(geom, dx) << "\n";
		std::cout << geom;
		geom.clear();
	}
	return 0;
}

int main() {
//	compare_noisy_remaps();
	return coarse();
}
