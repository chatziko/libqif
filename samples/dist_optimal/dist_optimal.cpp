
#include <string>
#include <iostream>     // std::cout, std::fixed

using std::cout;

#include "qif"
using namespace qif;
using namespace qif::measure;


void print_mech(std::string name, chan& C, double eps, prob& pi, Metric<double, uint> d, Metric<double, uint> loss) {
	cout << "\n---- " << name << " ----\n";
//	cout << C;
	cout << "eps deviation: " << (mechanism::smallest_epsilon(C, d) - eps) << "\n";
	cout << "proper: " << channel::is_proper(C) << "\n";
	cout << "util with remap : " << l_uncert::posterior(loss, pi, C) << "\n";
	cout << "util without remap : " << utility::expected_distance(loss, pi, C) << "\n";
}

// read prior from Ehab's files.
// Start from bottom-left corner and go first right and then up
//
prob load_prior(std::string file) {
	mat M;
	M.load(file);

	prob pi(M.n_elem);
	uint c = 0;
	for(int i = M.n_rows - 1; i >= 0; i--)
		for(uint j = 0; j < M.n_cols; j++)
			pi(c++) = M(i,j);
	std::cout << "prior sum: " << arma::accu(pi) << "\n";
	return pi;
}

int main() {
	double
		alpha = 1.2,
		eps = std::log(alpha),
		cell_size = 1;
	uint
		width = 7,
		height = width,
		n_inputs = width * height,
		n_outputs = n_inputs;

	std::cout << "eps: ln(" << alpha << ")\n";
	std::cout << "size: " << width << "x" << height << "\n";

	auto geo_d = cell_size * metric::euclidean<double, Point<uint>>();
//	auto geo_d = cell_size * metric::manhattan<double, Point<uint>>();

	auto dx = metric::grid(width, geo_d);		// distance on grid from distance on points
	auto loss = dx;

	prob pi = probab::uniform<double>(n_inputs);
//	prob pi = load_prior("/home/vagabond/Desktop/ehab_priors/p1.pr");


	chan opt = mechanism::optimal_exp_loss::under_d_privacy(pi, n_outputs, eps * dx, loss);
	print_mech("optimal", opt, eps, pi, dx, loss);

	chan dist_opt = mechanism::optimal_exp_loss::under_d_privacy(pi, n_outputs, eps * dx, loss, "dist");
	print_mech("dist-optimal", dist_opt, eps, pi, dx, loss);

	chan dist_opt_strict = mechanism::optimal_exp_loss::under_d_privacy(pi, n_outputs, eps * dx, loss, "dist_strict");
	print_mech("dist-optimal-strict", dist_opt_strict, eps, pi, dx, loss);

	chan laplace = mechanism::planar_laplace_grid<double>(width, height, cell_size, eps);
	print_mech("laplace", laplace, eps, pi, dx, loss);
}
