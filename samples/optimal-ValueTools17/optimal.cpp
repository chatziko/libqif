#include "qif"

using std::cout;
using std::vector;
using std::string;

using namespace qif;
using namespace qif::measure;


void compute_optimal(string method) {
	if(method != "direct" && method != "kostas" && method != "geometric")
		throw std::runtime_error("invalid method: " + method);

	lp::Defaults::method = lp::Method::SIMPLEX_DUAL;
	lp::Defaults::msg_level = lp::MsgLevel::ALL;
	lp::Defaults::presolve = true;

	double unit = 1;
	double eps = std::log(2) / (2*unit);
	uint width = 8,
		 n_inputs = width * width,
		 n_outputs = n_inputs;
	double R = 2.2 * unit;								// the R and rho for
	double rho = unit * std::sqrt(2)/2;					// Kostas' method

	auto d_grid = unit * metric::grid<double>(width);	// d_grid(i, j) = euclidean distance between cell indices i,j
	auto d_loss = d_grid;								// loss metric = euclidean
	auto dx = d_grid;									// privacy metric

	// this disables the removal of constraints that are reduntant due to transitivity
	dx.chainable = [](const uint&, const uint&) -> bool { return false; };

	if(method == "kostas") {
		// use more relexed dx and update eps
		//
		dx = metric::threshold_inf(d_grid, R);			// inf if above threshold R

		// Catuscia's crazy formula
		double c = R / rho;
		double c2 = c * c;
		double delta = (c - 1) / std::sqrt(
			c2 - 2*c - 3 + 4/c + 1/c2
		);
		cout << "delta: " << delta << "\n";

		eps /= delta;									// update eps to compensate for the use of a bigger dx
	}

	prob pi = probab::uniform<double>(n_inputs);		// uniform prior
	chan C = method == "geometric"
		? mechanism::geo::planar_geometric_grid<double>(width, width, unit, eps)
		: mechanism::optimal_exp_loss::under_d_privacy(pi, n_outputs, eps * dx, d_loss);

	cout << "size: " << C.n_rows << "\n";
	cout << "proper: " << channel::is_proper(C) << "\n";
	cout << "util : " << l_uncert::posterior(d_loss, pi, C) << "\n";	// expected d_loss between x and y

	// apply optimal remap (probably does nothing for uniform priors)
	mat Loss = l_uncert::metric_to_mat(d_loss, n_inputs);
	chan Remap = channel::deterministic<double>(l_uncert::strategy(Loss, pi, C), Loss.n_rows);	// compute remap
	cout << "util after remap : " << l_uncert::posterior(d_loss, pi, (chan)(C*Remap)) << "\n";
}

int main(int argc, char *argv[]) {
	vector<string> args(argv+1, argv + argc);


//	compute_optimal("direct");
//	compute_optimal("geometric");
	compute_optimal("kostas");
}


