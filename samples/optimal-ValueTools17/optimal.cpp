#include "qif"

using std::cout;
using std::vector;
using std::string;

using namespace qif;


void compute_optimal(string method) {
	lp::Defaults::method = lp::method_t::simplex_dual;
	lp::Defaults::glp_msg_level = lp::msg_level_t::all;
	lp::Defaults::glp_presolve = true;

	double unit = 1;
	double eps = std::log(2) / (2*unit);
	uint width = 8,
		 n_inputs = width * width,
		 n_outputs = n_inputs;
	double R = 2.2 * unit;								// the R and rho for
	double rho = unit * std::sqrt(2)/2;					// Kostas' method

	auto d_grid = unit * metric::grid<double>(width);	// d_grid(i, j) = euclidean distance between cell indices i,j
	auto d_loss = d_grid;								// loss metric = euclidean

	Metric<double, uint> dx;							// privacy metric, depends on the method

	if(method == "direct") {
		dx = d_grid;									// the real deal

	} else {
		dx = metric::threshold_inf(d_grid, R);			// inf if above threshold R

		// Catuscia's crazy formula
		double R2 = R * R;
		double rho2 = rho * rho;
		double delta = (R - rho) / std::sqrt(
			R2 - 2*rho*R - 3*rho2 + 4*rho*rho2/R + rho2*rho2/R2
		);
		cout << "delta: " << delta << "\n";

		eps /= delta;									// update eps to compensate for the use of a bigger dx
	}

	prob pi = probab::uniform<double>(n_inputs);		// uniform prior

	chan opt = mechanism::optimal_utility(pi, n_outputs, eps * dx, d_loss);

	cout << "size: " << opt.n_rows << "\n";
	cout << "proper: " << channel::is_proper(opt) << "\n";
	cout << "util : " << l::post_entropy(d_loss, pi, opt) << "\n";	// expected d_loss between x and y
}

int main(int argc, char *argv[]) {
	vector<string> args(argv+1, argv + argc);


//	compute_optimal("direct");
	compute_optimal("kostas");
}


