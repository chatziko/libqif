#include <qif>
#include <string>

using namespace qif;
using namespace qif::measure;
using namespace std;

int usage() {
	cout
		<< "Usage: geo-optimal <goal> <width> <height> <cell_size> <constraint> <hard_max_loss> <prior> <solver>\n\n"
		<< "goal:          one of\n"
		<< "                 min_loss_given_min_bayesrisk\n"
		<< "                 min_loss_given_min_georisk\n"
		<< "                 max_bayesrisk_given_max_loss\n"
		<< "                 max_georisk_given_max_loss\n"
		<< "width:         width of grid\n"
		<< "height:        height of grid\n"
		<< "cell_size:     length of each cell\n"
		<< "constraint:    the bayesrisk, georisk or loss constraint, depending on the goal\n"
		<< "hard_max_loss: C_xy is forced to 0 when loss(x,y) > hard_max_loss (can greatly reduce the problem size)\n"
		<< "prior:         file with a single line containing the prior | uniform | random\n"
		<< "solver:        CLP | GLOP | GLPK\n";
	return -1;
}

int main(int argc, char* argv[]) {
	qif::rng::set_seed_random();		// RNG initialization

	if(argc != 9) return usage();

	string goal = argv[1];
	uint width = std::stoi(argv[2]);
	uint height = std::stoi(argv[3]);
	uint n = width * height;
	double cell_size = std::stod(argv[4]);

	auto euclid = metric::euclidean<double, point>();
	auto loss = cell_size * metric::compose(euclid, geo::cell_to_point(width));

	double constraint = std::stod(argv[5]);
	double hard_max_loss = std::stod(argv[6]);
	string pi_file = argv[7];
	string solver = argv[8];

	prob pi;
	if(pi_file == "uniform")
		pi = probab::uniform(n);
	else if(pi_file == "random")
		pi = probab::randu(n);
	else {
		if(!pi.load(pi_file)) {
			cerr << "cannot open " << pi_file << "\n";
			return 1;
		}
		if(pi.n_elem != n) {
			cerr << "prior should have size " << n << "\n";
			return 1;
		}
		if(!probab::is_proper(pi))
			cerr << "warning: prior is not a proper distribution\n";
	}

	lp::Defaults::msg_level = lp::MsgLevel::ALL;
	lp::Defaults::solver = solver;

	cout << "computing optimal mechanism\n";

	chan C;

	if(goal == "min_loss_given_min_bayesrisk")
		C = mechanism::bayes_risk::min_loss_given_min_risk(pi, n, constraint, loss, hard_max_loss);
	else if(goal == "min_loss_given_min_georisk")
		C = mechanism::l_risk::min_loss_given_min_risk(pi, n, n, constraint, loss, loss, hard_max_loss);
	else if(goal == "max_bayesrisk_given_max_loss")
		C = mechanism::bayes_risk::max_risk_given_max_loss(pi, n, constraint, loss, hard_max_loss);
	else if(goal == "max_georisk_given_max_loss")
		C = mechanism::l_risk::max_risk_given_max_loss(pi, n, n, constraint, loss, loss, hard_max_loss);
	else
		return usage();

	cout << "done\n";

	if(C.is_empty()) {
		cout << "no solution for given constraints\n";
		return -1;
	}

	C.save("C", arma::raw_ascii);

	cout << "Channel size: " << C.n_rows << "x" << C.n_cols << "\n";
	cout << "BayesRisk: " << bayes_risk::posterior(pi, C) << "\n";
	cout << "GeoRisk: " << l_risk::posterior(loss, pi, C) << "\n";
	cout << "Exp Util Loss: " << utility::expected_distance(loss, pi, C) << "\n";
}
