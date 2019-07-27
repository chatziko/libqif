#include <qif>
#include <string>


using namespace qif;
using namespace std;


// args: <width> <cell_size> <min_bayes_risk>

int main(int argc, char* argv[]) {
	qif::rng::set_seed_random();		// RNG initialization

	if(argc != 6) {
		cout << "args should be: <width> <height> <cell_size> <min_bayes_risk> <max_loss>\n";
		return -1;
	}

	uint width = std::stoi(argv[1]);
	uint height = std::stoi(argv[2]);
	uint n = width * height;
	double cell_size = std::stod(argv[3]);
	prob pi = probab::uniform(n);

	auto euclid = metric::euclidean<double, point>();
	auto loss = cell_size * metric::compose(euclid, geo::cell_to_point(width));
	double max_vuln = 1 - std::stod(argv[4]);
	double max_loss = std::stod(argv[5]);

	lp::Defaults::msg_level = lp::MsgLevel::ALL;

	cout << "computing optimal mechanism\n";

	chan C = mechanism::optimal_exp_loss::given_max_bayes_vulnerability(pi, n, max_vuln, loss, max_loss);

	cout << "done\n";

	if(C.is_empty()) {
		cout << "no solution for given constraints\n";
		return -1;
	}

	C.save("C", arma::raw_ascii);

	cout << "Channel size: " << C.n_rows << "x" << C.n_cols << "\n";
	cout << "Bayes Risk: " << bayes_risk::posterior(pi, C) << "\n";
	cout << "Exp Loss: " << utility::expected_distance(loss, pi, C) << "\n";
}
