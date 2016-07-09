#include "qif"
#include <math.h>       /* log */
using namespace qif;

using std::cout;
using std::string;

const uint grid_size = 100;		// 100x100 grid
const double cell_width = 100;	// each cell is 100x100meters


double utility_to_epsilon(string dataset, string util_metric) {

	/*prob prior_global;
	prior_global.load("marco/priors/paris-center-"+dataset+"-global.dat");

	arma::rowvec util;
	util.load("marco/utility/paris-center-"+util_metric+".dat");
	double utility = arma::cdot(prior_global, util);

	double eps = util_metric == "euclidean"
		? 2.0 / utility
		: mechanism::inverse_cumulative_gamma(stod(util_metric.replace(0, 7, "")), utility);

	cout << "utility for " << dataset << "/" << util_metric << ": " << utility << ", eps: " << eps << "\n";*/
	double eps = log(4);

	return eps;
}


double x_privacy(const chan& C, mat& L, const arma::ucolvec strategy, uint i) {
	arma::colvec c = L.col(i);
	return arma::cdot(C.row(i), c(strategy));
}


void create_l(string area, string dataset, string priv_metric, mat& L) {
	string areadataset = area + "-" + dataset;

	if(priv_metric == "binary") {

		mat pois = arma::ones(1, L.n_rows);
		pois.load("marco/poi/"+areadataset+".dat");
		for(auto& e : pois)
			if(!equal(e, 0.0))
				e = 1-1/e;
			else
				e = 1;

		L.ones();
		L.diag() = pois;

	} else {

		Metric<double, uint> euclid = cell_width * metric::grid<double>(grid_size);
		for(uint i = 0; i < L.n_rows; i++)
			for(uint j = i; j < L.n_cols; j++)
				L.at(i, j) = L.at(j, i) = euclid(i, j);

	}
}


void compute_laplace_privacy(string area, string dataset, string priv_metric, double eps) {
	string areadataset = area + "-" + dataset;

	chan laplace;
	if(!laplace.load("temp/laplace-"+std::to_string(eps)+".bin")) {

		laplace = mechanism::planar_laplace_grid<double>(grid_size, grid_size, cell_width, eps);
		laplace.save("temp/laplace-"+std::to_string(eps)+".bin");
	}
	uint n = laplace.n_rows;

	// load workes with csv also 
	mat priors;
	priors.load(areadataset+"-user.csv");

	prob prior_global;
	prior_global.load(areadataset+"-global.csv");

	mat L;
	L.resize(n, n);
	create_l(area, dataset, priv_metric, L);

	auto strategy = l::strategy(L, prior_global, laplace);

	std::ofstream myfile;
	myfile.open("generated_data/laplace-" + areadataset + "-" + priv_metric);

	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);

		double privacy = 0.0;
		for(uint i = 0; i < n; i++)
			privacy += pi(i) * x_privacy(laplace, L, strategy, i);
//		cout << arma::accu(pi > 0) << ": ";
		myfile << privacy << "\n";
	}
	myfile.close();
}





int main() {
	arma::arma_rng::set_seed_random();


	double eps_sum = 0.0;
	int n = 0;

	for(string dataset     : { "gowalla"})
		for(string util_metric : { "euclidean"}) {
			eps_sum += utility_to_epsilon(dataset, util_metric);
			n++;
		}

	double eps = eps_sum / n;
	cout << "using eps: " << eps << "\n";


	for(string area        : { "paris"})
		for(string dataset     : { "gowalla"})
			for(string priv_metric : {"euclidean"}) {
				//compute_elastic_privacy(area, dataset, priv_metric);
				compute_laplace_privacy(area, dataset, priv_metric, eps);
			}

}
