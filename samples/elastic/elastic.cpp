#include "qif"

using namespace qif;
using namespace qif::measure;

using std::cout;
using std::string;

const uint grid_size = 100;		// 100x100 grid
const double cell_width = 100;	// each cell is 100x100meters

double utility_to_epsilon(string dataset, string util_metric) {

	prob prior_global;
	prior_global.load("marco/priors/paris-center-"+dataset+"-global.dat");

	arma::rowvec util;
	util.load("marco/utility/paris-center-"+util_metric+".dat");
	double utility = arma::cdot(prior_global, util);

	double eps = util_metric == "euclidean"
		? 2.0 / utility
		: mechanism::geo::inverse_cumulative_gamma(stod(util_metric.replace(0, 7, "")), utility);

	cout << "utility for " << dataset << "/" << util_metric << ": " << utility << ", eps: " << eps << "\n";

	return eps;
}


double x_privacy(const chan& C, mat& L, const arma::ucolvec strategy, uint i) {
	arma::colvec c = L.col(i);
	return arma::cdot(C.row(i), c(strategy));
}


//void output_x_privacy(const MinEntropy<double>& me, const mat& pois, const prob& pi) {

//	auto str = me.strategy(pi % pois);

//	for(uint i = 0; i < pi.n_elem; i++) {
//		if(pi(i) > 0)
//			cout << 1 - pois(i) * x_privacy(me.C, str, i) << "\n";
//	}
//}


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


void compute_elastic_privacy(string area, string dataset, string priv_metric) {
	string areadataset = area + "-" + dataset;

	chan elastic;
	if(!elastic.load("temp/elastic-"+area+".bin")) {

		mat dist;
		if(!dist.load("temp/metric-"+area+".bin")) {
			dist.load("marco/metric/"+area+".dat");
			dist.save("temp/metric-"+area+".bin");
		}
		Metric<double, uint> d = metric::from_distance_matrix(dist);


		elastic = mechanism::d_priv::exponential(dist.n_rows, d);
		elastic.save("temp/elastic-"+area+".bin");

		dist.reset();
	}
	uint n = elastic.n_rows;


	mat priors;
	priors.load("marco/priors/"+areadataset+"-user.dat");

	prob prior_global;
	prior_global.load("marco/priors/"+areadataset+"-global.dat");

	mat L;
	L.resize(n, n);
	create_l(area, dataset, priv_metric, L);

	auto strategy = l_uncert::strategy(L, prior_global, elastic);

	std::ofstream myfile;
	myfile.open("generated_data/elastic-" + areadataset + "-" + priv_metric);

	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);

		double privacy = 0.0;
		for(uint i = 0; i < n; i++)
			privacy += pi(i) * x_privacy(elastic, L, strategy, i);
//		cout << arma::accu(pi > 0) << ": ";
		myfile << privacy << "\n";
	}
	myfile.close();
}


void compute_laplace_privacy(string area, string dataset, string priv_metric, double eps) {
	string areadataset = area + "-" + dataset;

	chan laplace;
	if(!laplace.load("temp/laplace-"+std::to_string(eps)+".bin")) {

		laplace = mechanism::geo::planar_laplace_grid<double>(grid_size, grid_size, cell_width, eps);
		laplace.save("temp/laplace-"+std::to_string(eps)+".bin");
	}
	uint n = laplace.n_rows;


	mat priors;
	priors.load("marco/priors/"+areadataset+"-user.dat");

	prob prior_global;
	prior_global.load("marco/priors/"+areadataset+"-global.dat");

	mat L;
	L.resize(n, n);
	create_l(area, dataset, priv_metric, L);

	auto strategy = l_uncert::strategy(L, prior_global, laplace);

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

	for(string dataset     : { "gowalla",   "brightkite"  })
	for(string util_metric : { "euclidean", "binrad-1200", "binrad-1500", "binrad-1800" }) {
		eps_sum += utility_to_epsilon(dataset, util_metric);
		n++;
	}

	double eps = eps_sum / n;
	cout << "using eps: " << eps << "\n";


	for(string area        : { "paris-center", "paris-nanterre" })
	for(string dataset     : { "gowalla",      "brightkite"     })
	for(string priv_metric : { "binary",       "euclidean"      }) {
		compute_elastic_privacy(area, dataset, priv_metric);
		compute_laplace_privacy(area, dataset, priv_metric, eps);
	}





}



	// normalize priors, keep top X locations
	//
//	uint keep = 10;
//	for(uint i = 0; i < priors.n_rows; i++) {
//		prob pi = priors.row(i);
//		arma::umat remove = arma::sort_index(pi).cols(0, n - 1 - keep);
//		pi(remove) *= 0;
//		pi /= arma::accu(pi);
//		priors.row(i) = pi;
//	}
//	cout << priors.row(0) << "\n";
//	cout << arma::accu(priors.row(0) > 0) << "\n";
//	return 0;
//
//
//
//
//	cout << "user, prior privacy, elastic privacy, elastic utility\n";
//	for(uint i = 0; i < priors.n_rows; i++) {
//		prob pi = priors.row(i);
//		if(!is_proper<prob>(pi, 1e-1)) {
//			cout << i << " skip, " << arma::accu(pi) << "\n";
//			continue;
//		}


//		cout << i << ", ";
//		cout << (1 - me.vulnerability(pi % pois)) << ", ";

//		cout << (1 - me.cond_vulnerability(pi % pois)) << ", ";
//		cout << utility::expected_distance<double>(Euclid, pi, elastic) << "\n";
//	}
