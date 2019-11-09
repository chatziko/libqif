#include "qif"
#include <math.h>       /* log */
using namespace qif;
using namespace qif::measure;

using std::cout;
using std::string;

const uint grid_size = 100;		// 100x100 grid
const double cell_width = 0.5;	// each cell is 0.5x0.5 KM


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
	/*laplace = mechanism::geo::planar_laplace_grid<double>(3,3,1, 0.5);
	cout << laplace << "\n";
	laplace.save("temp/laplace-"+std::to_string(0.5)+".mat");*/
	/*laplace = mechanism::geo::planar_laplace_grid<double>(grid_size, grid_size, cell_width, eps);		
	cout << "ok" << "\n";
	laplace.save("temp/laplace-"+std::to_string(eps)+".mat");*/
	
	if(!laplace.load("temp_foursquare/laplace-"+std::to_string(eps)+".mat")) {

		laplace = mechanism::geo::planar_laplace_grid<double>(grid_size, grid_size, cell_width, eps);
		
		cout << "ok created laplace" << "\n";
		laplace.save("temp_foursquare/laplace-"+std::to_string(eps)+".mat");
	}
	uint n = laplace.n_rows;
	//cout << typeid(laplace).name() << "\n";

	cout << laplace.n_cols << "\n";

	// load workes with csv also 
	mat priors;
	priors.load("preproc_foursquare/temp/"+areadataset+"-user.csv");
	cout << priors.n_rows<< "\n";
	prob prior_global;
	prior_global.load("preproc_foursquare/temp/"+areadataset+"-global.csv");

	cout << "Loading priors done!" << "\n";

	mat L;
	L.resize(n, n);
	create_l(area, dataset, priv_metric, L);

	arma::ucolvec strategy;

	if (!strategy.load("temp_foursquare/strategy-"+std::to_string(eps)+".mat"))
	{
		cout << "Creating strategy" << "\n";

		auto strategy = l_uncert::strategy(L, prior_global, laplace);

		cout << "Strategy done!" << "\n";
		strategy.save("temp_foursquare/strategy-"+std::to_string(eps)+".mat");
		
	}

	cout << "strategy loaded" << "\n";

	std::ofstream myfile;
	myfile.open("generated_data_foursquare/laplace-"+std::to_string(eps) +"-"+ areadataset + "-" + priv_metric);

	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);

		double privacy = 0.0;
		for(uint i = 0; i < n; i++){
			
			privacy += pi(i) * x_privacy(laplace, L, strategy, i);
		}
		std::cout << "iter done ----------------------------------------" << ": \n";
		std::cout << i << ": of "<< priors.n_rows<< "\n";
		std::cout << privacy <<"\n";
		myfile << privacy << "\n";	
		
	}
	myfile.close();
}


void compute_laplace_privacy_seml(string area, string dataset, string priv_metric, string seml,double eps) {
	string areadataset = area + "-" + dataset;

	chan laplace;
	/*laplace = mechanism::geo::planar_laplace_grid<double>(3,3,1, 0.5);
	cout << laplace << "\n";
	laplace.save("temp/laplace-"+std::to_string(0.5)+".mat");*/
	/*laplace = mechanism::geo::planar_laplace_grid<double>(grid_size, grid_size, cell_width, eps);		
	cout << "ok" << "\n";
	laplace.save("temp/laplace-"+std::to_string(eps)+".mat");*/
	
	if(!laplace.load("temp_foursquare/laplace-"+std::to_string(eps)+".mat")) {

		laplace = mechanism::geo::planar_laplace_grid<double>(grid_size, grid_size, cell_width, eps);
		
		cout << "ok created laplace" << "\n";
		laplace.save("temp_foursquare/laplace-"+std::to_string(eps)+".mat");
	}
	uint n = laplace.n_rows;
	//cout << typeid(laplace).name() << "\n";

	cout << laplace.n_cols << "\n";

	// load workes with csv also 
	mat priors;
	priors.load("preproc_foursquare/temp/"+areadataset+"-user-"+seml+".csv");
	cout << priors.n_rows<< "\n";
	prob prior_global;
	prior_global.load("preproc_foursquare/temp/"+areadataset+"-global-"+seml+".csv");

	cout << "Loading priors done!" << "\n";

	mat L;
	L.resize(n, n);
	create_l(area, dataset, priv_metric, L);

	arma::ucolvec strategy;

	if (!strategy.load("temp_foursquare/strategy-"+std::to_string(eps)+"-"+seml+".mat"))
	{
		cout << "Creating strategy" << "\n";

		auto strategy = l_uncert::strategy(L, prior_global, laplace);

		cout << "Strategy done!" << "\n";
		strategy.save("temp_foursquare/strategy-"+std::to_string(eps)+"-"+seml+".mat");
		
	}

	cout << "strategy loaded" << "\n";

	std::ofstream myfile;
	myfile.open("generated_data_foursquare/laplace-"+seml+"-"+std::to_string(eps) +"-"+ areadataset + "-" + priv_metric);

	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);

		double privacy = 0.0;
		for(uint i = 0; i < n; i++){
			
			privacy += pi(i) * x_privacy(laplace, L, strategy, i);
		}
		std::cout << "iter done ----------------------------------------" << ": \n";
		std::cout << i << ": of "<< priors.n_rows<< "\n";
		std::cout << privacy <<"\n";
		myfile << privacy << "\n";	
		
	}
	myfile.close();
}

int main() {
	arma::arma_rng::set_seed_random();

	double eps = log(2)/0.5;
	cout << "using eps: " << eps << "\n";
	string seml = "Medical Center";


	for(string area        : { "nyc.csv"})
		for(string dataset     : { "foursquare"})
			for(string priv_metric : {"euclidean"}) {
				
				//compute_laplace_privacy(area, dataset, priv_metric, eps);
				compute_laplace_privacy_seml(area, dataset, priv_metric, seml,eps);
				

			}

}
