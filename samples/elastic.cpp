#include "types.h"
#include "aux.h"
#include "Metric.h"
#include "GLeakage.h"
#include "MinEntropy.h"
#include "Mechanism.h"
#include "Chan.h"
#include "PlanarLaplace.h"
#include "utility.h"



using std::cout;
using std::string;


#include <cstdlib>
#include <stdexcept>
#include <execinfo.h>
#include <iostream>
#include <fstream>
void my_terminate() {
    static bool tried_throw = false;

    try {
        // try once to re-throw currently active exception
        if (!tried_throw++) throw;
    }
    catch (const std::exception &e) {
        std::cerr << "\n" << __FUNCTION__ << " caught unhandled exception. what(): "
                  << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << __FUNCTION__ << " caught unknown/unhandled exception." 
                  << std::endl;
    }

    void * array[50];
    int size = backtrace(array, 50);    
    std::cerr << __FUNCTION__ << " backtrace returned " 
              << size << " frames\n\n";
    char ** messages = backtrace_symbols(array, size);

    for (int i = 0; i < size && messages != NULL; ++i) {
        std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
    }
    std::cerr << std::endl;
    free(messages);
    abort();
}


void compute_avg_error(prob pi, string area, string dataset) {

	arma::rowvec error;
	error.load("marco/error/"+area+".dat");

	cout << "average error for " << area << "/" << dataset << ": " << arma::cdot(pi, error) << "\n";
//	cout << "error of 5050: " << error(5050) << "\n";
//	cout << "error of 5850: " << error(5850) << "\n";
//	cout << "error of 5800: " << error(5800) << "\n";
}


double x_privacy(const chan& C, const arma::ucolvec strategy, uint i) {
	return arma::cdot(C.row(i), arma::conv_to<arma::colvec>::from(strategy == i));
}


void output_x_privacy(const MinEntropy<double>& me, const mat& pois, const prob& pi) {

	auto str = me.strategy(pi % pois);

	for(uint i = 0; i < pi.n_elem; i++) {
		if(pi(i) > 0)
			cout << 1 - pois(i) * x_privacy(me.C, str, i) << "\n";
	}
}



void compute_elastic_privacy(string area, string dataset) {
	string areadataset = area + "-" + dataset;

	Mechanism<double> elastic;
	if(!elastic.C.load("temp/elastic-"+area+".bin")) {	

		mat dist;
		if(!dist.load("temp/metric-"+area+".bin")) {
			dist.load("marco/metric/"+area+".dat");
			dist.save("temp/metric-"+area+".bin");
		}
		Metric<double, uint> d = metrics::from_distance_matrix(dist);


		elastic = mechanisms::exponential(dist.n_rows, d);
		elastic.C.save("temp/elastic-"+area+".bin");

		dist.reset();
	}
	uint n = elastic.C.n_rows;


	mat priors;
	priors.load("marco/priors/"+areadataset+"-user.dat");

	prob prior_global;
	prior_global.load("marco/priors/"+areadataset+"-global.dat");


	mat pois = arma::ones(1, n);
	pois.load("marco/poi/"+areadataset+".dat");
	for(auto& e : pois)
		if(!equal(e, 0.0))
			e = 1/e;


	MinEntropy<double> me;
	share_memory(me.C, elastic.C);

	auto strategy = me.strategy(prior_global % pois);

	std::ofstream myfile;
	myfile.open("generated_data/elastic-" + areadataset);

	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);

		double privacy = 1.0;
		for(uint i = 0; i < n; i++)
			privacy -= pi(i) * pois(i) * x_privacy(me.C, strategy, i);
//		cout << arma::accu(pi > 0) << ": ";
		myfile << privacy << "\n";
	}
	myfile.close();
}


void compute_laplace_privacy(string area, string dataset) {
	string areadataset = area + "-" + dataset;

	uint grid_size = 100;
	double cell_width = 100;

	Mechanism<double> laplace;
	double eps = 2.0/1483;		// eps that gives ... avg error		
	if(!laplace.C.load("temp/laplace-"+std::to_string(eps)+".bin")) {

		laplace = mechanisms::planar_laplace_grid<double>(grid_size, grid_size, cell_width, eps);
		laplace.C.save("temp/laplace-"+std::to_string(eps)+".bin");
	}
	uint n = laplace.C.n_rows;


	mat priors;
	priors.load("marco/priors/"+areadataset+"-user.dat");

	prob prior_global;
	prior_global.load("marco/priors/"+areadataset+"-global.dat");


	mat pois = arma::ones(1, n);
	pois.load("marco/poi/"+areadataset+".dat");
	for(auto& e : pois)
		if(!equal(e, 0.0))
			e = 1/e;


	MinEntropy<double> me;
	share_memory(me.C, laplace.C);

	auto strategy = me.strategy(prior_global % pois);

	std::ofstream myfile;
	myfile.open("generated_data/laplace-" + areadataset);

	for(uint i = 0; i < priors.n_rows; i++) {
		prob pi = priors.row(i);

		double privacy = 1.0;
		for(uint i = 0; i < n; i++)
			privacy -= pi(i) * pois(i) * x_privacy(me.C, strategy, i);
//		cout << arma::accu(pi > 0) << ": ";
		myfile << privacy << "\n";
	}
	myfile.close();
}


int main() {
	arma::arma_rng::set_seed_random();

	std::set_terminate( my_terminate );

	
	compute_elastic_privacy("paris-nanterre", "gowalla");
	compute_elastic_privacy("paris-nanterre", "brightkite");
	compute_laplace_privacy("paris-nanterre", "gowalla");
	compute_laplace_privacy("paris-nanterre", "brightkite");

	compute_elastic_privacy("paris-center", "gowalla");
	compute_elastic_privacy("paris-center", "brightkite");
	compute_laplace_privacy("paris-center", "gowalla");
	compute_laplace_privacy("paris-center", "brightkite");

	return 0;


//	compute_avg_error(prior_global, area, dataset);
//	return 0;

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
//		cout << utility::expected_distance<double>(elastic.C, pi, Euclid) << "\n";
//	}
