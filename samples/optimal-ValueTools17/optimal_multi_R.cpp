#include "qif"
#include <iostream>
#include <fstream>
#include <time.h>
#include <python2.7/Python.h>

using std::cout;
using std::cin;
using std::vector;
using std::string;
using std::sqrt;
using namespace qif;


void compute_optimal(string method, uint W, double radius, double delta) {
	if(method != "direct" && method != "kostas" && method != "geometric")
		throw std::runtime_error("invalid method: " + method);

	lp::Defaults::method = lp::method_t::simplex_dual;
	lp::Defaults::glp_msg_level = lp::msg_level_t::all;
	lp::Defaults::glp_presolve = true;

	double unit = 1;
	double eps = std::log(2) / (2*unit);
	uint width = W,
		 n_inputs = width * width,
		 n_outputs = n_inputs;
	double R = unit * radius;
	double rho = unit * sqrt(2)/2;	

	auto d_grid = unit * metric::grid<double>(width);	// d_grid(i, j) = euclidean distance between cell indices i,j
	auto d_loss = d_grid;                              // loss metric = euclidean
	auto dx = d_grid;                                 // privacy metric


	if(method == "kostas") {
		// use more relexed dx and update eps
		dx = metric::threshold_inf(d_grid, R);    // inf if above threshold R

		cout << "R: " << R << "\n";
		cout << "delta: " << delta << "\n";
		eps /= delta;     // update eps to compensate for the use of a bigger dx
		cout << "eps/delta: " << eps << "\n";
	}

	clock_t t_start = clock();
	dx.chainable = [](const uint&, const uint&) -> bool { return false; };
	prob pi = probab::uniform<double>(n_inputs);		// uniform prior
	chan C = method == "geometric"
		? mechanism::planar_geometric_grid<double>(width, width, unit, eps)
		: mechanism::optimal_utility(pi, n_outputs, eps * dx, d_loss);
        clock_t comp_time = (clock() - t_start)/CLOCKS_PER_SEC;
	double time_min = (double) comp_time/60.0;

//	cout << "Number of Locations: " << C.n_rows << "\n";
//	cout << "proper: " << channel::is_proper(C) << "\n";
//	cout << "Utility Loss : " << l::post_entropy(d_loss, pi, C) << "\n";	// expected d_loss between x and y

	// apply optimal remap (probably does nothing for uniform priors)
	mat Loss = l::metric_to_mat(d_loss, n_inputs);
	chan Remap = channel::deterministic<double>(l::strategy(Loss, pi, C), Loss.n_rows);	// compute remap
	cout << "Utility Loss after remap : " << l::post_entropy(d_loss, pi, (chan)(C*Remap)) << "\n";
       
	// append result to file
	std::ofstream myfile;
        myfile.open ("results.txt", std::ios::app);
        myfile << C.n_rows << "\t" << channel::is_proper(C) << "\t"  << l::post_entropy(d_loss, pi, C) << "\t\t" << time_min << "\n";
        myfile.close();
}

int main(int argc, char *argv[]) {
	vector<string> args(argv+1, argv + argc);
	
    string Labels[4] = {"c=2.8", "c=3", "c=3.4", "c=4.2"};
    double Rad[4] = {sqrt(98)/5, sqrt(18)/2, 1.7*sqrt(2), 2.1*sqrt(2)};
    double Del[4] = {1.243833713, 1.199065696, 1.149582276, 1.109489965};


    string m_name;
    cout << "Please, choose a method among: 'kostas', 'geometric', or 'direct'.\n";
    cin >> m_name;
  
    uint w_min, w_max;
    cout << "Enter the width of the smallest grid to consider.\n";
    cin >> w_min;
    cout << "Enter the width of the largest grid to consider.\n";
    cin >> w_max;

    // write attribute names into first line of "results.txt"
    std::ofstream myfile;
    myfile.open ("results.txt");
    myfile << "N\tProper\tUL\t\tT(m)\n";
    myfile.close();


	for (int i = 0; i < 4; i++)
    {
        // append label to "results.txt" to separate results
        std::ofstream myfile;
        myfile.open ("results.txt", std::ios::app);
        myfile << "flag\t"+Labels[i]+"\n";
        myfile.close();
		
	    // for the currecnt c value = $Labels[i], iterate over grid's width from $w_min to $w_max
	    for(uint n = w_min; n <= w_max; n++)
	    {
        	compute_optimal(m_name, n, Rad[i], Del[i]);
	    }
    }
	
	Py_SetProgramName(argv[0]);
	Py_Initialize();
	PyRun_SimpleString("exec(open('plot.py').read())");
	Py_Finalize();
}


