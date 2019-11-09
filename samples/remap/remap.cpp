
#include <string>
#include <iostream>     // std::cout, std::fixed

using std::cout;
using std::cerr;

#include "qif"
using namespace qif;
using namespace qif::measure;


int main(int argc, char** argv) {

	// parse and check input
	//
	if(argc != 4) {
		cerr << "Usage: remap <C> <pi> <L>\n";
		return 1;
	}

	chan C;
	prob pi;
	mat L;

	if(!C.load(argv[1])) {
		cerr << "Cannot load C from " << argv[1] << "\n";
		return 1;
	}
	if(!pi.load(argv[2])) {
		cerr << "Cannot load pi from " << argv[2] << "\n";
		return 1;
	}
	if(!L.load(argv[3])) {
		cerr << "Cannot load L from " << argv[3] << "\n";
		return 1;
	}

	if(!channel::is_proper(C))
		cerr << "Warning: C might not be a proper channel\n";
	if(!probab::is_proper(pi))
		cerr << "Warning: pi might not be a proper distribution\n";

	// compute remap and print
	//
	chan R = channel::deterministic<double>(l_uncert::strategy(L, pi, C), L.n_rows);
	chan CR = C * R;

	cout.setf(std::ios::fixed);
	cout.precision(6);
	CR.raw_print(cout);
}
