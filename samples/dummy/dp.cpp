#include <iostream>
#include <string>

#include "qif"

using namespace qif;

void dp() {
	std::cout << "Using QIF Library Example" << std::endl;


	std::string channel_elements = "1 0 0; 0 1 0; 0 0 1";
	chan C = channel_elements;
	auto d = metric::euclidean<double, uint>();

	//this example asks with epsilon=0.05
	std::cout << "is_private: " << measure::d_privacy::is_private(C, 0.05 * d) << "\n";
}
