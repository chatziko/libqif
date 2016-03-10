#include <iostream>
#include <string>

#include "qif"

using namespace qif;

void dp() {
	std::cout << "Using QIF Library Example" << std::endl;


	std::string channel_elements = "1 0 0; 0 1 0; 0 0 1";
	Mechanism<double> mechanism;
	mechanism.C = channel_elements;

	//this example asks with epsilon=0.05
	std::cout << "is_private: " << mechanism.is_private(0.05) << "\n";
}
