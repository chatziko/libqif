#include "tests_aux.h"

#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision


int main(int argc, char **argv) {
	qif::rng::set_seed_random();

	std::cout << std::setprecision(100);

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
