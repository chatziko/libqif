#include <armadillo>
#include "gtest/gtest.h"


int main(int argc, char **argv) {
	arma::arma_rng::set_seed_random();

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
