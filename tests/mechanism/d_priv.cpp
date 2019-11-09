#include "tests_aux.h"

using namespace mechanism::d_priv;
using namespace measure::d_priv;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MechDPrivTest : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(MechDPrivTest);


TYPED_TEST_P(MechDPrivTest, Reals) {
	typedef TypeParam eT;

	uint size = 10;
	eT step = 1.1;
	eT epsilon = 0.9;

	auto d = step * metric::euclidean<eT, uint>();

	Chan<eT> geom = geometric<eT>(size, epsilon * step);
	Chan<eT> expon = exponential<eT>(size, epsilon * d);
	Chan<eT> tc = tight_constraints<eT>(size, epsilon * d);

	// tight constraints with first/last and middle coeffs forced to be equal
	arma::uvec cols(size, arma::fill::ones);
	cols(0) = cols(size-1) = 0;
	Chan<eT> tc2 = tight_constraints<eT>(cols, epsilon * d);

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, geom, size, size);
	EXPECT_TRUE(is_private(geom, epsilon * d));
	EXPECT_FALSE(is_private(geom, (epsilon - eT(0.01)) * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(geom, d));

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, expon, size, size);
	EXPECT_TRUE(is_private(expon, epsilon * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0.64197307180467134, smallest_epsilon(expon, d));

	EXPECT_PRED_FORMAT2(chan_equal2<eT>, geom, tc);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, geom, tc2);

	// fat
	geom = geometric<eT>(size, epsilon * step, 2*size);
	expon = exponential<eT>(size, epsilon * d, 2*size);

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, geom, size, 2*size);
	if(!std::is_same<eT, float>::value) { // not-enough precision
		EXPECT_TRUE(is_private(geom, epsilon * d));
		EXPECT_FALSE(is_private(geom, (epsilon - eT(0.01)) * d));
		EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(geom, d));
	}

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, expon, size, 2*size);
	EXPECT_TRUE(is_private(expon, epsilon * d));

	// skinny
	geom = geometric<eT>(2*size, epsilon * step, size);
	expon = exponential<eT>(2*size, epsilon * d, size);

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, geom, 2*size, size);
	if(!std::is_same<eT, float>::value) { // not-enough precision
		EXPECT_TRUE(is_private(geom, epsilon * d));
		EXPECT_FALSE(is_private(geom, (epsilon - eT(0.01)) * d));
		EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(geom, d));
	}

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, expon, 2*size, size);
	EXPECT_TRUE(is_private(expon, epsilon * d));

	// geometric, different truncations
	// We create 4 channels with secrets 5..14 and obs 0..9 / 10..19 / 1..18 / 7..12
	// We first create by remapping a 0..19 x 0..19 channel, then create directly and compare
	//
	Chan<eT> G0 = geometric(20, epsilon);		// the 0..19 x 0..19 to remap

	Chan<eT> T1 = G0 * channel::deterministic<eT>([=](uint x) -> uint { return std::min(         x      ,  9u); }, 20, 20);	// 0..9
	Chan<eT> T2 = G0 * channel::deterministic<eT>([=](uint x) -> uint { return std::min(std::max(x, 10u), 19u); }, 20, 20);	// 10..19
	Chan<eT> T3 = G0 * channel::deterministic<eT>([=](uint x) -> uint { return std::min(std::max(x,  1u), 18u); }, 20, 20);	// 1..18
	Chan<eT> T4 = G0 * channel::deterministic<eT>([=](uint x) -> uint { return std::min(std::max(x,  7u), 12u); }, 20, 20);	// 7..12

	Chan<eT> G1 = T1.submat(5, 0,  14,  9);		// crop rows to
	Chan<eT> G2 = T2.submat(5, 10, 14, 19);		// 5..14
	Chan<eT> G3 = T3.submat(5, 1,  14, 18);		// and the columns to
	Chan<eT> G4 = T4.submat(5, 7,  14, 12);		// 0..9 / 10..19 / 1..18 / 7..12

	Chan<eT> D1 = geometric(10, epsilon, 10, 5, 0 );	// same channels
	Chan<eT> D2 = geometric(10, epsilon, 10, 5, 10);	// created
	Chan<eT> D3 = geometric(10, epsilon, 18, 5, 1 );	// directly by passing
	Chan<eT> D4 = geometric(10, epsilon, 6,  5, 7 );	// first_y / first_x

	EXPECT_PRED_FORMAT2(chan_equal2<eT>, G1, D1);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, G2, D2);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, G3, D3);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, G4, D4);
}

TYPED_TEST_P(MechDPrivTest, Discrete) {
	typedef TypeParam eT;

	uint size = 5;
	eT epsilon = 0.9;

	auto d = metric::discrete<eT, uint>();

	Chan<eT> tc = tight_constraints<eT>(size, epsilon * d);
	Chan<eT> expon = exponential<eT>(size, 2 * epsilon * d);
	Chan<eT> rr = randomized_response<eT>(size, epsilon);

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, tc, size, size);
	EXPECT_TRUE(is_private(tc, epsilon * d));
	EXPECT_FALSE(is_private(tc, (epsilon - eT(0.01)) * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(tc, d));

	// the exponential mechanism with 2*epsilon should be the same as the tight constraints and randomized response
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, expon, tc);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, expon, rr);
}

TYPED_TEST_P(MechDPrivTest, Grid) {
	typedef TypeParam eT;

	uint width = 3,
		 height = 3,
		 size = width * height;
	eT step = 1.1,
	   epsilon = 0.9;

	auto d = step * metric::grid<eT>(width);

	Chan<eT> tc = tight_constraints<eT>(size, epsilon * d);
	Chan<eT> expon = exponential<eT>(size, epsilon * d);

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, tc, size, size);
	EXPECT_TRUE(is_private(tc, epsilon * d));
	EXPECT_FALSE(is_private(tc, (epsilon - eT(0.01)) * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(tc, d));

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, expon, size, size);
	EXPECT_TRUE(is_private(expon, epsilon * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0.58945591528726249, smallest_epsilon(expon, d));
}



REGISTER_TYPED_TEST_SUITE_P(MechDPrivTest, Reals, Discrete, Grid);

INSTANTIATE_TYPED_TEST_SUITE_P(Mech, MechDPrivTest, NativeTypes);

