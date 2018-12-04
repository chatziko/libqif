#include "tests_aux.h"

using namespace mechanism;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MechTest : public BaseTest<eT> {};

TYPED_TEST_CASE_P(MechTest);


TYPED_TEST_P(MechTest, Is_private) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Chan<eT> C;
	auto d = metric::euclidean<eT, uint>();

	C = t.noint_10;
	EXPECT_TRUE(is_private(C, eT(0)*d));

	C = t.id_10;
	EXPECT_FALSE(is_private(C, eT(10000)*d));

	C = t.c1;
	EXPECT_TRUE (is_private(C, std::log(eT(7.0)/2)*d));
	EXPECT_FALSE(is_private(C, std::log(eT(6.9)/2)*d));
}

TYPED_TEST_P(MechTest, Smallest_epsilon) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Chan<eT> C;
	auto d = metric::euclidean<eT, uint>();

	C = t.noint_10;
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), smallest_epsilon(C, d));

	C = t.id_10;
	EXPECT_PRED_FORMAT2(equal2<eT>, infinity<eT>(), smallest_epsilon(C, d));

	C = t.c1;
	EXPECT_PRED_FORMAT2(equal2<eT>, std::log(7.0/2), smallest_epsilon(C, d));
}

TYPED_TEST_P(MechTest, Reals) {
	typedef TypeParam eT;

	uint size = 10;
	eT step = 1.1;
	eT epsilon = 0.9;

	auto d = step * metric::euclidean<eT, uint>();

	Chan<eT> geom = mechanism::geometric<eT>(size, epsilon * d);
	Chan<eT> expon = mechanism::exponential<eT>(size, epsilon * d);
	Chan<eT> tc = mechanism::tight_constraints<eT>(size, epsilon * d);

	// tight constraints with first/last and middle coeffs forced to be equal
	arma::uvec cols(size, arma::fill::ones);
	cols(0) = cols(size-1) = 0;
	Chan<eT> tc2 = mechanism::tight_constraints<eT>(cols, epsilon * d);

	expect_channel(size, size, geom);
	EXPECT_TRUE(is_private(geom, epsilon * d));
	EXPECT_FALSE(is_private(geom, (epsilon - eT(0.01)) * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(geom, d));

	expect_channel(size, size, expon);
	EXPECT_TRUE(is_private(expon, epsilon * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0.64197307180467134, smallest_epsilon(expon, d));

	expect_channel(geom, tc );
	expect_channel(geom, tc2);

	// fat
	geom = mechanism::geometric<eT>(size, epsilon * d, 2*size);
	expon = mechanism::exponential<eT>(size, epsilon * d, 2*size);

	expect_channel(size, 2*size, geom);
	if(!std::is_same<eT, float>::value) { // not-enough precision
		EXPECT_TRUE(is_private(geom, epsilon * d));
		EXPECT_FALSE(is_private(geom, (epsilon - eT(0.01)) * d));
		EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(geom, d));
	}

	expect_channel(size, 2*size, expon);
	EXPECT_TRUE(is_private(expon, epsilon * d));

	// skinny
	geom = mechanism::geometric<eT>(2*size, epsilon * d, size);
	expon = mechanism::exponential<eT>(2*size, epsilon * d, size);

	expect_channel(2*size, size, geom);
	if(!std::is_same<eT, float>::value) { // not-enough precision
		EXPECT_TRUE(is_private(geom, epsilon * d));
		EXPECT_FALSE(is_private(geom, (epsilon - eT(0.01)) * d));
		EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(geom, d));
	}

	expect_channel(2*size, size, expon);
	EXPECT_TRUE(is_private(expon, epsilon * d));
}

TYPED_TEST_P(MechTest, Discrete) {
	typedef TypeParam eT;

	uint size = 5;
	eT epsilon = 0.9;

	auto d = metric::discrete<eT, uint>();

	Chan<eT> tc = mechanism::tight_constraints<eT>(size, epsilon * d);
	Chan<eT> expon = mechanism::exponential<eT>(size, 2 * epsilon * d);

	expect_channel(size, size, tc);
	EXPECT_TRUE(is_private(tc, epsilon * d));
	EXPECT_FALSE(is_private(tc, (epsilon - eT(0.01)) * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(tc, d));

	// the exponential mechanism with 2*epsilon should be the same as the tight constraints
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, expon, tc);
}

TYPED_TEST_P(MechTest, Grid) {
	typedef TypeParam eT;

	uint width = 3,
		 height = 3,
		 size = width * height;
	eT step = 1.1,
	   epsilon = 0.9;

	auto d = step * metric::grid<eT>(width);

	Chan<eT> laplace = mechanism::planar_laplace_grid<eT>(width, height, step, epsilon);
	Chan<eT> geom = mechanism::planar_geometric_grid<eT>(width, height, step, epsilon);
	Chan<eT> tc = mechanism::tight_constraints<eT>(size, epsilon * d);
	Chan<eT> expon = mechanism::exponential<eT>(size, epsilon * d);

	expect_channel(size, size, tc);
	EXPECT_TRUE(is_private(tc, epsilon * d));
	EXPECT_FALSE(is_private(tc, (epsilon - eT(0.01)) * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, epsilon, smallest_epsilon(tc, d));

	// for the planar laplace, the accuracy that we get through numeric integration is not that great
	EXPECT_PRED_FORMAT2(chan_is_proper2<eT>, laplace, eT(1e-3));
	EXPECT_TRUE(is_private(laplace, epsilon * d));
	EXPECT_PRED_FORMAT4(equal4<eT>, 0.8844, smallest_epsilon(laplace, d), 0, eT(1e-3));

	// same (just a bit better) for planar_geometric
	expect_channel(size, size, geom);
	EXPECT_TRUE(is_private(geom, (epsilon+eT(1e-5)) * d));
	EXPECT_PRED_FORMAT4(equal4<eT>, epsilon, smallest_epsilon(geom, d), 0, eT(1e-6));

	expect_channel(size, size, expon);
	EXPECT_TRUE(is_private(expon, epsilon * d));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0.58945591528726249, smallest_epsilon(expon, d));
}



REGISTER_TYPED_TEST_CASE_P(MechTest, Is_private, Smallest_epsilon, Reals, Discrete, Grid);

INSTANTIATE_TYPED_TEST_CASE_P(Mech, MechTest, NativeTypes);

