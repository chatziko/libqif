#include "tests_aux.h"

using namespace mechanism;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MechTest : public BaseTest<eT> {};

TYPED_TEST_CASE_P(MechTest);


TYPED_TEST_P(MechTest, Is_private) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Mech<eT> mech;
	mech.d = metric::euclidean<eT, uint>();

	mech.C = t.noint_10;
	EXPECT_TRUE(is_private(mech, eT(0)));

	mech.C = t.id_10;
	EXPECT_FALSE(is_private(mech, eT(10000)));

	mech.C = t.c1;
	EXPECT_TRUE (is_private(mech, std::log(eT(7.0)/2)));
	EXPECT_FALSE(is_private(mech, std::log(eT(6.9)/2)));
}

TYPED_TEST_P(MechTest, Smallest_epsilon) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Mech<eT> mech;
	mech.d = metric::euclidean<eT, uint>();

	mech.C = t.noint_10;
	EXPECT_PRED2(equal2<eT>, eT(0), smallest_epsilon(mech));

	mech.C = t.id_10;
	EXPECT_PRED2(equal2<eT>, infinity<eT>(), smallest_epsilon(mech));

	mech.C = t.c1;
	EXPECT_PRED2(equal2<eT>, std::log(7.0/2), smallest_epsilon(mech));
}

TYPED_TEST_P(MechTest, Reals) {
	typedef TypeParam eT;

	uint size = 10;
	eT step = 1.1;
	eT epsilon = 0.9;

	auto d = step * metric::euclidean<eT, uint>();

	Mech<eT> geom = mechanism::geometric<eT>(size, step, epsilon);
	Mech<eT> expon = mechanism::exponential<eT>(size, d, epsilon);
	Mech<eT> tc = mechanism::tight_constraints<eT>(size, d, epsilon);

	expect_channel(size, size, geom.C);
	EXPECT_TRUE(is_private(geom, epsilon));
	EXPECT_FALSE(is_private(geom, epsilon - eT(0.01)));
	EXPECT_PRED2(equal2<eT>, epsilon, smallest_epsilon(geom));

	expect_channel(size, size, expon.C);
	EXPECT_TRUE(is_private(expon, epsilon));
	EXPECT_PRED2(equal2<eT>, 0.64197307180467134, smallest_epsilon(expon));

	expect_channel(geom.C, tc.C);
}

TYPED_TEST_P(MechTest, Discrete) {
	typedef TypeParam eT;

	uint size = 5;
	eT epsilon = 0.9;

	auto d = metric::discrete<eT, uint>();

	Mech<eT> tc = mechanism::tight_constraints<eT>(size, d, epsilon);
	Mech<eT> expon = mechanism::exponential<eT>(size, d, 2*epsilon);

	expect_channel(size, size, tc.C);
	EXPECT_TRUE(is_private(tc, epsilon));
	EXPECT_FALSE(is_private(tc, epsilon - eT(0.01)));
	EXPECT_PRED2(equal2<eT>, epsilon, smallest_epsilon(tc));

	// the exponential mechanism with 2*epsilon should be the same as the tight constraints
	EXPECT_PRED2(chan_equal2<eT>, expon.C, tc.C);
}

TYPED_TEST_P(MechTest, Grid) {
	typedef TypeParam eT;

	uint width = 3,
		 height = 3,
		 size = width * height;
	eT step = 1.1,
	   epsilon = 0.9;

	auto d = step * metric::grid<eT, Point<eT>>(width);

	Mech<eT> laplace = mechanism::planar_laplace_grid<eT>(width, height, step, epsilon);
	Mech<eT> tc = mechanism::tight_constraints<eT>(size, d, epsilon);
	Mech<eT> expon = mechanism::exponential<eT>(size, d, epsilon);

	expect_channel(size, size, tc.C);
	EXPECT_TRUE(is_private(tc, epsilon));
	EXPECT_FALSE(is_private(tc, epsilon - eT(0.01)));
	EXPECT_PRED2(equal2<eT>, epsilon, smallest_epsilon(tc));

	// for the planar laplace, the accuracy that we get through numeric integration is not that great
	EXPECT_PRED2(channel::is_proper<eT>, laplace.C, eT(1e-3));
	EXPECT_TRUE(is_private(laplace, epsilon));
	EXPECT_PRED4(equal<eT>, 0.8844, smallest_epsilon(laplace), 0, eT(1e-3));

	expect_channel(size, size, expon.C);
	EXPECT_TRUE(is_private(expon, epsilon));
	EXPECT_PRED2(equal2<eT>, 0.58945591528726249, smallest_epsilon(expon));
}



REGISTER_TYPED_TEST_CASE_P(MechTest, Is_private, Smallest_epsilon, Reals, Discrete, Grid);

INSTANTIATE_TYPED_TEST_CASE_P(Mech, MechTest, NativeTypes);

