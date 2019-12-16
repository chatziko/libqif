#include "tests_aux.h"

using namespace mechanism::d_privacy;
using namespace measure::d_privacy;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MechGeoTest : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(MechGeoTest);


TYPED_TEST_P(MechGeoTest, Grid) {
	typedef TypeParam eT;

	uint width = 3,
		 height = 3,
		 size = width * height;
	eT step = 1.1,
	   epsilon = 0.9;

	auto d = step * metric::grid<eT>(width);

	Chan<eT> laplace = mechanism::geo::planar_laplace_grid<eT>(width, height, step, epsilon);
	Chan<eT> geom = mechanism::geo::planar_geometric_grid<eT>(width, height, step, epsilon);

	// for the planar laplace, the accuracy that we get through numeric integration is not that great
	EXPECT_PRED_FORMAT2(chan_is_proper2<eT>, laplace, eT(1e-3));
	EXPECT_TRUE(is_private(laplace, epsilon * d));
	EXPECT_PRED_FORMAT4(equal4<eT>, 0.8844, smallest_epsilon(laplace, d), 0, eT(1e-3));

	// same (just a bit better) for planar_geometric
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, geom, size, size);
	EXPECT_TRUE(is_private(geom, (epsilon+eT(1e-5)) * d));
	EXPECT_PRED_FORMAT4(equal4<eT>, epsilon, smallest_epsilon(geom, d), 0, eT(1e-6));
}



REGISTER_TYPED_TEST_SUITE_P(MechGeoTest, Grid);

INSTANTIATE_TYPED_TEST_SUITE_P(Mech, MechGeoTest, NativeTypes);

