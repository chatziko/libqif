#include "tests_aux.h"

using namespace measure::d_privacy;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MeasureDPrivTest : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(MeasureDPrivTest);


TYPED_TEST_P(MeasureDPrivTest, Is_private) {
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

TYPED_TEST_P(MeasureDPrivTest, Smallest_epsilon) {
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


REGISTER_TYPED_TEST_SUITE_P(MeasureDPrivTest, Is_private, Smallest_epsilon);

INSTANTIATE_TYPED_TEST_SUITE_P(Mech, MeasureDPrivTest, NativeTypes);

