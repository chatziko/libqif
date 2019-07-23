#include "tests_aux.h"

using namespace mechanism;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MechOptTest : public BaseTest<eT> {};

TYPED_TEST_CASE_P(MechOptTest);



TYPED_TEST_P(MechOptTest, OptimalUtility) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	uint size = 10;
	eT step = 1.1;
	eT epsilon = 0.9;

	// with euclidean metric we get the geometric mechanism
	auto d = step * metric::euclidean<eT, uint>();
	auto bin_loss = metric::discrete<eT, uint>();

	Chan<eT> geom = mechanism::geometric<eT>(size, epsilon * step);
	Chan<eT> opt = mechanism::optimal_utility<eT>(t.unif_10, size, epsilon * d, bin_loss);
	Chan<eT> dist_opt = mechanism::dist_optimal_utility<eT>(t.unif_10, size, epsilon * d, bin_loss);
	Chan<eT> dist_opt2 = mechanism::dist_optimal_utility_strict<eT>(t.unif_10, size, epsilon * d, bin_loss);

	EXPECT_PRED_FORMAT2(chan_equal2<eT>, geom, opt);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, geom, dist_opt);

	EXPECT_PRED_FORMAT2(chan_is_proper2<eT>, dist_opt2, eT(1e-3));		// not sure why we don't get good precision

	// with inf constraints we get the identity channel
	eT inf = infinity<eT>();
	d = inf * metric::discrete<eT, uint>();

	opt = mechanism::optimal_utility<eT>(t.unif_10, size, epsilon * d, bin_loss);

	EXPECT_PRED_FORMAT2(chan_equal2<eT>, t.id_10, opt);
}



REGISTER_TYPED_TEST_CASE_P(MechOptTest, OptimalUtility);

INSTANTIATE_TYPED_TEST_CASE_P(MechOpt, MechOptTest, NativeTypes);

