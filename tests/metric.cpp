#include "tests_aux.h"


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class MetricTest : public BaseTest<eT> {};
template <typename eT>
class MetricTestReals : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(MetricTest);
TYPED_TEST_SUITE_P(MetricTestReals);		// tests that run only on double/float


TYPED_TEST_P(MetricTest, Euclidean_uint) {
	typedef TypeParam eT;

	auto euclid = metric::euclidean<eT, uint>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), euclid(3, 3));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(5), euclid(0, 5));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(5), euclid(5, 0));

	EXPECT_FALSE(euclid.chainable(0, 1));
	EXPECT_FALSE(euclid.chainable(5, 6));
	EXPECT_TRUE(euclid.chainable(0, 5));
}

TYPED_TEST_P(MetricTest, Discrete) {
	typedef TypeParam eT;

	auto disc = metric::discrete<eT, uint>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), disc(3, 3));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1), disc(0, 2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1), disc(5, 0));
	EXPECT_FALSE(disc.chainable(0, 1));
	EXPECT_FALSE(disc.chainable(0, 2));
}

TYPED_TEST_P(MetricTest, Scale) {
	typedef TypeParam eT;

	auto scaled_euclid = eT(10) * metric::euclidean<eT, uint>();
	auto scaled_disc   = eT(10) * metric::discrete <eT, uint>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), scaled_euclid(3, 3));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(50), scaled_euclid(0, 5));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(50), scaled_euclid(5, 0));
	EXPECT_FALSE(scaled_euclid.chainable(0, 1));
	EXPECT_FALSE(scaled_euclid.chainable(5, 6));
	EXPECT_TRUE(scaled_euclid.chainable(0, 5));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), scaled_disc(3, 3));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(10), scaled_disc(0, 3));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(10), scaled_disc(5, 0));
	EXPECT_FALSE(scaled_disc.chainable(0, 1));
	EXPECT_FALSE(scaled_disc.chainable(0, 2));
}

TYPED_TEST_P(MetricTest, Threshold) {
	typedef TypeParam eT;

	auto thres_euclid = metric::threshold_bin(metric::euclidean<eT, uint>(), eT(10)) ;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), thres_euclid(3, 3));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), thres_euclid(0, 5));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1), thres_euclid(0, 10));
	EXPECT_FALSE(thres_euclid.chainable(0, 1));
	EXPECT_TRUE(thres_euclid.chainable(0, 2));
	EXPECT_FALSE(thres_euclid.chainable(0, 10));
}

TYPED_TEST_P(MetricTestReals, Euclidean_point) {
	typedef TypeParam eT;
	typedef Point<eT> P;

	auto euclid = metric::euclidean<eT, P>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),            euclid(P(1, 1), P(1, 1)));
	EXPECT_PRED_FORMAT2(equal2<eT>, std::sqrt(eT(2)), euclid(P(0, 0), P(1, 1)));
	EXPECT_PRED_FORMAT2(equal2<eT>, std::sqrt(eT(5)), euclid(P(2, 3), P(1, 1)));

	EXPECT_FALSE(euclid.chainable(P(0, 0), P(1, 1)));
	EXPECT_FALSE(euclid.chainable(P(2, 3), P(1, 1)));
}

TYPED_TEST_P(MetricTest, Manhattan_point) {
	typedef TypeParam eT;
	typedef Point<eT> P;

	auto manh = metric::manhattan<eT, P>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), manh(P(1, 1), P(1, 1)));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(2), manh(P(0, 0), P(1, 1)));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(3), manh(P(2, 3), P(1, 1)));

	EXPECT_FALSE(manh.chainable(P(0, 0), P(1, 1)));
	EXPECT_FALSE(manh.chainable(P(2, 3), P(1, 1)));
}

TYPED_TEST_P(MetricTestReals, Grid_point) {
	typedef TypeParam eT;

	auto grid_euclid = metric::grid<eT>(4);
	auto grid_manh   = metric::grid<eT>(4, metric::manhattan<eT, Point<uint>>());

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),            grid_euclid(5, 5));
	EXPECT_PRED_FORMAT2(equal2<eT>, std::sqrt(eT(2)), grid_euclid(0, 5));		// cell  5 is (1,1)
	EXPECT_PRED_FORMAT2(equal2<eT>, std::sqrt(eT(5)), grid_euclid(11, 5));		// cell 11 is (2,3)

	EXPECT_FALSE (grid_euclid.chainable(0, 5));								// cell  5 is (1,1)
	EXPECT_FALSE (grid_euclid.chainable(11, 5));								// cell 11 is (2,3)
	EXPECT_TRUE(grid_euclid.chainable(0, 2));								// cell  2 is (0,2)
	EXPECT_TRUE(grid_euclid.chainable(0, 10));								// cell 10 is (2,2)
	EXPECT_TRUE(grid_euclid.chainable(0, 8));								// cell  8 is (2,0)

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), grid_manh(5, 5));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(2), grid_manh(0, 5));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(3), grid_manh(11, 5));

	EXPECT_FALSE (grid_manh.chainable(0, 5));									// cell  5 is (1,1)
	EXPECT_TRUE(grid_manh.chainable(11, 5));									// cell 11 is (2,3)
	EXPECT_TRUE(grid_manh.chainable(0, 2));									// cell  2 is (0,2)
	EXPECT_TRUE(grid_manh.chainable(0, 10));									// cell 10 is (2,2)
	EXPECT_TRUE(grid_manh.chainable(0, 8));									// cell  8 is (2,0)
}

TYPED_TEST_P(MetricTest, Total_variation) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	auto tv = metric::total_variation<eT, Prob<eT>>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), tv(t.unif_4, t.unif_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), tv(t.dirac_4, t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(3)/4, tv(t.unif_4, t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(3)/4, tv(t.dirac_4, t.unif_4));
}

TYPED_TEST_P(MetricTest, Convex_separation) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	auto bed = metric::convex_separation<eT, Prob<eT>>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), bed(t.unif_4, t.unif_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), bed(t.dirac_4, t.dirac_4));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1), bed(t.unif_4, t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1), bed(t.dirac_4, t.unif_4));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(9)/eT(14), bed(t.unif_4, t.pi5));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(9)/eT(14), bed(t.pi5, t.unif_4));
}

TYPED_TEST_P(MetricTestReals, Multiplicative_distance) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	auto mtv = metric::mult_total_variation<eT, Prob<eT>>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), mtv(t.unif_4, t.unif_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0), mtv(t.dirac_4, t.dirac_4));

	EXPECT_PRED_FORMAT2(equal2<eT>, infinity<eT>(), mtv(t.unif_4, t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, infinity<eT>(), mtv(t.dirac_4, t.unif_4));

	EXPECT_PRED_FORMAT2(equal2<eT>, std::log(0.7/0.25), mtv(t.unif_4, t.pi5));
	EXPECT_PRED_FORMAT2(equal2<eT>, std::log(0.7/0.25), mtv(t.pi5, t.unif_4));
}

TYPED_TEST_P(MetricTest, Kantorovich) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	bool is_double		= std::is_same<eT, double>::value;
	auto disc			= metric::discrete<eT, uint>(),
		 euclid			= metric::euclidean<eT, uint>();
	auto kant_disc		= metric::kantorovich<eT, Prob<eT>>(disc),
		 kant_euclid	= metric::kantorovich<eT, Prob<eT>>(euclid),
		 kant_lp_disc	= metric::kantorovich_lp<eT, Prob<eT>>(disc),
		 kant_lp_euclid	= metric::kantorovich_lp<eT, Prob<eT>>(euclid),
		 tv				= metric::total_variation<eT, Prob<eT>>();

	for (bool use_lp : { false, true }) {

		auto kant_cur_disc   = use_lp ? kant_lp_disc   : kant_disc;
		auto kant_cur_euclid = use_lp ? kant_lp_euclid : kant_euclid;

		// metric::kantorovich runs kantorovich_fastemd for double and kantorovich_lp for the rest
		// fastemd has worse precision, so we need to set a larger mrd
		eT mrd = is_double && !use_lp ? 1e-5 : def_mrd<eT>;

		EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),   kant_cur_disc(t.unif_4, t.unif_4));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),   kant_cur_disc(t.dirac_4, t.dirac_4));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(3)/4, kant_cur_disc(t.unif_4, t.dirac_4));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(3)/4, kant_cur_disc(t.dirac_4, t.unif_4));

		EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),   kant_cur_euclid(t.unif_4, t.unif_4));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),   kant_cur_euclid(t.dirac_4, t.dirac_4));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(3)/2, kant_cur_euclid(t.unif_4, t.dirac_4));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(3)/2, kant_cur_euclid(t.dirac_4, t.unif_4));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(9)/2, kant_cur_euclid(t.unif_10, t.dirac_10));
		EXPECT_PRED_FORMAT2(equal2<eT>, eT(9)/2, kant_cur_euclid(t.dirac_10, t.unif_10));

		for(uint i = 0; i < 3; i++) {
			// distance of dirac dists is the same as the distance between the corresponding elements
			EXPECT_PRED_FORMAT2(equal2<eT>, disc  (0, i), kant_cur_disc  (t.dirac_4, probab::dirac<eT>(4, i)));
			EXPECT_PRED_FORMAT4(equal4<eT>, euclid(0, i), kant_cur_euclid(t.dirac_4, probab::dirac<eT>(4, i)), 0, mrd);
		}

		// kantorovich over the discrete metric = total variation
		//
		// Note: When eT=float, this fails under GLOP which is now the default solver so we disable.
		//       It also fails when randu uses the "naif normalize" algorithm,
		//       probably the sum of each dist is slighly diffent than 1.0, so the transportation problem is infeasible.
		//       When randu uses the "differences of sorted list" algorithm the test always passes (whith CLP).
		//
		if(is_double) {
			for(uint i = 0; i < 10; i++) {
				auto p1 = probab::randu<eT>(10),
					p2 = probab::randu<eT>(10);
				EXPECT_PRED_FORMAT4(equal4<eT>, tv(p1, p2), kant_cur_disc(p1, p2), 0, mrd);
			}
		}
	}

	// kantorovich and kantorovich_lp should produce the same result
	if(is_double) {
		for(uint i = 0; i < 10; i++) {
			auto p1 = probab::randu<eT>(10),
				 p2 = probab::randu<eT>(10);
			EXPECT_PRED_FORMAT4(equal4<eT>, kant_disc  (p1, p2), kant_lp_disc  (p1, p2), 0, 1e-5);
			EXPECT_PRED_FORMAT4(equal4<eT>, kant_euclid(p1, p2), kant_lp_euclid(p1, p2), 0, 1e-4);
		}
	}
}

TYPED_TEST_P(MetricTestReals, Mult_kantorovich) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	auto disc			= infinity<eT>() * metric::discrete<eT, uint>(),
		 euclid			= metric::euclidean<eT, uint>();
	auto mkant_disc		= metric::mult_kantorovich<eT, Prob<eT>>(disc),
		 mkant_euclid	= metric::mult_kantorovich<eT, Prob<eT>>(euclid),
		 mtv			= metric::mult_total_variation<eT, Prob<eT>>();

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),              mkant_disc(t.unif_4,  t.unif_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),              mkant_disc(t.dirac_4, t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, infinity<eT>(),     mkant_disc(t.unif_4,  t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, infinity<eT>(),     mkant_disc(t.dirac_4, t.unif_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, std::log(0.7/0.25), mkant_disc(t.unif_4,  t.pi5));
	EXPECT_PRED_FORMAT2(equal2<eT>, std::log(0.7/0.25), mkant_disc(t.pi5,     t.unif_4));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),                  mkant_euclid(t.unif_4,   t.unif_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),                  mkant_euclid(t.dirac_4,  t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(2.0538953374413045), mkant_euclid(t.unif_4,   t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(2.0538953374413045), mkant_euclid(t.dirac_4,  t.unif_4));
	EXPECT_PRED_FORMAT4(equal4<eT>, eT(7.156044651432), mkant_euclid(t.unif_10,  t.dirac_10), 1e-12, 0);
	EXPECT_PRED_FORMAT4(equal4<eT>, eT(7.156044651432), mkant_euclid(t.dirac_10, t.unif_10), 1e-12, 0);

	for(uint i = 0; i < 3; i++) {
		// distance of dirac dists is the same as the distance between the corresponding elements
		EXPECT_PRED_FORMAT2(equal2<eT>, disc  (0, i), mkant_disc  (t.dirac_4, probab::dirac<eT>(4, i)));
		EXPECT_PRED_FORMAT2(equal2<eT>, euclid(0, i), mkant_euclid(t.dirac_4, probab::dirac<eT>(4, i)));
	}

	// mult kantorovich over the discrete metric (with inf value) = mult total variation
	auto p1 = probab::randu<eT>(10),
		 p2 = probab::randu<eT>(10);
	EXPECT_PRED_FORMAT2(equal2<eT>, mtv(p1, p2), mkant_disc(p1, p2));
}


// run the MetricTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_SUITE_P(MetricTest, Euclidean_uint, Scale, Threshold, Discrete, Manhattan_point, Total_variation, Convex_separation, Kantorovich);
REGISTER_TYPED_TEST_SUITE_P(MetricTestReals, Euclidean_point, Grid_point, Multiplicative_distance, Mult_kantorovich);

INSTANTIATE_TYPED_TEST_SUITE_P(Metric, MetricTest, AllTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Metric, MetricTestReals, NativeTypes);

