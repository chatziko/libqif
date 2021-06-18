#include "tests_aux.h"

using namespace probab;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ProbTest : public BaseTest<T> {};

TYPED_TEST_SUITE_P(ProbTest);


TYPED_TEST_P(ProbTest, Construct) {
	typedef TypeParam eT;

	std::string s = format_num<eT>("0.5 0.25 0.25");
	Prob<eT> pi = { eT(1)/2, eT(1)/4, eT(1)/4 };

	EXPECT_PRED_FORMAT2(prob_equal2<eT>, Prob<eT>(s.c_str()),  pi); // char*
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, Prob<eT>(s),          pi); // std::string
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, Prob<eT>(pi),         pi); // copy

	// malformed prob
	//
	std::string s2 = format_num<eT>("0.1 0.2 0.3");
	Prob<eT> pi2 = { eT(1)/10, eT(2)/10, eT(3)/10 };

	EXPECT_ANY_THROW( assert_proper(Prob<eT>(s2.c_str())); ); // char*
	EXPECT_ANY_THROW( assert_proper(Prob<eT>(s2));         ); // std::string
	EXPECT_ANY_THROW( assert_proper(Prob<eT>(pi2));        ); // Prob
}

TYPED_TEST_P(ProbTest, Uniform) {
	typedef TypeParam eT;

	Prob<eT> pi(1);
	uniform(pi);
	EXPECT_PRED_FORMAT2(prob_is_proper_size2<eT>, pi, 1);

	pi = uniform<eT>(4);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, format_num<eT>("0.25 0.25 0.25 0.25"));
}

TYPED_TEST_P(ProbTest, Randu) {
	typedef TypeParam eT;

	Prob<eT> pi(200);
	randu(pi);
	EXPECT_PRED_FORMAT2(prob_is_proper_size2<eT>, pi, 200);

	pi = probab::randu<eT>(5);
	EXPECT_PRED_FORMAT2(prob_is_proper_size2<eT>, pi, 5);
}

TYPED_TEST_P(ProbTest, Dirac) {
	typedef TypeParam eT;

	Prob<eT> pi(4);
	probab::point(pi);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, "1 0 0 0");

	pi = probab::point<eT>(4, 2);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, "0 0 1 0");
}


// run the ProbTest test-case for double, float, urat
//
REGISTER_TYPED_TEST_SUITE_P(ProbTest, Construct, Uniform, Randu, Dirac);

INSTANTIATE_TYPED_TEST_SUITE_P(Prob, ProbTest, AllTypes);

