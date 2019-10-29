#include "tests_aux.h"
using namespace qif::qp;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class QuadraticProgramTest : public BaseTest<eT> {};

TYPED_TEST_CASE_P(QuadraticProgramTest);


TYPED_TEST_P(QuadraticProgramTest, Optimal) {
	typedef TypeParam eT;

	QuadraticProgram<eT> qp;

	// from https://osqp.org/docs/examples/demo.html
	qp.from_matrix(
		format_num<eT>("4 1; 1 2"),		// P
		format_num<eT>("1 1"),			// c
		format_num<eT>("1 1; 1 0; 0 1"),// A
		format_num<eT>("1 0 0"),		// l
		format_num<eT>("1 0.7 0.7") 	// u
	);

	EXPECT_TRUE(qp.solve());
	EXPECT_EQ(qp::Status::OPTIMAL, qp.status);
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1.88), qp.objective());
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, qp.solution(), format_num<eT>("0.3; 0.7"));
}

TYPED_TEST_P(QuadraticProgramTest, Infeasible) {
	typedef TypeParam eT;

	QuadraticProgram<eT> qp;

	// same demo, with infeasible constraints
	qp.from_matrix(
		format_num<eT>("4 1; 1 2"),
		format_num<eT>("1 1"),
		format_num<eT>("1 0; 1 0"),
		format_num<eT>("0 2"),
		format_num<eT>("1 3")
	);

	EXPECT_FALSE(qp.solve());
	EXPECT_EQ(qp::Status::INFEASIBLE, qp.status);
}


REGISTER_TYPED_TEST_CASE_P(QuadraticProgramTest, Optimal, Infeasible);

INSTANTIATE_TYPED_TEST_CASE_P(QuadraticProgram, QuadraticProgramTest, NativeTypes);

