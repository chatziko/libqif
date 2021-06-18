#include "tests_aux.h"
using namespace qif::lp;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class LinearProgramTest : public BaseTest<eT> {
	public:
		void SetUp() override {
			// construct method/solver/presolve combinations to test
			std::list<string> solvers = { Solver::INTERNAL };
			#ifdef QIF_USE_ORTOOLS
				solvers.push_back(Solver::GLOP);
				solvers.push_back(Solver::CLP);
			#else
				std::cerr << "\nOR-tools not found, skipping GLOP and CLP tests\n\n";
			#endif
			#ifdef QIF_USE_GLPK
				solvers.push_back(Solver::GLPK);
			#else
				std::cerr << "\nGLPK not found, skipping GLPK tests\n\n";
			#endif

			for(string method : { Method::SIMPLEX_PRIMAL, Method::SIMPLEX_DUAL, Method::INTERIOR }) {
			for(string solver : solvers) {
			for(bool presolve : { false, true }) {
				// some combinations are not valid
				if(method == Method::INTERIOR && (presolve || solver == Solver::GLOP || solver == Solver::CLP)) continue; // interior: no presolver, no GLOP support, unstable with CLP
				if(solver == Solver::INTERNAL && (presolve || method != Method::SIMPLEX_PRIMAL)               ) continue; // internal solver: only simplex_primal/no presolve
				if(this->is_rat               && solver != Solver::INTERNAL                                   ) continue; // rat: only internal solver

				combs.push_back(std::tuple(method, solver, presolve));
			}}}
		}
		std::vector<std::tuple<string,string,bool>> combs;
};

TYPED_TEST_SUITE_P(LinearProgramTest);


TYPED_TEST_P(LinearProgramTest, Optimal) {
	typedef TypeParam eT;
	LinearProgramTest<eT>& t = *this;

	// the default acceptance range is too string for linear programs, we need a more permissive mrd
	eT md(def_md<eT>);
	eT mrd(def_mrd<float>);		// always use the mrd for floats

	for(auto comb : t.combs) {

		LinearProgram<eT> lp;
		std::tie(lp.method, lp.solver, lp.presolve) = comb;

		lp.from_matrix(
			format_num<eT>("1 2; 3 1"),
			format_num<eT>("1 2"),
			format_num<eT>("0.6 0.5")
		);

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(Status::OPTIMAL, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(46)/100, lp.objective(), md, mrd);
		EXPECT_PRED_FORMAT4(chan_equal4<eT>, lp.solution(), format_num<eT>("0.6; 0.2"), md, mrd);

		lp.from_matrix(
			format_num<eT>("1 1 0; 0 1 1"),
			format_num<eT>("1 1"),
			format_num<eT>("1 2 -1")
		);

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(Status::OPTIMAL, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(2), lp.objective(), md, mrd);
		EXPECT_PRED_FORMAT4(chan_equal4<eT>, lp.solution(), format_num<eT>("0; 1; 0"), md, mrd);

		lp.maximize = false;
		lp.from_matrix(
			format_num<eT>("3 -4; 1 2; 1 0"),
			format_num<eT>("12 4 1"),
			format_num<eT>("3 4"),
			{ '<', '>', '>' }
		);

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(Status::OPTIMAL, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(9), lp.objective(), md, mrd);
		EXPECT_PRED_FORMAT4(chan_equal4<eT>, lp.solution(), format_num<eT>("1; 1.5"), md, mrd);

		lp.maximize = false;
		lp.from_matrix(
			format_num<eT>("1 2 2; 2 1 2; 2 2 1"),
			format_num<eT>("20 20 20"),
			format_num<eT>("-10 -12 -12")
		);

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(Status::OPTIMAL, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(-136), lp.objective(), md, mrd);
		EXPECT_PRED_FORMAT4(chan_equal4<eT>, lp.solution(), format_num<eT>("4; 4; 4"), md, mrd);

		lp.clear();
		lp.maximize = false;
		auto v = lp.make_var(eT(-5), infinity<eT>());
		lp.make_con(eT(0), eT(0));		// glpk needs at least one constraint, add dummy one
		lp.set_obj_coeff(v, eT(1));

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(Status::OPTIMAL, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(-5), lp.objective(), md, mrd);
		EXPECT_PRED_FORMAT4(chan_equal4<eT>, lp.solution(), format_num<eT>("-5"), md, mrd);

		// the problematic program from the add. refin. metric: redundant constraints causing internal simplex to fail
		lp.maximize = false;
		lp.from_matrix(
			format_num<eT>(R"(
				-1       0        1        0        0        0       -1        0        0        0        0        0        0        0        0;
				-1       0        0        0        1        0        0       -1        0        0        0        0        0        0        0;
				0.75     0    -0.75        0        0        0        0        0       -1        0        0        0        0        0        0;
				0        0    -0.75        0     0.75        0        0        0        0       -1        0        0        0        0        0;
				0.25     0        0        0    -0.25        0        0        0        0        0       -1        0        0        0        0;
				0        0     0.25        0    -0.25        0        0        0        0        0        0       -1        0        0        0;
				1        0        0        0        0        0        0        0        0        0        0        0        1        0        0;
				0        0     0.75        0        0        0        0        0        0        0        0        0        0        1        0;
				0        0        0        0     0.25        0        0        0        0        0        0        0        0        0        1;
			)"),
			format_num<eT>(" 0 0 0 0 0 0 1 0.75 0.25"),
			format_num<eT>("-1 0 0.75 0 0.25 0 0 0 0 0 0 0 0 0 0"),
			"= = = = = = = = ="
		);

		EXPECT_TRUE(lp.solve());
		EXPECT_EQ(Status::OPTIMAL, lp.status);
		EXPECT_PRED_FORMAT4(equal4<eT>, eT(0), lp.objective(), md, mrd);
	}
}

TYPED_TEST_P(LinearProgramTest, Infeasible) {
	typedef TypeParam eT;
	LinearProgramTest<eT>& t = *this;

	for(auto comb : t.combs) {

		LinearProgram<eT> lp;
		std::tie(lp.method, lp.solver, lp.presolve) = comb;

		string status = lp.method == Method::INTERIOR
				? Status::INFEASIBLE_OR_UNBOUNDED	// sometimes we just know that the problem is infeasible OR unbounded
				: Status::INFEASIBLE;

		lp.from_matrix(
			format_num<eT>("1; 1"),
			format_num<eT>("3 2"),
			format_num<eT>("1"),
			{ '>', '<' }
		);

		EXPECT_FALSE(lp.solve());
		EXPECT_EQ(status, lp.status);

		lp.from_matrix(
			format_num<eT>("1; -1"),
			format_num<eT>("3 -2"),
			format_num<eT>("4"),
			{ '>', '>' }
		);

		EXPECT_FALSE(lp.solve());
		EXPECT_EQ(status, lp.status);
	}
}

TYPED_TEST_P(LinearProgramTest, Unbounded) {
	typedef TypeParam eT;
	LinearProgramTest<eT>& t = *this;

	for(auto comb : t.combs) {

		LinearProgram<eT> lp;
		std::tie(lp.method, lp.solver, lp.presolve) = comb;

		// EXTRA conditions only for unbounded
		if(lp.solver == Solver::GLOP || lp.solver == Solver::CLP) continue; // OR-tools/DUAL seems unstable with unbounded problems (TODO: investigae)

		auto status = lp.solver != Solver::GLPK || (lp.method == Method::SIMPLEX_PRIMAL && !lp.presolve)
				? Status::UNBOUNDED
				: Status::INFEASIBLE_OR_UNBOUNDED;	// sometimes we just know that the problem is infeasible OR unbounded

		lp.maximize = false;
		lp.from_matrix(
			format_num<eT>("1"),
			format_num<eT>("2"),
			format_num<eT>("-1"),
			{ '>' }
		);

		EXPECT_FALSE(lp.solve());
		EXPECT_EQ(status, lp.status);
	}
}


REGISTER_TYPED_TEST_SUITE_P(LinearProgramTest, Optimal, Infeasible, Unbounded);

INSTANTIATE_TYPED_TEST_SUITE_P(LinearProgram, LinearProgramTest, AllTypes);

