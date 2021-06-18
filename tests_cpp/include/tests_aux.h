#ifndef _QIF_tests_aux_h_
#define _QIF_tests_aux_h_

#include <string>
#include <regex>
#include "gtest/gtest.h"

#include "qif"

using namespace qif;
using std::string;


// for rat tests, change num strings like "0.5 0.1" to "5/10 1/10"
//
template<typename eT>
inline
string format_num(string s) {
	if constexpr (std::is_same<eT, rat>::value) {
		std::smatch m;
		std::regex e("(\\d+)\\.(\\d+)");

		while(std::regex_search(s, m, e)) {
			string num = m[1];
			string den = m[2];
			string r = num + den + "/1" + string(den.length(), '0');		// turn 1.5 to 15/10
			s = s.replace(m.position(), m.length(), r);
		}
	}
	return s;
}


// for tests, change num strings like "5/10 1/10" to "0.5 0.1"
//
template<typename eT>
inline
string format_rat(string s) {
	if constexpr (!std::is_same<eT, rat>::value) {
		std::smatch m;
		std::regex e("(\\d+)/(\\d+)");

		while(std::regex_search(s, m, e)) {
			string num = m[1];
			string den = m[2];
			string r = std::to_string(to_double(rat(std::stoi(num), std::stoi(den))));
			s = s.replace(m.position(), m.length(), r);
		}
	}
	return s;
}


typedef ::testing::Types<double, float, rat> AllTypes;
typedef ::testing::Types<double, float> NativeTypes;
typedef ::testing::Types<rat> RatTypes;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class BaseTest : public ::testing::Test {
	public:
		const bool is_rat = std::is_same<eT,rat>::value;

		const Prob<eT>
			unif_2   = probab::uniform<eT>(2),
			unif_4   = probab::uniform<eT>(4),
			unif_10  = probab::uniform<eT>(10),
			unif_100 = probab::uniform<eT>(100),
			point_2  = probab::point<eT>(2),
			point_4  = probab::point<eT>(4),
			point_10 = probab::point<eT>(10),
			pi1      = format_num<eT>("0.2 0.8"),
			pi2      = format_num<eT>("0.2 0.8 0 0 0 0 0 0 0 0"),
			pi3      = format_num<eT>("0.25 0.75"),
			pi4      = format_num<eT>("0.75 0.25"),
			pi5		 = format_num<eT>("0.1 0.1 0.1 0.7"),
			prand_10 = probab::randu<eT>(10);

		const Chan<eT>
			id_2      = channel::identity<eT>(2),
			id_4      = channel::identity<eT>(4),
			id_10     = channel::identity<eT>(10),
			noint_4   = channel::no_interference<eT>(4),
			noint_10  = channel::no_interference<eT>(10),
			c1        = format_num<eT>("0.8 0.2; 0.3 0.7"),
			crand_2   = channel::randu<eT>(2),
			crand_10  = channel::randu<eT>(10),
			crand_100 = channel::randu<eT>(100);
};


// for use with EXPECT_PRED_FORMATN
//
template<typename eT>
inline ::testing::AssertionResult equal2(const char* x_expr, const char* y_expr, const eT& x, const eT& y) {
	return equal<eT>(x, y)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << x_expr << " = " << x << " and " << y_expr << " = " << y << " are not equal";
}

template<typename eT>
inline ::testing::AssertionResult equal4(const char* x_expr, const char* y_expr, const char* /* md_expr */, const char* /* mrd_expr */, const eT& x, const eT& y, const eT& md, const eT& mrd) {
	return equal<eT>(x, y, md, mrd)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << x_expr << " = " << x << " and " << y_expr << " = " << y << " are not equal (md: " << md << ", mrd: " << mrd << ")";
}

// ----- channels -------
template<typename eT>
inline ::testing::AssertionResult chan_equal2(const char* x_expr, const char* y_expr, const Chan<eT>& x, const Chan<eT>& y) {
	return channel::equal<eT>(x, y)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << "channels " << x_expr << " and " << y_expr << " are not equal\n" << x_expr << ":\n" << x << y_expr << ":\n" << y;
}

template<typename eT>
inline ::testing::AssertionResult chan_equal4(const char* x_expr, const char* y_expr, const char* md_expr, const char* mrd_expr, const Chan<eT>& x, const Chan<eT>& y, const eT& md, const eT& mrd) {
	return channel::equal<eT>(x, y, md, mrd)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << "channels " << x_expr << " and " << y_expr << " are not equal with md: " << md_expr << ", mrd: " << mrd_expr << "\n" << x_expr << ":\n" << x << y_expr << ":\n" << y;
}

template<typename eT>
inline ::testing::AssertionResult chan_is_proper1(const char* x_expr, const Chan<eT>& x) {
	return channel::is_proper<eT>(x)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << x_expr << " is not a proper channel\n" << x_expr << ":\n" << x;
}

template<typename eT>
inline ::testing::AssertionResult chan_is_proper2(const char* x_expr, const char* mrd_expr, const Chan<eT>& x, const eT& mrd) {
	return channel::is_proper<eT>(x, mrd)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << x_expr << " is not a proper channel (mrd: " << mrd_expr << ")\n" << x_expr << ":\n" << x;
}

template<typename eT>
inline ::testing::AssertionResult chan_is_proper_size3(const char* x_expr, const char* n_rows_expr, const char* n_cols_expr, const Chan<eT>& x, const uint& n_rows, const uint& n_cols) {
	return channel::is_proper<eT>(x) && x.n_rows == n_rows && x.n_cols == n_cols
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << x_expr << " is not a proper " << n_rows_expr << " x " << n_cols_expr << " (" << n_rows << " x " << n_cols << ") channel\n" << x_expr << ":\n" << x;
}

// ----- probability distributions -------
template<typename eT>
inline ::testing::AssertionResult prob_equal2(const char* x_expr, const char* y_expr, const Prob<eT>& x, const Prob<eT>& y) {
	return probab::equal<eT>(x, y)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << "prob distributions " << x_expr << " and " << y_expr << " are not equal\n" << x_expr << ":" << x << "\n" << y_expr << ": " << y << "\n";
}

template<typename eT>
inline ::testing::AssertionResult prob_equal4(const char* x_expr, const char* y_expr, const char* md_expr, const char* mrd_expr, const Prob<eT>& x, const Prob<eT>& y, const eT& md, const eT& mrd) {
	return probab::equal<eT>(x, y, md, mrd)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << "prob distributions " << x_expr << " and " << y_expr << " are not equal with md: " << md_expr << ", mrd: " << mrd_expr << "\n" << x_expr << ": " << x << "\n" << y_expr << ": " << y;
}

template<typename eT>
inline ::testing::AssertionResult prob_is_proper1(const char* x_expr, const Prob<eT>& x) {
	return probab::is_proper<eT>(x)
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << x_expr << " is not a proper probability distribution\n" << x_expr << ": " << x << "\n";
}

template<typename eT>
inline ::testing::AssertionResult prob_is_proper_size2(const char* x_expr, const char* n_cols_expr, const Chan<eT>& x, const uint& n_cols) {
	return probab::is_proper<eT>(x) && x.n_cols == n_cols
		? ::testing::AssertionSuccess()
		: ::testing::AssertionFailure() << x_expr << " is not a proper probability distribution of size " << n_cols_expr << " == " << n_cols << "\n" << x_expr << ": " << x << "\n";
}


#endif
