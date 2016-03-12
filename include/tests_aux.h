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
	return s;
}
template<>
inline
string format_num<rat>(string s) {
	std::smatch m;
	std::regex e("(\\d+)\\.(\\d+)");

	while(std::regex_search(s, m, e)) {
		string num = m[1];
		string den = m[2];
		string r = num + den + "/1" + string(den.length(), '0');		// turn 1.5 to 15/10
		s = s.replace(m.position(), m.length(), r);
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
		Prob<eT>
			unif_2   = probab::uniform<eT>(2),
			unif_4   = probab::uniform<eT>(4),
			unif_10  = probab::uniform<eT>(10),
			unif_100 = probab::uniform<eT>(100),
			dirac_2  = probab::dirac<eT>(2),
			dirac_4  = probab::dirac<eT>(4),
			dirac_10 = probab::dirac<eT>(10),
			pi1      = format_num<eT>("0.2 0.8"),
			pi2      = format_num<eT>("0.2 0.8 0 0 0 0 0 0 0 0"),
			pi3      = format_num<eT>("0.25 0.75"),
			pi4      = format_num<eT>("0.75 0.25"),
			pi5		 = format_num<eT>("0.1 0.1 0.1 0.7");

		Chan<eT>
			id_2      = channel::identity<eT>(2),
			id_4      = channel::identity<eT>(4),
			id_10     = channel::identity<eT>(10),
			noint_4   = channel::no_interference<eT>(4),
			noint_10  = channel::no_interference<eT>(10),
			c1        = format_num<eT>("0.8 0.2; 0.3 0.7"),
			crand_10  = channel::randu<eT>(10),
			crand_100 = channel::randu<eT>(100);
};


// this is used with EXPECT_PRED2, cause 'equal' takes 4 arguments
//
template<typename eT>
inline bool equal2(const eT& x, const eT& y) {
	return equal<eT>(x, y);
}

template<typename eT>
inline bool chan_equal2(const Chan<eT>& x, const Chan<eT>& y) {
	return channel::equal<eT>(x, y);
}

template<typename eT>
inline bool prob_is_proper1(const Prob<eT>& x) {
	return probab::is_proper<eT>(x);
}

template<typename eT>
inline bool chan_is_proper1(const Chan<eT>& x) {
	return channel::is_proper<eT>(x);
}


template<typename eT>
void expect_channel(const Mat<eT>& m, const Chan<eT>& c) {
	EXPECT_PRED2(chan_equal2<eT>, m, c);
	EXPECT_PRED1(chan_is_proper1<eT>, c);
}

template<typename eT>
void expect_channel(const string& s, const Chan<eT>& c) {
	expect_channel(Mat<eT>(s), c);
}

template<typename eT>
void expect_channel(uint rn, uint cn, const Chan<eT>& c) {
	EXPECT_EQ(rn, c.n_rows);
	EXPECT_EQ(cn, c.n_cols);

	EXPECT_PRED1(chan_is_proper1<eT>, c);
}


template<typename eT>
void expect_prob(const Prob<eT>& m, const Prob<eT>& p) {
	EXPECT_EQ(m.n_cols, p.n_cols);

	for(uint j = 0; j < p.n_cols; j++)
		EXPECT_PRED2(equal2<eT>, m.at(j), p.at(j));

	EXPECT_PRED1(prob_is_proper1<eT>, p);
}

template<typename eT>
void expect_prob(const string& s, const Prob<eT>& p) {
	expect_prob(Prob<eT>(s), p);
}

template<typename eT>
void expect_prob(uint cn, const Prob<eT>& p) {
	EXPECT_EQ(cn, p.n_cols);

	EXPECT_PRED1(prob_is_proper1<eT>, p);
}


template<typename eT>
void expect_mat(const Mat<eT>& m, const Mat<eT>& c, const eT& md = def_max_diff<eT>(), const eT& mrd = def_max_rel_diff<eT>()) {
	EXPECT_EQ(m.n_rows, c.n_rows);
	EXPECT_EQ(m.n_rows, c.n_rows);

	for(uint i = 0; i < c.n_rows; i++)
		for(uint j = 0; j < c.n_cols; j++)
			EXPECT_PRED4(equal<eT>, m.at(i, j), c.at(i, j), md, mrd);
}

template<typename eT>
void expect_mat(const string& s, const Mat<eT>& c, const eT& md = def_max_diff<eT>(), const eT& mrd = def_max_rel_diff<eT>()) {
	expect_mat(Mat<eT>(s), c, md, mrd);
}

template<typename eT>
void expect_mat(uint rn, uint cn, const Mat<eT>& c) {
	EXPECT_EQ(rn, c.n_rows);
	EXPECT_EQ(cn, c.n_cols);
}

#endif
