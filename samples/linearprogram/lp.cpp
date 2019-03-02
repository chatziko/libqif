#include "qif"

#include <iostream>
#include <list>

using namespace qif::lp;
using namespace std;

typedef rat eT;
typedef qif::MatrixEntry<eT> ME;

int main() {


	/*
	maximize –2x + 5y, subject to: 
	x <= 200 
	x >= 100
	y <= 170 
	y >= 80
	x + y >= 200 
	Solution 650 given for (x,y) = (100,170)
	*/

	// LinearProrgam solves the linear program
	// {min/max} dot(c,x)
	// subject to A x {>=|==|<=} b
	// x >= 0
	//
	LinearProgram<eT> lp;
	lp.maximize = true;

	auto vars = lp.make_vars(2);
	// cost function
	lp.set_obj_coeff(vars[0], -2);
	lp.set_obj_coeff(vars[1], 5);

	eT inf = qif::infinity<eT>();

	// coefficients
	vector<uint> cons;
	cons.push_back(lp.make_con(-inf, 200));
	cons.push_back(lp.make_con(100, inf));
	cons.push_back(lp.make_con(-inf, 170));
	cons.push_back(lp.make_con(80, inf));
	cons.push_back(lp.make_con(200, inf));

	// matrix A is:
	// 1 0
	// 1 0
	// 0 1
	// 0 1
	// 1 1
	// batch insert (non-zero entries):
	// std::list<ME> entries;
	// entries.push_back(ME(0, 0, 1));	// A(0,0) = 1
	// entries.push_back(ME(1, 0, 1));
	// entries.push_back(ME(2, 1, 1));
	// entries.push_back(ME(3, 1, 1));
	// entries.push_back(ME(4, 0, 1));
	// entries.push_back(ME(4, 1, 1));
	// qif::fill_spmat(lp.A, 5, 2, entries);

	lp.set_con_coeff(cons[0], vars[0], 1);
	lp.set_con_coeff(cons[1], vars[0], 1);
	lp.set_con_coeff(cons[2], vars[1], 1);
	lp.set_con_coeff(cons[3], vars[1], 1);
	lp.set_con_coeff(cons[4], vars[0], 1);
	lp.set_con_coeff(cons[4], vars[1], 1);

	// same thing but much slower
//	lp.A.set_size(5, 2);
//	lp.A(0,0) = 1;
//	lp.A(1,0) = 1;
//	lp.A(2,1) = 1;
//	lp.A(3,1) = 1;
//	lp.A(4,0) = 1;
//	lp.A(4,1) = 1;

	// methods: simplex_primal, simplex_dual, interior
	lp.method = method_t::simplex_primal;
	lp.glp_msg_level = msg_level_t::on;

	bool solved = lp.solve();

	// cout << lp.A;
	cout << lp.status;
	cout
		<< "solved: " << solved
		<< "\nstatus: " << lp.status
		<< "\nmethod: " << lp.method
		<< "\nsolution:\n" << lp.x
		<< "\noptimum: " << lp.optimum()
		<< "\n";

	arma::mat m = "";
	cout << arma::cdot(m,m);
}

