#include "qif"

#include <iostream>
#include <list>

using namespace qif::lp;
using namespace std;

typedef double eT;

int main() {
	// force link with cblas (TODO: investigate)
	arma::mat m = "";
	eT d = arma::cdot(m,m);
	d=d;


	// LinearProrgam solves the linear program
	// {min/max} dot(c,x)
	// subject to A x {>=|==|<=} b
	// x >= 0
	//
	eT inf = qif::infinity<eT>();
	LinearProgram<eT> lp;


	// auto v = lp.make_var(-inf, inf);
	// lp.set_obj_coeff(v, 1);
	// auto c = lp.make_con(-5, 5);
	// lp.set_con_coeff(c, v, 1);
	// lp.maximize = false;
	// bool s = lp.solve();
	// cout
	// 	<< "solved: " << s << "\n"
	// 	<< "\nstatus: " << lp.status
	// 	<< "\nmethod: " << lp.method
	// 	<< "\nsolution:\n" << lp.solution()
	// 	<< "\nobjective: " << lp.objective()
	// 	<< "\n";
	// return 0;




	/*
	maximize –2x + 5y, subject to: 
	x <= 200 
	x >= 100
	y <= 170 
	y >= 80
	x + y >= 200 
	Solution 650 given for (x,y) = (100,170)
	*/
	lp.maximize = true;

	auto x = lp.make_var(100, 200);
	auto y = lp.make_var(80, 170);
	// cost function
	lp.set_obj_coeff(x, -2);
	lp.set_obj_coeff(y, 5);

	auto con = lp.make_con(200, inf);
	lp.set_con_coeff(con, x, 1);
	lp.set_con_coeff(con, y, 1);

	// methods: simplex_primal, simplex_dual, interior
	lp.method = Method::SIMPLEX_DUAL;
	lp.msg_level = MsgLevel::ON;

	lp.solver = Solver::GLOP;
	bool solved = lp.solve();
	cout
		<< "solved: " << solved
		<< "\nstatus: " << lp.status
		<< "\nmethod: " << lp.method
		<< "\nsolver: " << lp.solver
		<< "\nsolution:\n" << lp.solution()
		<< "\nobjective: " << lp.objective()
		<< "\n";
}

