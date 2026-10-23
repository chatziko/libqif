
// Wrap dependencies of libqif, so that they are already linked in libqif.so, and the user
// only needs -lqif and not -losqp etc
//

namespace wrapper {


// SoPlex methods ///////////////////////////////////////////////////////////////

// Solves the LP using SoPlex (implemented in src/soplex.cpp, so that soplex/boost headers are not exposed): exactly
// (rational mode) for rat, in floating point for double.
// Program and arguments are as in lp::LinearProgram: con_coeff is a map <col,row> => value, infinite bounds are
// +-infinity<eT>(). Returns one of lp::Status, and the solution in sol (if OPTIMAL).
std::string soplex_solve(bool maximize, const std::string& method, bool presolve, const std::string& msg_level,
	const std::vector<rat>& obj_coeff, const std::map<std::pair<uint,uint>,rat>& con_coeff,
	const std::vector<rat>& var_lb, const std::vector<rat>& var_ub,
	const std::vector<rat>& con_lb, const std::vector<rat>& con_ub,
	std::vector<rat>& sol);
std::string soplex_solve(bool maximize, const std::string& method, bool presolve, const std::string& msg_level,
	const std::vector<double>& obj_coeff, const std::map<std::pair<uint,uint>,double>& con_coeff,
	const std::vector<double>& var_lb, const std::vector<double>& var_ub,
	const std::vector<double>& con_lb, const std::vector<double>& con_ub,
	std::vector<double>& sol);


// OSQP methods /////////////////////////////////////////////////////////////////

void osqp_set_default_settings(OSQPSettings *settings);
OSQPInt osqp_setup(OSQPSolver **solverp, const OSQPCscMatrix *P, const OSQPFloat *q, const OSQPCscMatrix *A, const OSQPFloat *l, const OSQPFloat *u, OSQPInt m, OSQPInt n, const OSQPSettings *settings);
OSQPInt osqp_solve(OSQPSolver *solver);
OSQPInt osqp_cleanup(OSQPSolver *solver);


}
