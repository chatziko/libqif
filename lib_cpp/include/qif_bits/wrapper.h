
// Wrap dependencies of libqif, so that they are already linked in libqif.so, and the user
// only needs -lqif and not -losqp etc
//

namespace wrapper {


// GLPK methods /////////////////////////////////////////////////////////////////

#ifdef QIF_USE_GLPK
glp_prob *glp_create_prob(void);
int glp_add_rows(glp_prob *P, int nrs);
int glp_add_cols(glp_prob *P, int ncs);
void glp_set_obj_dir(glp_prob *P, int dir);
void glp_set_obj_coef(glp_prob *P, int j, double coef);
void glp_set_row_bnds(glp_prob *P, int j, int type, double lb, double ub);
void glp_load_matrix(glp_prob *P, int ne, const int ia[], const int ja[], const double ar[]);
void glp_set_col_bnds(glp_prob *P, int j, int type, double lb, double ub);
void glp_init_smcp(glp_smcp *parm);
int glp_simplex(glp_prob *P, const glp_smcp *parm);
int glp_get_status(glp_prob *P);
int glp_get_dual_stat(glp_prob *P);
double glp_get_col_prim(glp_prob *P, int j);
int glp_interior(glp_prob *P, const glp_iptcp *parm);
void glp_init_iptcp(glp_iptcp *parm);
int glp_ipt_status(glp_prob *P);
double glp_ipt_col_prim(glp_prob *P, int j);
void glp_delete_prob(glp_prob *P);
int glp_free_env(void);
#endif


// SoPlex methods ///////////////////////////////////////////////////////////////

// Solves the LP exactly using SoPlex (implemented in src/soplex.cpp, so that soplex/boost headers are not exposed).
// Program and arguments are as in lp::LinearProgram: con_coeff is a map <col,row> => value, infinite bounds are
// +-infinity<rat>(). Returns one of lp::Status, and the solution in sol (if OPTIMAL).
std::string soplex_solve(bool maximize, const std::string& method, bool presolve, const std::string& msg_level,
	const std::vector<rat>& obj_coeff, const std::map<std::pair<uint,uint>,rat>& con_coeff,
	const std::vector<rat>& var_lb, const std::vector<rat>& var_ub,
	const std::vector<rat>& con_lb, const std::vector<rat>& con_ub,
	std::vector<rat>& sol);


// OSQP methods /////////////////////////////////////////////////////////////////

void osqp_set_default_settings(OSQPSettings *settings);
OSQPInt osqp_setup(OSQPSolver **solverp, const OSQPCscMatrix *P, const OSQPFloat *q, const OSQPCscMatrix *A, const OSQPFloat *l, const OSQPFloat *u, OSQPInt m, OSQPInt n, const OSQPSettings *settings);
OSQPInt osqp_solve(OSQPSolver *solver);
OSQPInt osqp_cleanup(OSQPSolver *solver);


}
