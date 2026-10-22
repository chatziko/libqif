#include "qif"

namespace qif::wrapper {


// GLPK methods /////////////////////////////////////////////////////////////////

#ifdef QIF_USE_GLPK
glp_prob *glp_create_prob(void)																	{ return ::glp_create_prob(); }
int glp_add_rows(glp_prob *P, int nrs)															{ return ::glp_add_rows(P, nrs); }
int glp_add_cols(glp_prob *P, int ncs)															{ return ::glp_add_cols(P, ncs); }
void glp_set_obj_dir(glp_prob *P, int dir)														{ return ::glp_set_obj_dir(P, dir); }
void glp_set_obj_coef(glp_prob *P, int j, double coef)											{ return ::glp_set_obj_coef(P, j, coef); }
void glp_set_row_bnds(glp_prob *P, int j, int type, double lb, double ub)						{ return ::glp_set_row_bnds(P, j, type, lb, ub); }
void glp_set_col_bnds(glp_prob *P, int j, int type, double lb, double ub)						{ return ::glp_set_col_bnds(P, j, type, lb, ub); }
void glp_load_matrix(glp_prob *P, int ne, const int ia[], const int ja[], const double ar[])	{ return ::glp_load_matrix(P, ne, ia, ja, ar); }
void glp_init_smcp(glp_smcp *parm)																{ return ::glp_init_smcp(parm); }
int glp_simplex(glp_prob *P, const glp_smcp *parm)												{ return ::glp_simplex(P, parm); }
int glp_get_status(glp_prob *P)																	{ return ::glp_get_status(P); }
int glp_get_dual_stat(glp_prob *P)																{ return ::glp_get_dual_stat(P); }
double glp_get_col_prim(glp_prob *P, int j)														{ return ::glp_get_col_prim(P, j); }
int glp_interior(glp_prob *P, const glp_iptcp *parm)											{ return ::glp_interior(P, parm); }
void glp_init_iptcp(glp_iptcp *parm)															{ return ::glp_init_iptcp(parm); }
int glp_ipt_status(glp_prob *P)																	{ return ::glp_ipt_status(P); }
double glp_ipt_col_prim(glp_prob *P, int j)														{ return ::glp_ipt_col_prim(P, j); }
void glp_delete_prob(glp_prob *P)																{ return ::glp_delete_prob(P); }
int glp_free_env(void)																			{ return ::glp_free_env(); }
#endif


// OSQP methods /////////////////////////////////////////////////////////////////

void osqp_set_default_settings(OSQPSettings *settings)									{ return ::osqp_set_default_settings(settings); }
OSQPInt osqp_setup(OSQPSolver **solverp, const OSQPCscMatrix *P, const OSQPFloat *q, const OSQPCscMatrix *A, const OSQPFloat *l, const OSQPFloat *u, OSQPInt m, OSQPInt n, const OSQPSettings *settings)
																						{ return ::osqp_setup(solverp, P, q, A, l, u, m, n, settings); }
OSQPInt osqp_solve(OSQPSolver *solver)													{ return ::osqp_solve(solver); }
OSQPInt osqp_cleanup(OSQPSolver *solver)												{ return ::osqp_cleanup(solver); }


} // namespace qif::wrapper
