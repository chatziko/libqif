
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


// OSQP methods /////////////////////////////////////////////////////////////////

void osqp_set_default_settings(OSQPSettings *settings);
c_int osqp_setup(OSQPWorkspace **work, const OSQPData *data, OSQPSettings *settings);
c_int osqp_solve(OSQPWorkspace *work);
c_int osqp_cleanup(OSQPWorkspace *work);
csc* csc_matrix(c_int m, c_int n, c_int nzmax, c_float *x, c_int *i, c_int *p);


}
