
// Wrap dependencies of libqif, so that they are already linked in libqif.so, and the user
// only needs -lqif and not -losqp etc
//

namespace wrapper {


// OSQP methods /////////////////////////////////////////////////////////////////

using namespace osqp;

void osqp_set_default_settings(OSQPSettings *settings);
OSQPWorkspace* osqp_setup(const OSQPData *data, OSQPSettings *settings);
c_int osqp_solve(OSQPWorkspace *work);
c_int osqp_cleanup(OSQPWorkspace *work);
csc* csc_matrix(c_int m, c_int n, c_int nzmax, c_float *x, c_int *i, c_int *p);


}
