#include "qif"

namespace qif::wrapper {


// OSQP methods /////////////////////////////////////////////////////////////////

void osqp_set_default_settings(OSQPSettings *settings)									{ return ::osqp_set_default_settings(settings); }
OSQPInt osqp_setup(OSQPSolver **solverp, const OSQPCscMatrix *P, const OSQPFloat *q, const OSQPCscMatrix *A, const OSQPFloat *l, const OSQPFloat *u, OSQPInt m, OSQPInt n, const OSQPSettings *settings)
																						{ return ::osqp_setup(solverp, P, q, A, l, u, m, n, settings); }
OSQPInt osqp_solve(OSQPSolver *solver)													{ return ::osqp_solve(solver); }
OSQPInt osqp_cleanup(OSQPSolver *solver)												{ return ::osqp_cleanup(solver); }


} // namespace qif::wrapper
