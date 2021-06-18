/*
The simplex method is adapted from
https://github.com/IainNZ/RationalSimplex.jl
Original code is under the MIT licence.

The MIT License (MIT)
Copyright (c) 2014 Iain Dunning
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace lp {

using std::string;

// Use a map <col, row> => value for sparse arrays. Put column first so that
// the lexicographic order used by map is top-to-bottom/left-to-right
// (this is for uniformity with QuadraticProgram, not really needed here)
template<typename eT> using Sparse = std::map<std::pair<uint,uint>,eT>;

// emulate enums with strings, it's easier
#undef ERROR	// MSVC adds this
namespace Status { const auto OPTIMAL = "OPTIMAL", INFEASIBLE = "INFEASIBLE", UNBOUNDED = "UNBOUNDED", INFEASIBLE_OR_UNBOUNDED = "INFEASIBLE_OR_UNBOUNDED", ERROR = "ERROR"; }
namespace Method { const auto AUTO = "AUTO", SIMPLEX_PRIMAL = "SIMPLEX_PRIMAL", SIMPLEX_DUAL = "SIMPLEX_DUAL", INTERIOR = "INTERIOR"; }					// AUTO: whatever is best
namespace Solver { const auto AUTO = "AUTO", INTERNAL = "INTERNAL", GLPK = "GLPK", GLOP = "GLOP", CLP = "CLP", GUROBI = "GUROBI", CPLEX = "CPLEX"; }	// for the each application
namespace MsgLevel { const auto OFF = "OFF", ERR = "ERR", ON = "ON", ALL = "ALL"; }


class Defaults {
	public:
		static bool presolve;
		static string msg_level;
		static string method;
		static string solver;
};

// Solve the linear program
// {min/max} dot(c,x)
// subject to lb <= A x <= ub
// lb <= x <= ub

template<typename eT>
class LinearProgram {
	typedef uint Var;
	typedef uint Con;

	const bool is_rat = std::is_same<eT, rat>::value;

	public:
		bool maximize = true;
		string method = Defaults::method;
		string solver = Defaults::solver;
		bool presolve = Defaults::presolve;
		string status;
		string msg_level = Defaults::msg_level;

		bool solve();
		string to_mps();

		eT objective();
		inline eT solution(Var x)	{ return sol(x); }
		inline Col<eT> solution()	{ return sol; };

		void clear();
		void from_matrix(const arma::SpMat<eT>& A, const Col<eT>& b, const Col<eT>& c, const Col<char>& sense = {}, bool non_negative = true);

		// API for creating variables and constraints
		Var make_var(eT lb = -infinity<eT>(), eT ub = infinity<eT>());
		std::vector<Var> make_vars(uint n, eT lb = -infinity<eT>(), eT ub = infinity<eT>());
		std::vector<std::vector<Var>> make_vars(uint n, uint m, eT lb = -infinity<eT>(), eT ub = infinity<eT>());

		Con make_con(eT lb, eT ub);

		void set_obj_coeff(Var var, eT coeff, bool add = false);
		void set_con_coeff(Con cons, Var var, eT coeff, bool add = false);

	protected:
		Col<eT> sol;			// solution

		bool glpk();
		bool ortools();
		bool internal_solver();
		bool simplex();

		void to_canonical_form();
		Col<eT> original_solution();
		void dump();

		uint n_var = 0,						// number of variables
			 n_con = 0;						// number of constraints

		std::vector<eT> obj_coeff;			// coefficients for the objective function
		Sparse<eT> con_coeff;				// coefficients for the constraints: <row,col> => val
		std::vector<eT> var_lb, var_ub,		// variables lower/upper bounds
						con_lb, con_ub;		// constraints lower/upper

		// info for transforming solution to the original one (see canonical_form, original_solution)
		std::vector< std::tuple<int,eT,eT> > var_transform;
};


template<typename eT>
inline
eT LinearProgram<eT>::objective() {
	if(sol.empty())
		throw std::runtime_error("no solution");

	eT res(0);
	for(uint var = 0; var < n_var; var++)
		res += obj_coeff[var] * sol(var);
	return res;
}

template<typename eT>
Col<eT> LinearProgram<eT>::original_solution() {
	Col<eT> res(var_transform.size());

	for(uint v = 0; v < res.n_elem; v++) {
		// A tuple (var, coeff, add) for v means that the original value of v is:  v*coeff + add - var
		auto& [v2, coeff, add] = var_transform[v];

		res(v) = sol(v) * coeff + add - (v2 >= 0 ? sol(v2) : eT(0));
	}

	return res;
}

template<typename eT>
inline
typename LinearProgram<eT>::Var LinearProgram<eT>::make_var(eT lb, eT ub) {
	var_lb.push_back(lb);
	var_ub.push_back(ub);
	obj_coeff.push_back(eT(0));
	return n_var++;
}

template<typename eT>
inline
std::vector<typename LinearProgram<eT>::Var> LinearProgram<eT>::make_vars(uint n, eT lb, eT ub) {
	std::vector<Var> res(n);
	for(uint i = 0; i < n; i++)
		res[i] = make_var(lb, ub);
	return res;
}

template<typename eT>
inline
std::vector<std::vector<typename LinearProgram<eT>::Var>> LinearProgram<eT>::make_vars(uint n, uint m, eT lb, eT ub) {
	std::vector<std::vector<Var>> res(n);
	for(uint i = 0; i < n; i++) {
		res[i].resize(m);
		for(uint j = 0; j < m; j++)
			res[i][j] = make_var(lb, ub);
	}
	return res;
}

template<typename eT>
inline
typename LinearProgram<eT>::Con LinearProgram<eT>::make_con(eT lb, eT ub) {
	if(ub == infinity<eT>() && lb == -ub)
		throw std::runtime_error("trying to add unconstrained constraint");

	con_lb.push_back(lb);
	con_ub.push_back(ub);
	return n_con++;
}

template<typename eT>
inline
void LinearProgram<eT>::set_obj_coeff(Var var, eT coeff, bool add) {
	if(add)
		obj_coeff[var] += coeff;
	else
		obj_coeff[var] = coeff;
}

template<typename eT>
inline
void LinearProgram<eT>::set_con_coeff(Con con, Var var, eT coeff, bool add) {
	if(equal<eT>(coeff, eT(0)))
		return;

	auto key = std::pair(var, con);	// <col, row>
	if(add && con_coeff.count(key))
		con_coeff[key] += coeff;
	else
		con_coeff[key] = coeff;
}


template<typename eT>
void LinearProgram<eT>::clear() {
	obj_coeff.clear();
	con_coeff.clear();
	var_lb.clear();
	var_ub.clear();
	con_lb.clear();
	con_ub.clear();
	n_var = n_con = 0;
}

// input program in matrix form
// A: constraint coefficients
// b: constraint constants
// c: objective constants
// sense: < = > for each constraint
// non_negative: if true assume x >= 0
//
template<typename eT>
void LinearProgram<eT>::from_matrix(const arma::SpMat<eT>& A, const Col<eT>& b, const Col<eT>& c, const Col<char>& sense, bool non_negative) {

	if(A.n_rows != b.n_rows || A.n_cols != c.n_rows)
		throw std::runtime_error("invalid size");

	clear();

	n_var = c.n_elem;
	n_con = A.n_rows;

	var_lb.resize(n_var, non_negative ? eT(0) : -infinity<eT>());
	var_ub.resize(n_var, infinity<eT>());

	for(uint con = 0; con < n_con; con++) {
		char s = sense.n_elem > con ? sense(con) : '<';
		con_lb.push_back(s == '<' ? -infinity<eT>() : b[con]);
		con_ub.push_back(s == '>' ?  infinity<eT>() : b[con]);
	}

	for(uint var = 0; var < n_var; var++) {
		obj_coeff.push_back(c[var]);
	}
	
	auto end = A.end();
	for(auto c = A.begin(); c != end; ++c) {		// c++ throws weird warning, ++c doesn't!
		set_con_coeff(c.row(), c.col(), *c);
	}
}

template<typename eT>
bool LinearProgram<eT>::solve() {
	// AUTO: internal for rat, GLPK for Interior, CLP if available, otherwise GLPK
	auto s = solver;
	if(s == Solver::AUTO) {
		if(is_rat) {
			s = Solver::INTERNAL;
		} else if(method == Method::INTERIOR) {
			#ifdef QIF_USE_GLPK
			s = Solver::GLPK;
			#else
			s = Solver::INTERNAL;
			#endif
		} else {
			#ifdef QIF_USE_ORTOOLS
			s = Solver::GLOP;
			#elif defined(QIF_USE_GLPK)
			s = Solver::GLPK;
			#else
			s = Solver::INTERNAL;
			#endif
		}
	}
	
	if(msg_level != MsgLevel::OFF)
		std::cerr << "Solving LP with solver: " << s << "\n";

	return
		s == Solver::GLPK ? glpk() :
		s == Solver::INTERNAL ? internal_solver() :
		ortools();		// make sure that AUTO in ortools() is treated in the same way as here!
}

// internal solver, uses simplex() after cloning and converting to canonical form. Mostly useful for rats
//
template<typename eT>
inline
bool LinearProgram<eT>::internal_solver() {

	if(method != Method::AUTO && method != Method::SIMPLEX_PRIMAL)
		throw std::runtime_error("internal solver only supports simplex_primal method");

	LinearProgram<eT> lp(*this);	// clone
	lp.to_canonical_form();

	bool res = lp.simplex();
	status = lp.status;

	if(res)
		sol = lp.original_solution();

	return res;
}

// solve the program using GLPK (directly, not through ortools)
template<typename eT>
bool LinearProgram<eT>::glpk() {
	if(is_rat)
		std::cout << "\nWARNING: using GLPK with rat. This will convert to double so it's not exact.\n\n";

#ifndef QIF_USE_GLPK
	throw std::runtime_error("glpk not available");
#else

	// create problem
	glp_prob *lp = wrapper::glp_create_prob();

	wrapper::glp_set_obj_dir(lp, maximize ? GLP_MAX : GLP_MIN);

	eT  inf = infinity<eT>();
	eT minf = -inf;

	// add variables
	// CAREFULL: all glp indexes are 1-based
	//
	wrapper::glp_add_cols(lp, n_var);
	for(uint j = 0; j < n_var; j++) {
		eT lb = var_lb[j];
		eT ub = var_ub[j];

		int type =
			lb == minf && ub == inf ? GLP_FR :	// free
			lb == minf              ? GLP_UP :	// upper only
			ub ==  inf              ? GLP_LO :	// lower only
			lb ==   ub              ? GLP_FX :	// fixed value
			                          GLP_DB ;	// both bounds

		wrapper::glp_set_col_bnds(lp, j+1, type, to_double(lb), to_double(ub));

		wrapper::glp_set_obj_coef(lp, j+1, to_double(obj_coeff[j]));				// coefficient in the cost functoin
	}

	// add constraints. glpk uses a "sparse" way of entering the rows, using
	// flat arrays ia, ja, ar. ia[z],ja[z] are the indexes of A for the value to
	// set, and ar[z] = A[ ia[z], ja[z] ]
	// we add entries only for non-zero elements, it's much faster!
	//
	wrapper::glp_add_rows(lp, n_con);

	int size = con_coeff.size();

	std::vector<int>	ia(size+1),
						ja(size+1);
	std::vector<double> ar(size+1);

	for(uint i = 0; i < n_con; i++) {
		eT lb = con_lb[i];
		eT ub = con_ub[i];

		int type =
			lb == minf && ub == inf ? GLP_FR :	// free
			lb == minf              ? GLP_UP :	// upper only
			ub ==  inf              ? GLP_LO :	// lower only
			lb ==   ub              ? GLP_FX :	// fixed value
			                          GLP_DB ;	// both bounds

		wrapper::glp_set_row_bnds(lp, i+1, type, to_double(lb), to_double(ub));
	}

	// loop over non-zero elements of sparse array
	int index = 1;
	for(auto& [key, val] : con_coeff) {
		auto& [col, row] = key;
		ia[index] = row + 1;
		ja[index] = col + 1;
		ar[index] = to_double(val);
		index++;
	}

	wrapper::glp_load_matrix(lp, size, &ia[0], &ja[0], &ar[0]);

	std::map<string,int> glp_msg_levs = { { MsgLevel::OFF, GLP_MSG_OFF }, { MsgLevel::ERR, GLP_MSG_ERR }, { MsgLevel::ON, GLP_MSG_ON }, { MsgLevel::ALL, GLP_MSG_ALL } };
	const int msg_lev = glp_msg_levs[msg_level];

	// solve
	const bool is_interior = method == Method::INTERIOR;
	if(!is_interior) {	// simplex primal/dual
		glp_smcp opt;
		wrapper::glp_init_smcp(&opt);
		opt.meth = method == Method::SIMPLEX_PRIMAL ? GLP_PRIMAL : GLP_DUALP;	// DUALP: use dual, switch to primal if it fails. DUALP is also set if method == AUTO
		opt.msg_lev = msg_lev;							// debug info sent to terminal, default off
		opt.presolve = presolve ? GLP_ON : GLP_OFF;		// use presolver

		//glp_scale_prob(lp, GLP_SF_AUTO);	// scaling is done by the presolver
		int glp_res = wrapper::glp_simplex(lp, &opt);

		int glp_status = wrapper::glp_get_status(lp);
		int glp_dual_st = wrapper::glp_get_dual_stat(lp);

		// Note:
		// - what we care about is the status if the primal problem
		// - if the presolver is used, glp_simplex might return GLP_ENOPFS/GLP_ENODFS (no feas primal/dual)
		//   while all statuses are GLP_UNDEF
		// - if we know that the dual problem is infeasible, then the primal has to be infeasible OR unbounded
		//   although we might not know which one
		//
		status =
			glp_status == GLP_OPT								? Status::OPTIMAL :
			glp_status == GLP_NOFEAS || glp_res == GLP_ENOPFS	? Status::INFEASIBLE :
			glp_status == GLP_UNBND								? Status::UNBOUNDED :
			glp_dual_st == GLP_NOFEAS || glp_res == GLP_ENODFS	? Status::INFEASIBLE_OR_UNBOUNDED :
			Status::ERROR;

	} else {
		glp_iptcp opt;
		wrapper::glp_init_iptcp(&opt);
		opt.msg_lev = msg_lev;	// debug info sent to terminal, default off

		wrapper::glp_interior(lp, &opt);

		// NOTE: glpk's interior point returns GLP_NOFEAS also for unbounded problems,
		//       not sure how we can check for boundedness
		//
		int glp_status = wrapper::glp_ipt_status(lp);
		// std::cout << "interior status: " << (status == GLP_OPT ? "GLP_OPT" : status == GLP_NOFEAS  ? "GLP_NOFEAS" : status == GLP_INFEAS ? "GLP_INFEAS" : status == GLP_UNDEF ? "GLP_UNDEF" : "XXX") << "\n";
		status =
			glp_status == GLP_OPT	? Status::OPTIMAL :
			glp_status == GLP_NOFEAS? Status::INFEASIBLE_OR_UNBOUNDED :
									  Status::ERROR;
	}

	// get optimal solution
	if(status == Status::OPTIMAL) {
		sol.set_size(n_var);
		for(uint j = 0; j < n_var; j++)
			sol.at(j) = is_interior ? wrapper::glp_ipt_col_prim(lp, j+1) : wrapper::glp_get_col_prim(lp, j+1);
	}

	// clean
	wrapper::glp_delete_prob(lp);
	wrapper::glp_free_env();

	return status == Status::OPTIMAL;

#endif // QIF_USE_GLPK
}

template<typename eT>
bool LinearProgram<eT>::ortools() {
	if(is_rat)
		std::cout << "\nWARNING: using OR-tools with rat. This will convert to double so it's not exact.\n\n";

#ifndef QIF_USE_ORTOOLS
	throw std::runtime_error("ortools not available");
#else
	using namespace operations_research;

	// ortools uses flags USE_<SOLVER> to enable the various solvers, but they don't seem to be "saved" in the headers
	// installed in the systme, instead the user is expected to pass them. Here we assume that USE_GLOP/USE_CLP are always
	// set (infact they are set in the qif header if USE_ORTOOLS is set). For the others (GUROBI/CPLEX), the user should
	// probably set -DUSE_<SOLVER> when compiling (not tested).
	//
	MPSolver::OptimizationProblemType ptype;
	if(solver ==  Solver::AUTO || solver ==  Solver::GLOP) {
		ptype = MPSolver::GLOP_LINEAR_PROGRAMMING;

	} else if(solver == Solver::CLP) {
		ptype = MPSolver::CLP_LINEAR_PROGRAMMING;

	} else if(solver == Solver::GUROBI) {
		#ifdef USE_GUROBI
		ptype = MPSolver::GUROBI_LINEAR_PROGRAMMING;
		break;
		#else
		throw std::runtime_error("GUROBI not enabled, please compile ortools with GUROBI and use -DUSE_GUROBI");
		#endif

	} else if(solver == Solver::CPLEX) {
		#ifdef USE_CPLEX
		ptype = MPSolver::CP:EX_LINEAR_PROGRAMMING;
		break;
		#else
		throw std::runtime_error("CPLEX not enabled, please compile ortools with GUROBI and use -DUSE_CPLEX");
		#endif

	} else {
		// GLPK/INTERNAL are handled by other methods
		throw std::runtime_error("shouldn't arrive here");
	}
	MPSolver orsolver("libqif", ptype);

	eT inf = infinity<eT>();
	const double sol_inf = orsolver.infinity();
	auto val = [&](eT a) -> double { return a == -inf ? -sol_inf : a == inf ? sol_inf : to_double(a); };

	MPObjective* const objective = orsolver.MutableObjective();
	if(maximize)
		objective->SetMaximization();
	else
		objective->SetMinimization();

	std::vector<MPVariable*> vars(n_var);
	std::vector<MPConstraint*> cons(n_con);

	for(uint x = 0; x < n_var; x++) {
		vars[x] = orsolver.MakeNumVar(val(var_lb[x]), val(var_ub[x]), "x"+std::to_string(x));

		objective->SetCoefficient(vars[x], val(obj_coeff[x]));
	}

	for(uint c = 0; c < n_con; c++)
		cons[c] = orsolver.MakeRowConstraint(val(con_lb[c]), val(con_ub[c]));

	for(auto& [key, value] : con_coeff) {
		auto& [col, row] = key;
		cons[row]->SetCoefficient(vars[col], val(value));
	}

	// set params
	MPSolverParameters param;
	param.SetIntegerParam(MPSolverParameters::LP_ALGORITHM,
		method == Method::SIMPLEX_PRIMAL ? MPSolverParameters::PRIMAL :
		method == Method::INTERIOR ? MPSolverParameters::BARRIER :
		MPSolverParameters::DUAL			// AUTO | SIMPLEX_DUAL
	);
	param.SetIntegerParam(MPSolverParameters::PRESOLVE, presolve ? MPSolverParameters::PRESOLVE_ON : MPSolverParameters::PRESOLVE_OFF);

	if(msg_level == MsgLevel::OFF)
		orsolver.SuppressOutput();
	else
		orsolver.EnableOutput();

	// go
	auto result_status = orsolver.Solve(param);

	status =
		result_status == MPSolver::OPTIMAL    ? Status::OPTIMAL :
		result_status == MPSolver::INFEASIBLE ? Status::INFEASIBLE :
		result_status == MPSolver::UNBOUNDED  ? Status::UNBOUNDED :
		Status::ERROR;

	if(status == Status::OPTIMAL) {
		sol.set_size(n_var);
		for(uint x = 0; x < n_var; x++)
			sol(x) = vars[x]->solution_value();
	}

	return status == Status::OPTIMAL;

#endif // QIF_USE_ORTOOLS
}


// transform the progarm in canonical form:
//        min  dot(c,x)
// subject to  A x == b
//               x >= 0
// b must be >= 0.
//
template<typename eT>
void LinearProgram<eT>::to_canonical_form() {

	uint n_var_orig = n_var;
	eT inf = infinity<eT>();

	if(var_transform.size() != 0)
		throw std::runtime_error("var_transform already set");

	// in canonical form, all variable bounds should be [0, infty]
	// we need to do various transformations, in the following we denote by x* the value of x in the original program
	//
	for(uint x = 0; x < n_var_orig; x++) {
		eT lb = var_lb[x];
		eT ub = var_ub[x];

		var_lb[x] = 0;
		var_ub[x] = inf;

		if(lb == -inf && ub == inf) {
			// An unbounded variable becomes two variables x* = x - xnew
			auto xnew = make_var(eT(0), inf);

			// recover x* as x - xnew
			var_transform.push_back(std::tuple(xnew, eT(1), eT(0)));

			// c*x* becomes c(x-xnew), so we need to update the obj
			set_obj_coeff(xnew, -obj_coeff[x]);
			
			// and for every coeff c of x in constraints, we need to add -c to xnew
			for(auto& [key, val] : con_coeff) {
				auto& [col, row] = key;
				if(col == x)
					set_con_coeff(row, xnew, -val);
			}

		} else if(lb == -inf) {
			// upper bounded variable, we set x = ub - x* (x* = ub - x)
			var_transform.push_back(std::tuple(-1, eT(-1), ub));

			obj_coeff[x] *= -1;

			for(auto& [key, val] : con_coeff) {
				auto& [col, row] = key;
				// a <= c x* <= b becomes a <= c(ub-x) <= b, so we need to subtract c*ub from lower/upper bound and then change sign of c
				if(col == x) {
					if(con_lb[row] != -inf)
						con_lb[row] -= val * ub;
					if(con_ub[row] != inf)
						con_ub[row] -= val * ub;

					val *= -1;
				}
			}

		} else { // lb != inf
			// lower or doubly bounded variable, we set x = x* - lb  (x* = x + lb)
			var_transform.push_back(std::tuple(-1, eT(1), lb));

			for(auto& [key, val] : con_coeff) {
				auto& [col, row] = key;
				// a <= c x* <= b becomes a <= c(x+lb) <= b, so we need to subtract c*lb from lower/upper bound
				if(col == x) {
					if(con_lb[row] != -inf)
						con_lb[row] -= val * lb;
					if(con_ub[row] != inf)
						con_ub[row] -= val * lb;
				}
			}

			// if an upper bound x* <= ub exists, it bocomes x+lb <= ub so we add a new constraint x <= ub - lb
			if(ub != inf) {
				uint con = make_con(-inf, ub - lb);
				set_con_coeff(con, x, eT(1));
			}
		}
	}

	// for every constraint lb <= cx <= ub with lb != ub, change ub to infty and add a separate consraint cx <= ub
	for(uint c = 0; c < n_con; c++) {
		eT lb = con_lb[c];
		eT ub = con_ub[c];

		if(lb != -inf && ub != inf && lb != ub) {
			con_ub[c] = inf;

			uint newc = make_con(-inf, ub);

			for(auto& [key, val] : con_coeff) {
				auto& [col, row] = key;
				if(row == c)
					set_con_coeff(newc, col, val);
			}
		}
	}

	// for every non-equality constraint, add slack variable
	for(uint c = 0; c < n_con; c++) {
		eT lb = con_lb[c];
		eT ub = con_ub[c];

		if(lb == -inf) {
			// upper bound, add slack newx to make equal
			con_lb[c] = ub;

			uint xnew = make_var(eT(0), inf);
			set_con_coeff(c, xnew, eT(1));

		} else if(ub == inf) {
			// lower bound, subtract slack newx to make equal
			con_ub[c] = lb;

			uint xnew = make_var(eT(0), inf);
			set_con_coeff(c, xnew, eT(-1));
		}
	}

	// invert constraints with negative constants
	for(uint c = 0; c < n_con; c++) {

		if(less_than(con_ub[c], eT(0))) {
			con_lb[c] *= -1;
			con_ub[c] *= -1;

			for(auto& [key, val] : con_coeff) {
				auto row = std::get<1>(key);
				if(row == c)
					val *= -1;
			}
		}
	}

	// canonical form is minimizing
	if(maximize) {
		maximize = false;

		for(auto& c : obj_coeff)
			c *= -1;
	}
}

// simplex
// Solve the linear program in canonical form
//        min  dot(c,x)
// subject to  A x == b
//               x >= 0
// b must be >= 0.
// 
// This is mainly to be used with rats
//
// The algorithm is the two-phase primal revised simplex method.
// In the first phase auxiliaries are created which we eliminate
// until we have a basis consisting solely of actual variables.
// This is pretty much the "textbook algorithm", and shouldn't
// be used for anything that matters. It doesn't exploit sparsity
// at all. You could use it with floating points but it wouldn't 
// work for anything except the most simple problem due to accumulated
// errors and the comparisons with zero.
//
// Summary of the algorithm:
// https://ocw.mit.edu/courses/sloan-school-of-management/15-093j-optimization-methods-fall-2009/lecture-notes/MIT15_093J_F09_lec04.pdf
//
template<typename eT>
bool LinearProgram<eT>::simplex() {
	using arma::zeros;
	using arma::ones;
	using arma::eye;
	using arma::uvec;

	// write program in matrix form
	const Col<eT> b = con_lb;
	const Row<eT> c = obj_coeff;
	Mat<eT> A = zeros<Mat<eT>>(n_con, n_var);	// use a dense matrix. The current algorithm doesn't use sparsity anyway, and operations on SpMat are much slower

	for(auto& [key, val] : con_coeff) {
		auto& [col, row] = key;
		A(row, col) = val;
	}

	assert(!maximize);
	for(uint i = 0; i < n_con; i++)
		assert(!less_than(b.at(i), eT(0)));
	
	uint n_con = this->n_con;		// to update locally

	Row<char> is_basic	= zeros<Row<char>>(n_var + n_con);
	uvec basic			= zeros<uvec>(n_con);				// indices of current basis
	Mat<eT> Binv		= eye<Mat<eT>>(n_con, n_con);		// inverse of basis matrix
	Row<eT> cB			= ones<Row<eT>>(n_con);				// costs of basic variables
	sol					= zeros<Col<eT>>(n_var + n_con);	// current solution

	Col<eT> BinvAs;

	// Intialize phase one of RSM by setting basis = auxiliaries.
	basic = arma::linspace<uvec>(n_var, n_var + n_con - 1, n_con); // n_var..n_var+n_con-1
	is_basic.tail(n_con).fill(1);
	sol.tail(n_con) = b;
	bool phase_one = true;

	auto pivot = [&](uint leaving, uint entering) -> void {
		// Update solution.
		eT ratio = sol(basic(leaving)) / BinvAs(leaving);
		sol(basic) -= ratio * BinvAs;
		sol(entering) = ratio;

		// Update basis inverse. Pivot on the `leaving' row, last col of [ Binv | BinvAs ]
		// Add to each row a multiple of Binv.row(leaving), so that the last col
		// becomes the identity vector with 1 on 'leaving'-row.
		//
		Binv.row(leaving) /= BinvAs(leaving);		// now the leaving row is normalized, there's 1 in the last col
		for(uint r = 0; r < n_con; r++)
			if(r != leaving)		// all rows except leaving row
				Binv.row(r) -= BinvAs(r) * Binv.row(leaving);

		// ... and variable status flags
		is_basic(basic(leaving)) = 0;
		is_basic(entering) = 1;
		cB(leaving) = phase_one ? eT(0) : c(entering);
		basic(leaving) = entering;
	};

	// Begin simplex iterations
	while(true) {
		// Calculate dual solution...
		Row<eT> pi = cB * Binv;

		// Use it to calculate the reduced costs of the variables. Don't
		// calculate for auxiliaries - they can't re-enter the basis.
		Row<eT> rc = (phase_one ? zeros<Row<eT>>(n_var) : c) - pi * A;
		uint entering;
		for(entering = 0; entering < n_var; entering++)					// crazy: using arma::index_min causes linker to complain about gmp symbols, even if rat is not used
			if(!is_basic(entering) && less_than(rc(entering), eT(0)))	// NOTE: Bland's rule, use the first index, to guarantee no cycles
				break;

		// If we couldn't find a variable with a negative reduced cost, 
		// we terminate this phase because we are at optimality for this
		// phase - not necessarily optimal for the actual problem.
		// std::cout << entering << "\n";
		if(entering == n_var) {
			if(phase_one) {
				phase_one = false;
				// Check objective - if 0, we are OK
				if(less_than(eT(0), accu(sol.tail(n_con)))) {
					// If any auxiliary still nonzero, we couldn't find a feasible
					// basis without auxiliaries.
					status = Status::INFEASIBLE;
					break;
				}

				// IMPORTANT: this fixes a bug in RationalSimplex.jl
				// if an artificial variable is still left in the basis (with 0 value) we need to drive it out
				// before starting Phase 2. Otherwise the artificial variable might become > 0 in the future!
				uvec keep = ones<uvec>(n_con);

				for(uint i = 0; i < n_con; i++) {
					if(basic(i) < n_var) continue; // non-artificial

					// Rule: check for non-zero elements in the i-th row of Binv*A
					Row<eT> t = Binv.row(i) * A;
					for(entering = 0; entering < n_var; entering++)
						if(!is_basic(entering) && !equal(t(entering), eT(0)))
							break;

					if(entering == n_var) {
						// all elements are zero so constraint i is redundant, mark to delete
						keep(i) = 0;
					} else {
						// found non-zero element (entering), pivot on that
						BinvAs = Binv * A.col(entering);
						pivot(i, entering);
					}
				}

				// remove marked constraints
				auto keep_i = arma::find(keep);
				if(size(keep_i).n_rows != n_con) {
					A = A.rows(keep_i);
					Binv = Binv.submat(keep_i, keep_i);
					basic = basic(keep_i);
					n_con = A.n_rows;
				}

				cB = c(basic).t();
				continue; // start phase 2

			} else {
				status = Status::OPTIMAL;
				break;
			}
		}

		// Calculate how the solution will change when our new
		// variable enters the basis and increases from 0
		BinvAs = Binv * A.col(entering);

		// Perform a "ratio test" on each variable to determine
		// which will reach 0 first
		uint leaving = n_con;
		eT min_ratio(0);
		for(uint j = 0; j < n_con; j++) {
			if(less_than(eT(0), BinvAs(j))) {
				eT ratio = sol(basic(j)) / BinvAs(j);
				if(less_than(ratio, min_ratio) || leaving == n_con) {
					min_ratio = ratio;
					leaving = j;
				}
			}
		}

		// If no variable will leave basis, then we have an 
		// unbounded problem.
		if(leaving == n_con) {
			status = Status::UNBOUNDED;
			break;
		}

		// ready to pivot
		pivot(leaving, entering);
	}

	sol = sol.head(n_var);		// the solution are the first n vars

	return status == Status::OPTIMAL;
}


template<typename eT>
void LinearProgram<eT>::dump() {
	std::cerr
		<< "var_lb: " << Row<eT>(var_lb)
		<< "var_ub: " << Row<eT>(var_ub)
		<< "con_lb: " << Row<eT>(con_lb)
		<< "con_ub: " << Row<eT>(con_ub)
		<< "obj_coeff: " << Row<eT>(obj_coeff)
		<< "con_coeff:\n";

	for(auto me : con_coeff)
		std::cout << "\t" << me.row << ", " << me.row << ", " << me.val << "\n";

	std::cout  << "\n";
}


template<typename eT>
string LinearProgram<eT>::to_mps() {
	throw std::runtime_error("not implemented");

	// using std::to_string;

	// string s;
	// s = "NAME PROG\n";

	// // rows
	// s += "ROWS\n";
	// s += " N  OBJ\n";	// objective function
	// for(uint i = 0; i < n_con; i++) {
	// 	char sense_i = sense.n_rows > i ? sense.at(i) : '<';	// default sense is <
	// 	string s_sense = sense_i == '<' ? "L" :
	// 					 sense_i == '>' ? "G" :
	// 						 			  "E";
	// 	s += " " + s_sense + " ROW" + to_string(i+1) + "\n";
	// }

	// // columns
	// s += "COLUMNS\n";
	// for(uint j = 0; j < n_var; j++) {
	// 	s += " X" + to_string(j+1) + " OBJ " + to_string(obj_coeff[j]) + "\n";

	// 	for(uint i = 0; i < n_var; i++)
	// 		s += " X" + to_string(j+1) + " ROW" + to_string(i+1) + " " + to_string(A.at(i, j)) + "\n";
	// }

	// // RHS
	// s += "RHS\n";
	// for(uint i = 0; i < n_con; i++) {
	// 	s += " RHS ROW" + to_string(i+1) + " " + to_string(con_ub[i]) + "\n";
	// }

	// // BOUNDS
	// eT inf = infinity<eT>();
	// s += "\nBOUNDS\n";
	// for(uint j = 0; j < n_var; j++) {
	// 	eT lb = var_lb[j];
	// 	eT ub = var_ub[j];

	// 	if(lb == -inf && ub == inf)
	// 		s += " FR BND X" + to_string(j+1) + "\n";
	// 	else(lb == ub)
	// 		s += " FX BND X" + to_string(j+1) + " " + to_string(lb) + "\n";
	// 	else {
	// 		if(ub != inf)
	// 			s += " UP BND X" + to_string(j+1) + " " + to_string(ub) + "\n";
	// 		else(lb != -inf)
	// 			s += " LO BND X" + to_string(j+1) + " " + to_string(lb) + "\n";
	// 	}
	// }

	// s += "ENDATA\n";

	// return s;
}

}

