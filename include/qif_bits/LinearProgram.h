/*
The canonical_form and simplex methods are adapted from
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

// for some reason this needs to be included here, not enough to do it in "qif"
#include <cassert>

namespace lp {

using std::string;

enum class status_t { optimal, infeasible, unbounded, infeasible_or_unbounded, error };
enum class method_t { simplex_primal, simplex_dual, simplex_dualp, interior };
enum class msg_level_t { off, err, on, all };

std::ostream& operator<<(std::ostream& os, const status_t& status);
std::ostream& operator<<(std::ostream& os, const method_t& method);
std::ostream& operator<<(std::ostream& os, const msg_level_t& level);


class Defaults {
	public:
		static bool        glp_presolve ;
		static msg_level_t glp_msg_level;
		static method_t    method;
};

// Solve the linear program
// {min/max} dot(c,x)
// subject to lb <= A x <= ub
// lb <= x <= ub

template<typename eT>
class LinearProgram {
	typedef uint Var;
	typedef uint Con;

	// workaround a bug in g++ 4.9 (jessie), std::is_same<A,B> fails when A,B
	// are aliases,  so we use __gmp_expr<mpq_t, mpq_t> instead of rat
	const bool is_rat = std::is_same<eT, __gmp_expr<mpq_t, mpq_t>>::value;

	public:
		bool maximize = true;
		method_t method = is_rat ? method_t::simplex_primal : Defaults::method;	// rat only supports simplex_primal
		bool glp_presolve = is_rat ? false : Defaults::glp_presolve;
		status_t status;
		msg_level_t glp_msg_level = Defaults::glp_msg_level;

		bool solve();
		string to_mps();

		eT optimum();
		inline eT solution(Var x)	{ return sol(x); }
		inline Col<eT> solution()	{ return sol; };

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
		bool simplex();

		void to_canonical_form();
		Col<eT> original_solution();

		uint n_var = 0,							// number of variables
			 n_con = 0;							// number of constraints

		std::vector<eT> obj_coeff;				// coefficients for the objective function
		std::list<MatrixEntry<eT>> con_coeff;	// coefficients for the constraints
		std::vector<eT> var_lb, var_ub,			// variables lower/upper bounds
						con_lb, con_ub;			// constraints lower/upper

		// info for transforming solution to the original one (see canonical_form, original_solution)
		std::vector< std::tuple<int,eT,eT> > var_transform;
};


template<typename eT>
inline
eT LinearProgram<eT>::optimum() {
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
		int v2;
		eT coeff, add;
		std::tie(v2, coeff, add) = var_transform[v];

		res(v) = sol(v) * coeff + add - (v2 >= 0 ? sol(v2) : 0);
	}

	return res;
}

template<typename eT>
inline
typename LinearProgram<eT>::Var LinearProgram<eT>::make_var(eT lb, eT ub) {
	var_lb.push_back(lb);
	var_ub.push_back(ub);
	obj_coeff.push_back(0);
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
	if(add) {
		for(auto& me : con_coeff)
			if(me.row == con && me.col == var) {
				me.val += coeff;
				return;
			}
	}
	con_coeff.push_back(MatrixEntry<eT>(con, var, coeff));
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

	obj_coeff.clear();
	con_coeff.clear();
	var_lb.clear();
	var_ub.clear();
	con_lb.clear();
	con_ub.clear();

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
	return glpk();
}

// for rats, we use the simplex() method after transforming to canonical form
//
template<>
inline
bool LinearProgram<rat>::solve() {

	if(method != method_t::simplex_primal)
		throw std::runtime_error("not supported");

	LinearProgram<rat> lp(*this);	// clone
	lp.to_canonical_form();

	bool res = lp.simplex();
	status = lp.status;

	if(res)
		sol = lp.original_solution();

	return res;
}


template<typename eT>
bool LinearProgram<eT>::glpk() {
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

		wrapper::glp_set_col_bnds(lp, j+1, type, lb, ub);

		wrapper::glp_set_obj_coef(lp, j+1, obj_coeff[j]);				// coefficient in the cost functoin
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

		wrapper::glp_set_row_bnds(lp, i+1, type, lb, ub);
	}

	// loop over non-zero elements of sparse array
	int index = 1;
	for(auto me : con_coeff) {
		ia[index] = me.row + 1;
		ja[index] = me.col + 1;
		ar[index] = me.val;
		index++;
	}

	wrapper::glp_load_matrix(lp, size, &ia[0], &ja[0], &ar[0]);

	const int glp_msg_levs[] = { GLP_MSG_OFF, GLP_MSG_ERR, GLP_MSG_ON, GLP_MSG_ALL };
	const int msg_lev = glp_msg_levs[static_cast<uint>(glp_msg_level)];

	// solve
	const bool is_interior = method == method_t::interior;
	if(!is_interior) {	// simplex primal/dual/dualp
		glp_smcp opt;
		wrapper::glp_init_smcp(&opt);
		opt.meth =
			method == method_t::simplex_primal ? GLP_PRIMAL :
			method == method_t::simplex_dual   ? GLP_DUAL :
			GLP_DUALP;
		opt.msg_lev = msg_lev;							// debug info sent to terminal, default off
		opt.presolve = glp_presolve ? GLP_ON : GLP_OFF;	// use presolver

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
			glp_status == GLP_OPT								? status_t::optimal :
			glp_status == GLP_NOFEAS || glp_res == GLP_ENOPFS	? status_t::infeasible :
			glp_status == GLP_UNBND								? status_t::unbounded :
			glp_dual_st == GLP_NOFEAS || glp_res == GLP_ENODFS	? status_t::infeasible_or_unbounded :
		  status_t::error;

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
			glp_status == GLP_OPT	? status_t::optimal :
			glp_status == GLP_NOFEAS? status_t::infeasible_or_unbounded :
									  status_t::error;
	}

	// get optimal solution
	if(status == status_t::optimal) {
		sol.set_size(n_var);
		for(uint j = 0; j < n_var; j++)
			sol.at(j) = is_interior ? wrapper::glp_ipt_col_prim(lp, j+1) : wrapper::glp_get_col_prim(lp, j+1);
	}

	// clean
	wrapper::glp_delete_prob(lp);
	wrapper::glp_free_env();

	return status == status_t::optimal;
}

template<>
inline
bool LinearProgram<rat>::glpk() {
	throw std::runtime_error("not available for rat");
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
			uint xnew = make_var(0, inf);

			// recover x as x1 - x2
			var_transform.push_back(std::make_tuple(xnew, 1, 0));

			// "c * x*" becomes "c * (x - xnew)", so we need to update the obj
			set_obj_coeff(xnew, -obj_coeff[x]);
			
			// and for every coeff c of x in constraints, we need to add -c to x2
			for(auto me : con_coeff)
				if(me.col == x)
					set_con_coeff(me.row, xnew, -me.val);

		} else if(lb == inf) {
			// upper bounded variable, we set x* = ub - x
			var_transform.push_back(std::make_tuple(-1, -1, ub));

			obj_coeff[x] *= -1;

			for(auto& me : con_coeff) {
				if(me.col == x) {
					if(con_lb[me.row] != -inf)
						con_lb[me.row] += me.val * -ub;
					else if(con_ub[me.row] != inf)
						con_ub[me.row] += me.val * -ub;
				}
				me.val *= -1;
			}

		} else {
			// lower or doubly bounded variable, we set x* = x - l
			var_transform.push_back(std::make_tuple(-1, 1, -lb));

			for(auto me : con_coeff)
				if(me.col == x) {
					if(con_lb[me.row] != -inf)
						con_lb[me.row] += me.val * lb;
					else if(con_ub[me.row] != inf)
						con_ub[me.row] += me.val * lb;
				}

			// if an upper bound exists, we add a new constraint x <= lb + ub
			if(ub != inf) {
				uint con = make_con(-inf, lb+ub);
				set_con_coeff(con, x, 1);
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

			for(auto me : con_coeff)
				if(me.row == c)
					set_con_coeff(newc, me.col, me.val);
		}
	}

	// for every non-equality constraint, add slack variable
	for(uint c = 0; c < n_con; c++) {
		eT lb = con_lb[c];
		eT ub = con_ub[c];

		if(lb == -inf) {
			// upper bound, add slack newx to make equal
			con_lb[c] = ub;

			uint xnew = make_var(0, inf);
			set_con_coeff(c, xnew, 1);

		} else if(ub == inf) {
			// lower bound, subtract slack newx to make equal
			con_ub[c] = lb;

			uint xnew = make_var(0, inf);
			set_con_coeff(c, xnew, -1);
		}
	}

	// invert constraints with negative constants
	for(uint c = 0; c < n_con; c++) {

		if(con_ub[c] < eT(0)) {
			con_lb[c] *= -1;
			con_ub[c] *= -1;

			for(auto& me : con_coeff)
				if(me.row == c)
					me.val *= -1;
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
template<typename eT>
bool LinearProgram<eT>::simplex() {
	using arma::zeros;
	using arma::ones;
	using arma::eye;
	using arma::umat;

	// write program in matrix form
	Row<eT> b = con_lb,
			c = obj_coeff;
	Mat<eT> Adense = zeros<Mat<eT>>(n_con, n_var);	// use a dense matrix. The current algorithm doesn't use sparsity anyway, and operations on SpMat are much slower

	for(auto me : con_coeff)
		Adense(me.row, me.col) = me.val;

	uint m = n_con,
		 n = n_var;

	assert(!maximize);
	for(uint i = 0; i < m; i++)
		assert(!less_than(b.at(i), eT(0)));

	Mat<char> is_basic	= zeros<Mat<char>>(n + m);
	umat basic			= zeros<umat>(m);				// indices of current basis
	Mat<eT> Binv		= eye<Mat<eT>>(m, m);			// inverse of basis matrix
	Row<eT> cB			= ones<Row<eT>>(m);				// costs of basic variables
	sol					= zeros<Col<eT>>(n + m);		// current solution

	// Intialize phase 1
	for(uint i = 0; i < m; i++) {
		basic(i) = i + n;
		is_basic(i + n) = 1;
		sol(i + n) = b(i);
	}
	bool phase_one = true;

	// Begin simplex iterations
	while(true) {
		// Calculate dual solution...
		Row<eT> pi_T = cB * Binv;

		// ... and thus the reduced costs
		int entering = -1;
		for(uint j = 0; j < n; j++) {
			if(is_basic(j)) continue;
			eT rc = (phase_one ? eT(0) : c(j)) - dot(pi_T, Adense.col(j));
			if(less_than(rc, eT(0))) {
				entering = j;
				break;
			}
		}

		// If we couldn't find a variable with a negative reduced cost, 
		// we terminate this phase because we are at optimality for this
		// phase - not necessarily optimal for the actual problem.
		if(entering == -1) {
			if(phase_one) {
				phase_one = false;
				// Check objective - if 0, we are OK
				for(uint j = n; j < n + m; j++) {
					if(less_than(eT(0), sol(j))) {
						// It couldn't reduce objective to 0 which is equivalent
						// to saying a feasible basis with no artificials could
						// not be found
						status = status_t::infeasible;
						// std::cout << "-- infeasible\n";
						goto EXIT;
					}
				}
				// Start again in phase 2 with our nice feasible basis
				for(uint i = 0; i < m; i++) {
					cB(i) = basic(i) > n ? eT(0) : c(basic(i));
				}
				continue;
			} else {
				status = status_t::optimal;
				// std::cout << "-- optimal\n";
				goto EXIT;
			}
		}

		// Calculate how the solution will change when our new
		// variable enters the basis and increases from 0
		Col<eT> BinvAs = Binv * Adense.col(entering);

		// Perform a "ratio test" on each variable to determine
		// which will reach 0 first
		int leaving = -1;
		eT min_ratio = eT(0);
		for(uint i = 0; i < m; i++) {
			if(less_than(eT(0), BinvAs(i))) {
				eT ratio = sol(basic(i)) / BinvAs(i);
				if(less_than(ratio, min_ratio) || leaving == -1) {
					min_ratio = ratio;
					leaving = i;
				}
			}
		}

		// If no variable will leave basis, then we have an 
		// unbounded problem.
		if(leaving == -1) {
			status = status_t::unbounded;
			// std::cout << "--unbounded\n";
			goto EXIT;
		}

		// Now we update solution...
		for(uint i = 0; i < m; i++) {
			sol(basic(i)) -= min_ratio * BinvAs(i);
		}
		sol(entering) = min_ratio;

		// ... and the basis inverse...
		// Our tableau is ( Binv b | Binv | BinvAs )
		// and we doing a pivot on the leaving row of BinvAs
		eT pivot_value = BinvAs(leaving);
		for(uint i = 0; i < m; i++) {  // all rows except leaving row
			if(i == static_cast<uint>(leaving)) continue;
			eT factor = BinvAs(i) / pivot_value;
			for(uint j = 0; j < m; j++)
				Binv(i, j) -= factor * Binv(leaving, j);
		}
		for(uint j = 0; j < m; j++)
			Binv(leaving, j) /= pivot_value;

		// ... and variable status flags
		is_basic(basic(leaving)) = 0;
		is_basic(entering) = 0;
		cB(leaving) = phase_one ? eT(0) : c(entering);
		basic(leaving) = entering;
	}

EXIT:
	sol = sol.subvec(0, n-1);		// the solution are the first n vars

	return status == status_t::optimal;
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

template<>
inline
string LinearProgram<rat>::to_mps() {
	// TODO: make to_mps work for rat
	throw std::runtime_error("not supported");
}

}

