namespace qp {

using std::string;

// Use a map <col, row> => value for sparse arrays. Put column first so that
// the lexicographic order used by map is top-to-bottom/left-to-right
template<typename eT> using Sparse = std::map<std::pair<uint,uint>,eT>;

enum class Status { OPTIMAL, INFEASIBLE, ERROR, STOPPED };
enum class Method { ADDM };

std::ostream& operator<<(std::ostream& os, const Status& status);
std::ostream& operator<<(std::ostream& os, const Method& method);


class Defaults {
	public:
		static bool		osqp_polish;
		static bool		osqp_verbose;
		static double	osqp_alpha;
		static double	osqp_eps_abs;
		static double	osqp_eps_rel;
		static Method	method;
};

// Solve the quadratic program
// min 1/2 x^T P x + c^T x
// subject to l <= A x <= u
// x >= 0

template<typename eT>
class QuadraticProgram {
	typedef uint Var;
	typedef uint Con;

	public:
		bool non_negative = false;
		Method method = Defaults::method;
		Status status;

		bool osqp_polish = Defaults::osqp_polish;
		bool osqp_verbose = Defaults::osqp_verbose;
		double osqp_alpha = Defaults::osqp_alpha;
		double osqp_eps_abs = Defaults::osqp_eps_abs;
		double osqp_eps_rel = Defaults::osqp_eps_rel;

		// Optional early stopping: called every early_stop_interval iterations with the current (approximate)
		// solution. If it returns true the solver stops with status STOPPED, and solution() returns that solution.
		std::function<bool(const Col<eT>&)> early_stop;
		uint early_stop_interval = 10;

		QuadraticProgram() {}

		bool solve();

		inline eT objective()				{ return obj; }
		inline eT solution(Var x)			{ return sol(x); }
		inline Col<eT> solution()			{ return sol; };

		void clear();
		void from_matrix(const arma::SpMat<eT>&P, const Col<eT>& c, const arma::SpMat<eT>& A, const Col<eT>& l, const Col<eT>& u);

		// API for creating variables and constraints
		Var make_var();
		std::vector<Var> make_vars(uint n);
		std::vector<std::vector<Var>> make_vars(uint n, uint m);

		Con make_con(eT lb, eT ub);

		void set_obj_coeff(Var var, eT coeff, bool add = false);				// linear part
		void set_obj_coeff(Var var1, Var var2, eT coeff, bool add = false);		// quadratic part
		void set_con_coeff(Con cons, Var var, eT coeff, bool add = false);

	protected:
		Col<eT> sol;			// solution
		eT obj;					// objective

		uint n_var = 0,							// number of variables
			 n_con = 0;							// number of constraints

		std::vector<OSQPFloat> obj_coeff_lin;		// coefficients for the objective function, linear part
		Sparse<eT> obj_coeff_quad;				// coefficients for the objective function, quadratic part
		Sparse<eT> con_coeff;					// coefficients for the constraints
		std::vector<OSQPFloat> con_lb, con_ub;	// constraints lower/upper

		bool osqp();
};

template<typename eT>
inline
typename QuadraticProgram<eT>::Var QuadraticProgram<eT>::make_var() {
	obj_coeff_lin.push_back(0);
	return n_var++;
}

template<typename eT>
inline
std::vector<typename QuadraticProgram<eT>::Var> QuadraticProgram<eT>::make_vars(uint n) {
	std::vector<Var> res(n);
	for(uint i = 0; i < n; i++)
		res[i] = make_var();
	return res;
}

template<typename eT>
inline
std::vector<std::vector<typename QuadraticProgram<eT>::Var>> QuadraticProgram<eT>::make_vars(uint n, uint m) {
	std::vector<std::vector<Var>> res(n);
	for(uint i = 0; i < n; i++) {
		res[i].resize(m);
		for(uint j = 0; j < m; j++)
			res[i][j] = make_var();
	}
	return res;
}

template<typename eT>
inline
typename QuadraticProgram<eT>::Con QuadraticProgram<eT>::make_con(eT lb, eT ub) {
	if(ub == infinity<eT>() && lb == -ub)
		throw std::runtime_error("trying to add unconstrained constraint");

	con_lb.push_back(lb);
	con_ub.push_back(ub);
	return n_con++;
}

template<typename eT>
inline
void QuadraticProgram<eT>::set_obj_coeff(Var var, eT coeff, bool add) {
	if(add)
		obj_coeff_lin[var] += coeff;
	else
		obj_coeff_lin[var] = coeff;
}

template<typename eT>
inline
void QuadraticProgram<eT>::set_obj_coeff(Var var1, Var var2, eT coeff, bool add) {
	if(equal<eT>(coeff, eT(0)))
		return;

	// maintain upper-triangular
	if(var1 > var2)
		std::swap(var1, var2);

	auto key = std::pair(var2, var1);		// <col, row>
	if(add && obj_coeff_quad.count(key))
		obj_coeff_quad[key] += coeff;
	else
		obj_coeff_quad[key] = coeff;
}

template<typename eT>
inline
void QuadraticProgram<eT>::set_con_coeff(Con con, Var var, eT coeff, bool add) {
	if(equal<eT>(coeff, eT(0)))
		return;

	auto key = std::pair(var, con);	// <col, row>
	if(add && con_coeff.count(key))
		con_coeff[key] += coeff;
	else
		con_coeff[key] = coeff;
}

template<typename eT>
void QuadraticProgram<eT>::clear() {
	obj_coeff_lin.clear();
	obj_coeff_quad.clear();
	con_coeff.clear();
	con_lb.clear();
	con_ub.clear();
	n_var = n_con = 0;
}

// input program in matrix form
// P: objective coeff, quadratic part
// c: objective constants, linear part
// A: constraint coefficients
// l: lower bounds for constraints
// u: upper bounds for constraints
//
template<typename eT>
void QuadraticProgram<eT>::from_matrix(const arma::SpMat<eT>&P, const Col<eT>& c, const arma::SpMat<eT>& A, const Col<eT>& l, const Col<eT>& u) {

	clear();
	n_var = A.n_cols;
	n_con = A.n_rows;

	if(P.n_rows != n_var || P.n_cols != n_var || c.n_elem != n_var || l.n_elem != n_con || u.n_elem != n_con)
		throw std::runtime_error("invalid size");
	
	auto Pend = P.end();
	for(auto c = P.begin(); c != Pend; ++c) {		// c++ throws weird warning, ++c doesn't!
		set_obj_coeff(c.row(), c.col(), *c);
	}

	using cfv = std::vector<OSQPFloat>;
	obj_coeff_lin = arma::conv_to<cfv>::from(c);
	
	auto Aend = A.end();
	for(auto c = A.begin(); c != Aend; ++c) {		// c++ throws weird warning, ++c doesn't!
		set_con_coeff(c.row(), c.col(), *c);
	}

	con_lb = arma::conv_to<cfv>::from(l);
	con_ub = arma::conv_to<cfv>::from(u);
}

template<typename eT>
bool QuadraticProgram<eT>::solve() {
	return osqp();
}

// Compressed Sparse Column (CSC) format, the arrays are owned by the object (OSQPCscMatrix only points to them).
// https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
//
struct CscMatrix {
	std::vector<OSQPFloat> val;		// array of non-zero values, in top-to-bottom, left-to-right order
	std::vector<OSQPInt> row_ind;	// for every value, the row index of that value
	std::vector<OSQPInt> col_ptr;	// for every column, the val index of the first element of that column.
									// Plus an extra element at the end, with the size of val
	OSQPCscMatrix mat;

	template<typename eT>
	CscMatrix(uint n_rows, uint n_cols, const Sparse<eT>& entries) {
		uint n_nonzero = entries.size();
		val.reserve(n_nonzero);
		row_ind.reserve(n_nonzero);
		col_ptr.reserve(n_cols + 1);

		for(auto& [key, value] : entries) {
			auto& [col, row] = key;

			// first time we see a column, update col_ptr
			while(col_ptr.size() <= col)	// possibly update previous empty columns
				col_ptr.push_back(val.size());

			val.push_back(OSQPFloat(value));
			row_ind.push_back(row);
		}
		while(col_ptr.size() <= n_cols)
			col_ptr.push_back(n_nonzero);

		mat = { OSQPInt(n_rows), OSQPInt(n_cols), col_ptr.data(), row_ind.data(), val.data(), OSQPInt(n_nonzero), -1, 0 };
	}

	CscMatrix(const CscMatrix&) = delete;	// mat points to our own vectors
};

template<typename eT>
bool QuadraticProgram<eT>::osqp() {

	// TODO
	if(non_negative)
		throw std::runtime_error("not_implemented");

	CscMatrix P(n_var, n_var, obj_coeff_quad);		// cost function, quadratic part
	CscMatrix A(n_con, n_var, con_coeff);			// constraints

	// settings
	OSQPSettings settings;
	wrapper::osqp_set_default_settings(&settings);
	settings.polishing = osqp_polish;
	settings.verbose = osqp_verbose;
	if(osqp_alpha   >= 0.0) settings.alpha = osqp_alpha;
	if(osqp_eps_abs >= 0.0) settings.eps_abs = osqp_eps_abs;
	if(osqp_eps_rel >= 0.0) settings.eps_rel = osqp_eps_rel;

	// With early_stop, we solve in chunks of early_stop_interval iterations. Each osqp_solve call continues from
	// the previous iterates (warm start), and the current solution is passed to early_stop after each chunk.
	const OSQPInt max_iter = settings.max_iter;
	if(early_stop) {
		settings.max_iter = early_stop_interval;
		settings.warm_starting = 1;
	}

	// solve
	OSQPSolver* solver = nullptr;
	OSQPInt err = wrapper::osqp_setup(&solver, &P.mat, obj_coeff_lin.data(), &A.mat, con_lb.data(), con_ub.data(), n_con, n_var, &settings);
	bool stopped = false;

	auto get_solution = [&]() {
		sol.set_size(n_var);
		for(uint j = 0; j < n_var; j++)
			sol.at(j) = solver->solution->x[j];
	};

	for(OSQPInt iter = 0; !err; ) {
		// at the end of each chunk OSQP checks for approximate termination, giving OSQP_SOLVED_INACCURATE.
		// Without early_stop this only happens after max_iter, so we treat it as OSQP_MAX_ITER_REACHED.
		err = wrapper::osqp_solve(solver);
		OSQPInt st = err ? 0 : solver->info->status_val;
		if(err || !early_stop || (st != OSQP_MAX_ITER_REACHED && st != OSQP_SOLVED_INACCURATE))
			break;

		get_solution();
		if(early_stop(sol)) {
			stopped = true;
			break;
		}

		iter += solver->info->iter;
		if(iter >= max_iter)
			break;
	}

	OSQPInt st = err ? 0 : solver->info->status_val;
	status =
		stopped														? Status::STOPPED :
		st == OSQP_SOLVED											? Status::OPTIMAL :
		st == OSQP_PRIMAL_INFEASIBLE || st == OSQP_DUAL_INFEASIBLE	? Status::INFEASIBLE :
		Status::ERROR;

	// get optimal solution
	if(status == Status::OPTIMAL) {
		get_solution();
		obj = solver->info->obj_val;
	}

	// cleanup
	if(solver)
		wrapper::osqp_cleanup(solver);

	return status == Status::OPTIMAL;
}


} // namespace qp

