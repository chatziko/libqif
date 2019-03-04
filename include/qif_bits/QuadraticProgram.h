namespace qp {

using std::string;

enum class Status { OPTIMAL, INFEASIBLE, ERROR };
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
	public:
		arma::SpMat<eT>
			P,				// cost function, quadratic part
			A;				// constraints
		Col<eT>
			x,				// solution
			l,				// lower-bound constants
			u,				// upper-bound constants
			c;				// cost function, linear part

		bool non_negative = false;
		Method method = Defaults::method;
		Status status;

		bool osqp_polish = Defaults::osqp_polish;
		bool osqp_verbose = Defaults::osqp_verbose;
		double osqp_alpha = Defaults::osqp_alpha;
		double osqp_eps_abs = Defaults::osqp_eps_abs;
		double osqp_eps_rel = Defaults::osqp_eps_rel;

		QuadraticProgram() {}

		bool solve();

		inline eT objective()				{ return eT(1)/2 * arma::dot(x.t() * P, x) + arma::dot(x, c); }

	protected:
		inline void check_sizes() {
			if(A.n_rows != l.n_rows ||
			   A.n_rows != u.n_rows ||
			   A.n_cols != P.n_rows ||
			   A.n_cols != P.n_cols ||
			   A.n_cols != c.n_rows
			) throw std::runtime_error("invalid size");
		}

		bool osqp();
};


template<typename eT>
bool QuadraticProgram<eT>::solve() {
	check_sizes();

	return osqp();
}

template<typename eT>
csc* to_csc(const arma::SpMat<eT>& M) {
	// Compressed Sparse Column (CSC) format.
	// https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
	// val:      array of non-zero values, in top-to-bottom, left-to-right order
	// row_ind:  for every value, the row index of that value
	// col_ptr:  for every column, the val index of the first element of that column.
	//           Plus an extra element at the end, with the size of val

	c_float* val = (c_float*) malloc(sizeof(c_float) * M.n_nonzero);
	c_int* row_ind = (c_int*) malloc(sizeof(c_int) * M.n_nonzero);
	c_int* col_ptr = (c_int*) malloc(sizeof(c_int) * M.n_cols + 1);

	uint i = 0;
	for(uint col = 0; col < M.n_cols; col++) {
		col_ptr[col] = i;

		for(auto it = M.begin_col(col); it != M.end_col(col); ++it) {
			row_ind[i] = it.row();
			val[i] = *it;
			i++;
		}
	}
	col_ptr[M.n_cols] = M.n_nonzero;

	return wrapper::csc_matrix(M.n_rows, M.n_cols, M.n_nonzero, val, row_ind, col_ptr);
}

inline void free_csc(csc* matrix) {
	c_free(matrix->p);
	c_free(matrix->i);
	c_free(matrix->x);
	c_free(matrix);
}

template<typename eT>
bool QuadraticProgram<eT>::osqp() {

	// TODO
	if(non_negative)
		throw std::runtime_error("not_implemented");

	// osqp's element type is c_float (double by default).
	// this might be different from our eT, so we might need to convert
	//
	Col<c_float> temp_c = arma::conv_to<Col<c_float>>::from(c);
	Col<c_float> temp_l = arma::conv_to<Col<c_float>>::from(l);
	Col<c_float> temp_u = arma::conv_to<Col<c_float>>::from(u);

	// populate data
	OSQPData* data = (OSQPData *)malloc(sizeof(OSQPData));
	data->n = A.n_cols;				// number of variables
	data->m = A.n_rows;				// number of constraints
	data->P = to_csc(P);			// cost function, quadratic part
	data->q = temp_c.memptr();		// cost function, linear part (we call it c, osqp calls it q)
	data->A = to_csc(A);			// constraints
	data->l = temp_l.memptr();		// lower bounds
	data->u = temp_u.memptr();		// upper bounds

	// settings
	OSQPSettings* settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
	wrapper::osqp_set_default_settings(settings);
	settings->polish = osqp_polish;
	settings->verbose = osqp_verbose;
	if(osqp_alpha   >= 0.0) settings->alpha = osqp_alpha;
	if(osqp_eps_abs >= 0.0) settings->eps_abs = osqp_eps_abs;
	if(osqp_eps_rel >= 0.0) settings->eps_rel = osqp_eps_rel;

	// solve
	OSQPWorkspace* work = wrapper::osqp_setup(data, settings);
	wrapper::osqp_solve(work);

	c_int st = work->info->status_val;
	status =
		st == OSQP_SOLVED											? Status::OPTIMAL :
		st == OSQP_PRIMAL_INFEASIBLE || st == OSQP_DUAL_INFEASIBLE	? Status::INFEASIBLE :
		Status::ERROR;

	// get optimal solution
	if(status == Status::OPTIMAL) {
		x.set_size(A.n_cols);
		for(uint j = 0; j < A.n_cols; j++)
			x.at(j) = work->solution->x[j];
	}

	// cleanup
	free_csc(data->P);
	free_csc(data->A);
	free(data);
	free(settings);
	wrapper::osqp_cleanup(work);

	return status == Status::OPTIMAL;
}


} // namespace qp

