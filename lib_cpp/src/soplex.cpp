// SoPlex wrapper. Separate from wrapper.cpp (and excluded from the unity build) so that only this file sees
// the soplex (and boost) headers.

#pragma GCC diagnostic push		// ignore warnings in soplex code
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#include <soplex.h>
#pragma GCC diagnostic pop
#include "qif"

namespace qif::wrapper {

using namespace lp;

// RAII array of mpq_t, so that it can be passed to soplex's GMP interface
class MpqArray {
	std::vector<__mpq_struct> v;
	public:
		MpqArray(size_t n) : v(n)	{ for(auto& q : v) mpq_init(&q); }
		~MpqArray()					{ for(auto& q : v) mpq_clear(&q); }
		MpqArray(const MpqArray&) = delete;

		mpq_t* data()				{ return reinterpret_cast<mpq_t*>(v.data()); }
		mpq_t& operator[](size_t i)	{ return data()[i]; }
		rat get(size_t i)			{ return rat(data()[i]); }

		void set(size_t i, const rat& x) {
			mpz_set(mpq_numref(data()[i]), x.get_num().get_mpz_view());
			mpz_set(mpq_denref(data()[i]), x.get_den().get_mpz_view());
		}
};

std::string soplex_solve(bool maximize, const std::string& method, bool presolve, const std::string& msg_level,
	const std::vector<rat>& obj_coeff, const std::map<std::pair<uint,uint>,rat>& con_coeff,
	const std::vector<rat>& var_lb, const std::vector<rat>& var_ub,
	const std::vector<rat>& con_lb, const std::vector<rat>& con_ub,
	std::vector<rat>& sol) {

	using soplex::SoPlex;
	SoPlex sp;

	// exact solving: rational mode (iterative refinement), zero tolerances
	sp.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
	sp.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
	sp.setIntParam(SoPlex::READMODE, SoPlex::READMODE_RATIONAL);
	sp.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);
	sp.setRealParam(SoPlex::FEASTOL, 0.0);
	sp.setRealParam(SoPlex::OPTTOL, 0.0);

	sp.setIntParam(SoPlex::OBJSENSE, maximize ? SoPlex::OBJSENSE_MAXIMIZE : SoPlex::OBJSENSE_MINIMIZE);
	sp.setIntParam(SoPlex::ALGORITHM, method == Method::SIMPLEX_PRIMAL ? SoPlex::ALGORITHM_PRIMAL : SoPlex::ALGORITHM_DUAL);
	sp.setIntParam(SoPlex::SIMPLIFIER, presolve ? SoPlex::SIMPLIFIER_INTERNAL : SoPlex::SIMPLIFIER_OFF);
	sp.setIntParam(SoPlex::VERBOSITY,
		msg_level == MsgLevel::ALL ? SoPlex::VERBOSITY_FULL :
		msg_level == MsgLevel::ON  ? SoPlex::VERBOSITY_NORMAL :
		                             SoPlex::VERBOSITY_ERROR);

	// soplex treats bounds with absolute value >= its INFTY param as infinite
	const rat inf = infinity<rat>();
	const rat sp_inf(sp.realParam(SoPlex::INFTY)), sp_minf = -sp_inf;
	auto bound = [&](const rat& x) -> const rat& { return x == inf ? sp_inf : x == -inf ? sp_minf : x; };

	const uint n_var = obj_coeff.size(),
			   n_con = con_lb.size();

	// add rows (with no coefficients)
	MpqArray lhs(1), rhs(1);
	for(uint i = 0; i < n_con; i++) {
		lhs.set(0, bound(con_lb[i]));
		rhs.set(0, bound(con_ub[i]));
		sp.addRowRational(&lhs[0], nullptr, nullptr, 0, &rhs[0]);
	}

	// add columns. con_coeff is ordered by <col,row> so all entries of each column are consecutive
	MpqArray obj(1), lb(1), ub(1), vals(n_con);
	std::vector<int> rows(n_con);
	auto it = con_coeff.begin();

	for(uint j = 0; j < n_var; j++) {
		int size = 0;
		for(; it != con_coeff.end() && it->first.first == j; ++it) {
			rows[size] = it->first.second;
			vals.set(size, it->second);
			size++;
		}

		obj.set(0, obj_coeff[j]);
		lb.set(0, bound(var_lb[j]));
		ub.set(0, bound(var_ub[j]));
		sp.addColRational(&obj[0], &lb[0], vals.data(), rows.data(), size, &ub[0]);
	}

	// solve
	auto sp_status = sp.optimize();

	using SPStatus = soplex::SPxSolverBase<soplex::Real>;
	std::string status =
		sp_status == SPStatus::OPTIMAL    ? Status::OPTIMAL :
		sp_status == SPStatus::INFEASIBLE ? Status::INFEASIBLE :
		sp_status == SPStatus::UNBOUNDED  ? Status::UNBOUNDED :
		sp_status == SPStatus::INForUNBD  ? Status::INFEASIBLE_OR_UNBOUNDED :
		                                    Status::ERROR;

	if(status == Status::OPTIMAL) {
		MpqArray primal(n_var);
		if(!sp.getPrimalRational(primal.data(), n_var))
			return Status::ERROR;

		sol.resize(n_var);
		for(uint j = 0; j < n_var; j++)
			sol[j] = primal.get(j);
	}

	return status;
}

} // namespace qif::wrapper
