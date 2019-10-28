#include "qif"

namespace qif::lp {


bool     Defaults::presolve  = true;
MsgLevel Defaults::msg_level = MsgLevel::OFF;
Method   Defaults::method    = Method::AUTO;
Solver   Defaults::solver    = Solver::AUTO;


std::ostream& operator<<(std::ostream& os, const Status& status) {
	std::string s[] = { "OPTIMAL", "INFEASIBLE", "UNBOUNDED", "INFEASIBLE_OR_UNBOUNDED", "ERROR" };
	return os << s[static_cast<uint>(status)];
}

std::ostream& operator<<(std::ostream& os, const Method& method) {
	std::string s[] = { "AUTO", "SIMPLEX_PRIMAL", "SIMPLEX_DUAL", "INTERIOR" };
	return os << s[static_cast<uint>(method)];
}

std::ostream& operator<<(std::ostream& os, const Solver& solver) {
	std::string s[] = { "AUTO", "INTERNAL", "GLPK", "GLOP", "CLP", "GUROBI", "CPLEX" };
	return os << s[static_cast<uint>(solver)];
}

std::ostream& operator<<(std::ostream& os, const MsgLevel& level) {
	std::string s[] = { "OFF", "ERR", "ON", "ALL" };
	return os << s[static_cast<uint>(level)];
}

} // namespace qif::lp
