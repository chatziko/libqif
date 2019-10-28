#include "qif"

namespace qif::qp {


Method    Defaults::method       = Method::ADDM;
bool        Defaults::osqp_polish  = true;
bool        Defaults::osqp_verbose = false;
double      Defaults::osqp_alpha   = -1.0;		// negative means
double      Defaults::osqp_eps_abs = -1.0;		// keep OSQP's
double      Defaults::osqp_eps_rel = -1.0;		// defaults


std::ostream& operator<<(std::ostream& os, const Status& status) {
	std::string s[] = { "optimal", "infeasible", "error" };
	return os << s[static_cast<uint>(status)];
}

std::ostream& operator<<(std::ostream& os, const Method& method) {
	std::string s[] = { "addm" };
	return os << s[static_cast<uint>(method)];
}

}
