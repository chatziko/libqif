#include "qif"

namespace qif {
namespace qp {


bool        Defaults::osqp_polish  = true;
bool        Defaults::osqp_verbose = false;
double      Defaults::osqp_alpha   = 1.0;
method_t    Defaults::method       = method_t::addm;


std::ostream& operator<<(std::ostream& os, const status_t& status) {
	std::string s[] = { "optimal", "infeasible", "error" };
	return os << s[static_cast<uint>(status)];
}

std::ostream& operator<<(std::ostream& os, const method_t& method) {
	std::string s[] = { "addm" };
	return os << s[static_cast<uint>(method)];
}

}}
