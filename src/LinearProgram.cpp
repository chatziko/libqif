#include "qif"

namespace qif {
namespace lp {


bool        Defaults::glp_presolve  = true;
msg_level_t Defaults::glp_msg_level = msg_level_t::off;
method_t    Defaults::method        = method_t::simplex_dual;


std::ostream& operator<<(std::ostream& os, const status_t& status) {
	std::string s[] = { "optimal", "infeasible", "unbounded", "infeasible_or_unbounded", "error" };
	return os << s[static_cast<uint>(status)];
}

std::ostream& operator<<(std::ostream& os, const method_t& method) {
	std::string s[] = { "simplex_primal", "simplex_dual", "simplex_dualp", "interior" };
	return os << s[static_cast<uint>(method)];
}

std::ostream& operator<<(std::ostream& os, const msg_level_t& level) {
	std::string s[] = { "off", "err", "on", "all" };
	return os << s[static_cast<uint>(level)];
}

}}
