#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_g_vuln_module(py::module m) {

	m.def("min_loss_given_max_vuln",	m::g_vuln::min_loss_given_max_vuln<double>, "pi"_a, "n_cols"_a, "n_guesses"_a, "max_vuln"_a, "gain"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("min_loss_given_max_vuln",	m::g_vuln::min_loss_given_max_vuln<rat>,    "pi"_a, "n_cols"_a, "n_guesses"_a, "max_vuln"_a, "gain"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

	m.def("min_vuln_given_max_loss",	m::g_vuln::min_vuln_given_max_loss<double>, "pi"_a, "n_cols"_a, "n_guesses"_a, "max_loss"_a, "gain"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("min_vuln_given_max_loss",	m::g_vuln::min_vuln_given_max_loss<rat>   , "pi"_a, "n_cols"_a, "n_guesses"_a, "max_loss"_a, "gain"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

}