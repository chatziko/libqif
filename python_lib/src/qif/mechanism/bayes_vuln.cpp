#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_bayes_vuln_module(py::module m) {

	m.def("min_loss_given_max_vuln",	m::bayes_vuln::min_loss_given_max_vuln<double>, "pi"_a, "n_cols"_a, "max_vuln"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("min_loss_given_max_vuln",	m::bayes_vuln::min_loss_given_max_vuln<rat>,    "pi"_a, "n_cols"_a, "max_vuln"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

	m.def("min_vuln_given_max_loss",	m::bayes_vuln::min_vuln_given_max_loss<double>, "pi"_a, "n_cols"_a, "max_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("min_vuln_given_max_loss",	m::bayes_vuln::min_vuln_given_max_loss<rat>,    "pi"_a, "n_cols"_a, "max_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

	m.def("min_vuln_for_row",			m::bayes_vuln::min_vuln_for_row<double>, "pi"_a, "p"_a, "C"_a);
	m.def("min_vuln_for_row",			m::bayes_vuln::min_vuln_for_row<rat>,    "pi"_a, "p"_a, "C"_a);

}