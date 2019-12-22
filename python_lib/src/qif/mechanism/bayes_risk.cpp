#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_bayes_risk_module(py::module m) {

	m.doc() = R"pbdoc(
		Mechanism construction for Bayes risk.
	)pbdoc";

	m.def("min_loss_given_min_risk",	m::bayes_risk::min_loss_given_min_risk<double>, "pi"_a, "n_cols"_a, "min_risk"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("min_loss_given_min_risk",	m::bayes_risk::min_loss_given_min_risk<rat>,    "pi"_a, "n_cols"_a, "min_risk"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

	m.def("max_risk_given_max_loss",	m::bayes_risk::max_risk_given_max_loss<double>, "pi"_a, "n_cols"_a, "max_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("max_risk_given_max_loss",	m::bayes_risk::max_risk_given_max_loss<rat>,    "pi"_a, "n_cols"_a, "max_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

	m.def("max_risk_for_row",			m::bayes_risk::max_risk_for_row<double>, "pi"_a, "p"_a, "C"_a);
	m.def("max_risk_for_row",			m::bayes_risk::max_risk_for_row<rat>,    "pi"_a, "p"_a, "C"_a);

}