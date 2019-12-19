#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_l_risk_module(py::module m) {

	m.def("min_loss_given_min_risk",	m::l_risk::min_loss_given_min_risk<double>, "pi"_a, "n_cols"_a, "n_guesses"_a, "min_risk"_a, "adv_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("min_loss_given_min_risk",	m::l_risk::min_loss_given_min_risk<rat>,    "pi"_a, "n_cols"_a, "n_guesses"_a, "min_risk"_a, "adv_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

	m.def("max_risk_given_max_loss",	m::l_risk::max_risk_given_max_loss<double>, "pi"_a, "n_cols"_a, "n_guesses"_a, "max_loss"_a, "adv_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<double>());
	m.def("max_risk_given_max_loss",	m::l_risk::max_risk_given_max_loss<rat>   , "pi"_a, "n_cols"_a, "n_guesses"_a, "max_loss"_a, "adv_loss"_a, "loss"_a, "hard_max_loss"_a = infinity<rat>   ());

}