#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_bayes_risk_module(py::module m) {

	m.def("prior",      			bayes_risk::prior<double>, "pi"_a);
	m.def("prior",      			bayes_risk::prior<rat>,    "pi"_a);

	m.def("posterior",     			bayes_risk::posterior<double>, "pi"_a, "C"_a);
	m.def("posterior",     			bayes_risk::posterior<rat>,    "pi"_a, "C"_a);

	m.def("add_leakage",   			bayes_risk::add_leakage<double>, "pi"_a, "C"_a);
	m.def("add_leakage",   			bayes_risk::add_leakage<rat>,    "pi"_a, "C"_a);

	m.def("mult_leakage",  			bayes_risk::mult_leakage<double>, "pi"_a, "C"_a);
	m.def("mult_leakage",  			bayes_risk::mult_leakage<rat>,    "pi"_a, "C"_a);

	m.def("mult_capacity",  		bayes_risk::mult_capacity<double>, "C"_a, "method"_a = "direct");
	m.def("mult_capacity",  		bayes_risk::mult_capacity<rat>,    "C"_a, "method"_a = "direct");

	m.def("strategy",  				bayes_risk::strategy<double>, "pi"_a, "C"_a);
	m.def("strategy",  				bayes_risk::strategy<rat>,    "pi"_a, "C"_a);

}