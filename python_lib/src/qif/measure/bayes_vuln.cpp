#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_bayes_vuln_module(py::module m) {

	m.doc() = R"pbdoc(
		Bayes vulnerability.
	)pbdoc";

	m.def("prior",      			bayes_vuln::prior<double>, "pi"_a);
	m.def("prior",      			bayes_vuln::prior<rat>,    "pi"_a);

	m.def("posterior",     			bayes_vuln::posterior<double>, "pi"_a, "C"_a);
	m.def("posterior",     			bayes_vuln::posterior<rat>,    "pi"_a, "C"_a);

	m.def("add_leakage",   			bayes_vuln::add_leakage<double>, "pi"_a, "C"_a);
	m.def("add_leakage",   			bayes_vuln::add_leakage<rat>,    "pi"_a, "C"_a);

	m.def("mult_leakage",  			bayes_vuln::mult_leakage<double>, "pi"_a, "C"_a);
	m.def("mult_leakage",  			bayes_vuln::mult_leakage<rat>,    "pi"_a, "C"_a);

	m.def("min_entropy_leakage",	bayes_vuln::min_entropy_leakage<double>, "pi"_a, "C"_a);

	m.def("mult_capacity",  		bayes_vuln::mult_capacity<double>, "C"_a);
	m.def("mult_capacity",  		bayes_vuln::mult_capacity<rat>,    "C"_a);

	m.def("strategy",  				bayes_vuln::strategy<double>, "pi"_a, "C"_a);
	m.def("strategy",  				bayes_vuln::strategy<rat>,    "pi"_a, "C"_a);

}