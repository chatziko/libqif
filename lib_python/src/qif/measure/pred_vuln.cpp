#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_pred_vuln_module(py::module m) {

	m.doc() = R"pbdoc(
		Vulnerability of guessing a predicate of the secret.
	)pbdoc";

	m.def("G_pred",			[](const arma::urowvec& P, double_c_t) { return pred_vuln::G_pred<double>(P); }, "P"_a, "type"_a = def_type());
	m.def("G_pred",			[](const arma::urowvec& P, rat_c_t   ) { return pred_vuln::G_pred<rat>   (P); }, "P"_a, "type"_a = def_type());

	m.def("prior",      	pred_vuln::prior<double>, "P"_a, "pi"_a);
	m.def("prior",      	pred_vuln::prior<rat>,    "P"_a, "pi"_a);

	m.def("posterior",     	pred_vuln::posterior<double>, "P"_a, "pi"_a, "C"_a);
	m.def("posterior",     	pred_vuln::posterior<rat>,    "P"_a, "pi"_a, "C"_a);

	m.def("add_leakage",   	pred_vuln::add_leakage<double>, "P"_a, "pi"_a, "C"_a);
	m.def("add_leakage",   	pred_vuln::add_leakage<rat>,    "P"_a, "pi"_a, "C"_a);

	m.def("mult_leakage",  	pred_vuln::mult_leakage<double>, "P"_a, "pi"_a, "C"_a);
	m.def("mult_leakage",  	pred_vuln::mult_leakage<rat>,    "P"_a, "pi"_a, "C"_a);

	m.def("mult_capacity",  pred_vuln::mult_capacity<double>, "P"_a, "C"_a, "method"_a = "direct");
	m.def("mult_capacity",  pred_vuln::mult_capacity<rat>,    "P"_a, "C"_a, "method"_a = "direct");

	m.def("binary_channel", pred_vuln::binary_channel<double>, "P"_a, "pi"_a, "C"_a);
	m.def("binary_channel",	pred_vuln::binary_channel<rat>,    "P"_a, "pi"_a, "C"_a);

}