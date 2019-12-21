#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_bayes_vuln_module(py::module);
void init_bayes_risk_module(py::module);
void init_g_vuln_module(py::module);
void init_l_risk_module(py::module);
void init_shannon_module(py::module);
void init_guessing_module(py::module);
void init_d_privacy_module(py::module);


void init_measure_module(py::module m) {

	m.doc() = R"pbdoc(
		.. autosummary::
			:toctree: _autosummary
			:template: template.rst

			bayes_vuln
			bayes_risk
			g_vuln
			l_risk
			shannon
			guessing
			d_privacy
	)pbdoc";

	init_bayes_vuln_module(m.def_submodule("bayes_vuln", ""));
	init_bayes_risk_module(m.def_submodule("bayes_risk", ""));
	init_g_vuln_module    (m.def_submodule("g_vuln",     ""));
	init_l_risk_module    (m.def_submodule("l_risk",     ""));
	init_shannon_module   (m.def_submodule("shannon",    ""));
	init_guessing_module  (m.def_submodule("guessing",   ""));
	init_d_privacy_module (m.def_submodule("d_privacy",  ""));

}