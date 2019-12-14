#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_channel_module(py::module);
void init_probab_module(py::module);
void init_metric_module(py::module);


PYBIND11_MODULE(qif, m) {
	arma::arma_rng::set_seed_random();

	// Init mp++'s pybind11 integration
	mppp_pybind11::init();

	init_channel_module(m.def_submodule("channel", ""));
	init_probab_module (m.def_submodule("probab", ""));
	init_metric_module (m.def_submodule("metric", ""));
}
