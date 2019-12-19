#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_shannon_module(py::module m) {

	m.def("max_entropy_given_same_loss",	m::shannon::max_entropy_given_same_loss<double>, "pi"_a, "out"_a, "loss"_a, "md"_a = def_md<double>, "mrd"_a = def_mrd<double>);

}