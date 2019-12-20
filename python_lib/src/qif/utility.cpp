#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::utility;


void init_utility_module(py::module m) {

	m.def("expected_distance",	overload<const  chan&,              const  prob&,const  chan&>(expected_distance<double>), "D"_a, "pi"_a, "C"_a);
	m.def("expected_distance",	overload<const Metric<double,uint>&,const  prob&,const  chan&>(expected_distance<double>), "d"_a, "pi"_a, "C"_a);
	m.def("expected_distance",	overload<const rchan&,              const rprob&,const rchan&>(expected_distance<rat>   ), "D"_a, "pi"_a, "C"_a);
	m.def("expected_distance",	overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(expected_distance<rat>   ), "d"_a, "pi"_a, "C"_a);

}