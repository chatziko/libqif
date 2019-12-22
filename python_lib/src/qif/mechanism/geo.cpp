#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_geo_module(py::module m) {

	m.doc() = R"pbdoc(
		Mechanism construction for location privacy.
	)pbdoc";

	m.def("planar_laplace_draw",	m::geo::planar_laplace_draw<double>, "epsilon"_a);

	m.def("planar_laplace_grid",	m::geo::planar_laplace_grid<double>, "width"_a, "height"_a, "step"_a, "epsilon"_a);

	m.def("planar_laplace_grid",	overload<double,double>     (m::geo::planar_geometric_draw<double>), "cell_size"_a, "epsilon"_a);
	m.def("planar_laplace_grid",	overload<double,double,uint>(m::geo::planar_geometric_draw<double>), "cell_size"_a, "epsilon"_a, "n_samples"_a);

	m.def("planar_geometric_grid",	m::geo::planar_geometric_grid<double>, "width"_a, "height"_a, "step"_a, "epsilon"_a);

}