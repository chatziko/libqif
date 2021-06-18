#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_geo_ind_module(py::module m) {

	m.doc() = R"pbdoc(
		Mechanism construction for geo-indistinguishability.
	)pbdoc";

	m.def("planar_laplace_sample",	m::geo_ind::planar_laplace_sample<double>, "epsilon"_a);

	m.def("planar_laplace_grid",	m::geo_ind::planar_laplace_grid<double>, "width"_a, "height"_a, "step"_a, "epsilon"_a);

	m.def("planar_geometric_sample",	overload<double,double>     (m::geo_ind::planar_geometric_sample<double>), "cell_size"_a, "epsilon"_a);
	m.def("planar_geometric_sample",	overload<double,double,uint>(m::geo_ind::planar_geometric_sample<double>), "cell_size"_a, "epsilon"_a, "n_samples"_a);

	m.def("planar_geometric_grid",	m::geo_ind::planar_geometric_grid<double>, "width"_a, "height"_a, "step"_a, "epsilon"_a);

}