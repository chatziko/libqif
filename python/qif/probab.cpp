#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_probab_module(py::module m) {

	m.def("uniform",		[](uint n, dtype_d) { return probab::uniform<double>(n); }, "n_elem"_a, "dtype"_a = float64);
	m.def("uniform",		[](uint n, dtype_r) { return probab::uniform<rat>   (n); }, "n_elem"_a, "dtype"_a);

	m.def("dirac",			[](uint n, uint x, dtype_d) { return probab::dirac<double>(n, x); }, "n_elem"_a, "x"_a = 0, "dtype"_a = float64);
	m.def("dirac",			[](uint n, uint x, dtype_r) { return probab::dirac<rat>   (n, x); }, "n_elem"_a, "x"_a = 0, "dtype"_a);

	m.def("randu",			[](uint n, dtype_d) { return probab::randu<double>(n); }, "n_elem"_a, "dtype"_a = float64);
	m.def("randu",			[](uint n, dtype_r) { return probab::randu<rat>   (n); }, "n_elem"_a, "dtype"_a);

	m.def("normalize",      overload<const  prob&>(probab::normalize<double>), "pi"_a);
	m.def("normalize",      overload<const rprob&>(probab::normalize<rat>   ), "pi"_a);

	m.def("draw",     		overload<const  prob&>(probab::draw<double>), "pi"_a);
	m.def("draw",     		overload<const rprob&>(probab::draw<rat>),    "pi"_a);

	m.def("draw",     		overload<const  prob&, uint>(probab::draw<double>), "pi"_a, "n_samples"_a);
	m.def("draw",     		overload<const rprob&, uint>(probab::draw<rat>),    "pi"_a, "n_samples"_a);

	m.def("is_uniform",     probab::is_uniform<double>, "pi"_a, "mrd"_a = def_mrd<double>);
	m.def("is_uniform",     probab::is_uniform<rat>,    "pi"_a, "mrd"_a = rat(0));

	m.def("is_proper",      probab::is_proper<double>, "pi"_a, "mrd"_a = def_mrd<double>);
	m.def("is_proper",      probab::is_proper<rat>,    "pi"_a, "mrd"_a = rat(0));

	m.def("assert_proper",  probab::assert_proper<double>, "pi"_a);
	m.def("assert_proper",  probab::assert_proper<rat>,    "pi"_a);

	m.def("project_to_simplex", overload<const  prob&>(probab::project_to_simplex<double>), "pi"_a);
	m.def("project_to_simplex", overload<const rprob&>(probab::project_to_simplex<rat>),    "pi"_a);

	m.def("to_grid",  		probab::to_grid<double>, "pi"_a, "width"_a);
	m.def("to_grid",  		probab::to_grid<rat>,    "pi"_a, "width"_a);

	m.def("from_grid",  	probab::from_grid<double>, "grid"_a);
	m.def("from_grid",  	probab::from_grid<rat>,    "grid"_a);

}
