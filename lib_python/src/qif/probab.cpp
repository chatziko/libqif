#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_probab_module(py::module m) {

	m.doc() = R"pbdoc(
		Probability distributions.
	)pbdoc";

	m.def("uniform",		[](uint n, double_c_t) { return probab::uniform<double>(n); }, "n_elem"_a, "type"_a = def_type());
	m.def("uniform",		[](uint n, rat_c_t   ) { return probab::uniform<rat>   (n); }, "n_elem"_a, "type"_a = def_type());

	m.def("point",			[](uint n, uint x, double_c_t) { return probab::point<double>(n, x); }, "n_elem"_a, "x"_a = 0, "type"_a = def_type());
	m.def("point",			[](uint n, uint x, rat_c_t   ) { return probab::point<rat>   (n, x); }, "n_elem"_a, "x"_a = 0, "type"_a = def_type());

	m.def("randu",			[](uint n, double_c_t) { return probab::randu<double>(n); }, "n_elem"_a, "type"_a = def_type());
	m.def("randu",			[](uint n, rat_c_t   ) { return probab::randu<rat>   (n); }, "n_elem"_a, "type"_a = def_type());

	m.def("normalize",      overload<const  prob&>(probab::normalize<double>), "pi"_a);
	m.def("normalize",      overload<const rprob&>(probab::normalize<rat>   ), "pi"_a);

	m.def("sample",     	overload<const  prob&>(probab::sample<double>), "pi"_a);
	m.def("sample",     	overload<const rprob&>(probab::sample<rat>),    "pi"_a);

	m.def("sample",    		overload<const  prob&, uint>(probab::sample<double>), "pi"_a, "n_samples"_a);
	m.def("sample",    		overload<const rprob&, uint>(probab::sample<rat>),    "pi"_a, "n_samples"_a);

	m.def("is_uniform",     probab::is_uniform<double>, "pi"_a, "mrd"_a = def_mrd<double>);
	m.def("is_uniform",     probab::is_uniform<rat>,    "pi"_a, "mrd"_a = rat(0));

	m.def("is_proper",      probab::is_proper<double>, "pi"_a, "mrd"_a = def_mrd<double>);
	m.def("is_proper",      probab::is_proper<rat>,    "pi"_a, "mrd"_a = rat(0));

	m.def("assert_proper",  probab::assert_proper<double>, "pi"_a);
	m.def("assert_proper",  probab::assert_proper<rat>,    "pi"_a);

	m.def("to_grid",  		probab::to_grid<double>, "pi"_a, "width"_a);
	m.def("to_grid",  		probab::to_grid<rat>,    "pi"_a, "width"_a);

	m.def("from_grid",  	probab::from_grid<double>, "grid"_a);
	m.def("from_grid",  	probab::from_grid<rat>,    "grid"_a);

}
