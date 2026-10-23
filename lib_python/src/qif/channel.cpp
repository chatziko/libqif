#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_channel_compose_module(py::module);


void init_channel_module(py::module m) {

	init_channel_compose_module(m.def_submodule("compose", ""));

	m.def("normalize",      overload<const  chan&>(channel::normalize<double>), "C"_a);
	m.def("normalize",      overload<const rchan&>(channel::normalize<rat>   ), "C"_a);

	m.def("identity",		[](uint n, double_c_t) { return channel::identity<double>(n); }, "n_rows"_a, "type"_a = def_type());
	m.def("identity",		[](uint n, rat_c_t   ) { return channel::identity<rat>   (n); }, "n_rows"_a, "type"_a = def_type());

	m.def("no_interference",[](uint n, uint m, double_c_t) { return channel::no_interference<double>(n, m); }, "n_rows"_a, "n_cols"_a = 1, "type"_a = def_type());
	m.def("no_interference",[](uint n, uint m, rat_c_t   ) { return channel::no_interference<rat>   (n, m); }, "n_rows"_a, "n_cols"_a = 1, "type"_a = def_type());

	m.def("randu",			[](uint n, uint m, double_c_t) { return channel::randu<double>(n, m); }, "n_rows"_a, "n_cols"_a = 0, "type"_a = def_type());
	m.def("randu",			[](uint n, uint m, rat_c_t   ) { return channel::randu<rat>   (n, m); }, "n_rows"_a, "n_cols"_a = 0, "type"_a = def_type());

	m.def("deterministic",	[](std::function<uint(uint)> f, uint n, uint m, double_c_t) { return channel::deterministic<double>(f, n, m); }, "map"_a, "n_rows"_a, "n_cols"_a, "type"_a = def_type());
	m.def("deterministic",	[](std::function<uint(uint)> f, uint n, uint m, rat_c_t   ) { return channel::deterministic<rat   >(f, n, m); }, "map"_a, "n_rows"_a, "n_cols"_a, "type"_a = def_type());

	m.def("is_proper",      channel::is_proper<double>, "C"_a, "mrd"_a = def_mrd<double>);
	m.def("is_proper",      channel::is_proper<rat>,    "C"_a, "mrd"_a = rat(0));

	m.def("assert_proper",  channel::assert_proper<double>, "C"_a);
	m.def("assert_proper",  channel::assert_proper<rat>,    "C"_a);

	m.def("posterior",      channel::posterior<double>, "C"_a, "pi"_a, "col"_a);
	m.def("posterior",      channel::posterior<rat>,    "C"_a, "pi"_a, "col"_a);

	m.def("posteriors",     channel::posteriors<double>, "C"_a, "pi"_a = prob());
	m.def("posteriors",     channel::posteriors<rat>,    "C"_a, "pi"_a);
	m.def("posteriors",     [](const rchan& C) { return channel::posteriors<rat>(C); }, "C"_a);	// pi = rprob() default causes a weird "vector out of range" error on windows.

	m.def("hyper",     		channel::hyper<double>, "C"_a, "pi"_a);
	m.def("hyper",     		channel::hyper<rat>,    "C"_a, "pi"_a);

	m.def("reduced",   		channel::reduced<double>, "C"_a);
	m.def("reduced",   		channel::reduced<rat>,    "C"_a);

	m.def("iterative_bayesian_update", channel::iterative_bayesian_update<double>, "C"_a, "out"_a, "start"_a = prob(),        "max_diff"_a = 1e-6, "max_iter"_a = 0);
	// For rat, use None defaults: start = rprob() causes a weird "vector out of range" error on windows, and a double max_diff
	// makes this overload fail to match, so calls with rat arguments would silently fall back to the double one.
	m.def("iterative_bayesian_update",
		[](const rchan& C, const rprob& out, py::object start, py::object max_diff, uint max_iter) {
			return channel::iterative_bayesian_update<rat>(C, out,
				start.is_none()    ? rprob()          : start.cast<rprob>(),
				max_diff.is_none() ? rat(1, 1000000)  : max_diff.cast<rat>(),
				max_iter);
		},
		"C"_a, "out"_a, "start"_a = py::none(), "max_diff"_a = py::none(), "max_iter"_a = 0
	);

	m.def("factorize",
		[](const chan& A, const chan& B, bool col_stoch, const std::string& method)   { return channel::factorize(A, B, col_stoch, method); },
		"A"_a, "B"_a, "col_stoch"_a = false, "method"_a = "auto");
	m.def("factorize",
		[](const rchan& A, const rchan& B, bool col_stoch, const std::string& method) { return channel::factorize(A, B, col_stoch, method); },
		"A"_a, "B"_a, "col_stoch"_a = false, "method"_a = "auto");

	m.def("left_factorize", channel::left_factorize<double>, "A"_a, "B"_a, "col_stoch"_a = false, "method"_a = "auto");
	m.def("left_factorize", channel::left_factorize<rat>,    "A"_a, "B"_a, "col_stoch"_a = false, "method"_a = "auto");

	m.def("sum_column_min", channel::sum_column_min<double>, "C"_a);
	m.def("sum_column_min", channel::sum_column_min<rat>,    "C"_a);

	m.def("sample",     		overload<const  chan&, const  prob&>(channel::sample<double>), "C"_a, "pi"_a);
	m.def("sample",     		overload<const rchan&, const rprob&>(channel::sample<rat>),    "C"_a, "pi"_a);

	m.def("sample",     		overload<const  chan&, const  prob&, uint>(channel::sample<double>), "C"_a, "pi"_a, "n_samples"_a);
	m.def("sample",     		overload<const rchan&, const rprob&, uint>(channel::sample<rat>),    "C"_a, "pi"_a, "n_samples"_a);

}