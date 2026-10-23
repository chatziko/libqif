# Changelog

All notable changes to libqif are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- `qif.lp`: SoPlex solver, solves programs exactly (both `rat` and `double`).
- `channel.factorize`, `left_factorize`: `"project_qp"` (Euclidean projection by quadratic programming) and `"project_lp"` (Manhattan projection by linear programming, exact for `rat`) methods.
- `refinement.refined_by_witness`, returns a counter-example gain function when refinement fails (replaces `refined_by(..., method="project")`). Works also for `rat` (exactly).
- `method` argument in `refinement.max_refined_by`.

### Changed
- `qif.lp`: `rat` programs are solved with SoPlex by default (instead of the internal simplex), which is much faster on larger programs.
- `channel.factorize`, `left_factorize`: methods renamed to `"solve_lp"` and `"project_subgrad"` (the old names `"lp"`, `"subgrad"` still work).
- `channel.factorize`, `left_factorize`, `refinement.refined_by`: for `double`, `method="auto"` uses the subgradient method with a budget of iterations, falling back to `"project_qp"` (or `"solve_lp"` for wide matrices). Previously it used the subgradient method without a budget for large matrices (which could be very slow) and LP otherwise.
- `refinement.refined_by`: `method` is any `channel.factorize` method, and it always returns a `bool`. `method="factorize"` and `"project"` are deprecated.
- `refinement.refined_by_witness` (previously `refined_by(..., "project")`): the counter-example is checked directly, so a false answer is always correct, and the QP stops as soon as one is found (~2x faster when refinement fails). `R` is empty when refinement fails.
- C++: `channel::factorize_lp`, `factorize_subgrad` renamed to `factorize_solve_lp`, `factorize_project_subgrad`.

### Fixed
- `rat` matrices: `max()` of negative values and `min()` of positive values returned 0.
- `channel.factorize` with the subgradient method could run forever on non-factorizable channels (wrong stopping bound).

## [1.5.0] - 2026-10-22

### Added
- New documentation site at [libqif.chatzi.org](https://libqif.chatzi.org), with one page per module.
- `qif.lp`: HiGHS solver.
- `method` argument in `channel.factorize`, `channel.left_factorize` (`"lp"`, `"subgrad"`) and `metric.kantorovich` (`"lp"`, `"fastemd"`).
- `bayes_vuln.mult_capacity_bound_cap` in python.
- `metric.manhattan` in python, `rat` support in `metric.l1`, `total_variation`, `convex_separation(_quasi)`, `kantorovich`, `from_distance_matrix` and `to_distance_matrix`.

### Changed
- `metric.euclidean(uint)` and `metric.discrete(uint)` return distances of the default type (`rat` after `set_default_type(rat)`, previously always `double`).
- `qif.lp`: the default solver is CLP (instead of GLOP), and HiGHS for interior point.
- `qif.metric`: the combinators work for metrics on any domain.
- `channel.hyper` removes zero-probability columns, `bayes_risk.mult_leakage` returns `inf` for zero posterior risk.
- Python: more consistent argument names.
- Upgraded Armadillo to 15.6.2 and OSQP to 1.0.0. CMake 4 works without extra flags.
- Tests are in python (`tests/`, pytest), `tests_cpp` only covers C++ internals (LP/QP solvers, `equal`).

### Removed
- C++: `refinement::refined_by(A, B, G)`, use `refined_by(A, B, G, R)` instead.
- C++: `metric::threshold_bin`.

### Fixed
- `mechanism.l_risk` ignored `n_guesses`.
- Windows: functions using linear system solving (`channel.factorize` with the subgradient method, `left_factorize` for large channels, `mechanism.d_privacy.tight_constraints`) failed, the wheels were built without LAPACK.
- `rat` infinity (e.g. `bayes_risk.mult_capacity` of channels with zero risk) was `2^31-1` on Windows, `2^63-1` elsewhere.
- `channel.left_factorize` failed for matrices with 1000 or more elements (and `channel.factorize` with `col_stoch` using the subgradient method).
- `channel.iterative_bayesian_update` and `channel.posteriors` silently computed in `double` for `rat` arguments, unless all optional arguments were given.

## [1.4.0] - 2026-10-20

### Added
- Wheels for macOS arm64.

### Changed
- Support for OR-Tools 9, with OR-Tools 9.15 bundled in the wheels on all platforms.
- Python is only needed when building the Python library.
- Updated the installation instructions for OR-Tools 9, Homebrew and CMake 4.

## [1.3.0] - 2026-10-20

### Changed
- Wheels for Python 3.11 to 3.15. Python 3.11 or later is required.
- An x86 CPU with AVX2 (Haswell or later) is required by default.

### Fixed
- Building with Python 3.12 or later.

## [1.2.4] - 2025-09-04

### Changed
- Dropped Python 3.8, added 3.12, 3.13 and 3.14.
- Upgraded pybind11 and mp++.

### Fixed
- Removed the deprecated `std::iterator`, and fixed a warning in `QuadraticProgram`.
- The sample code in the README.

## [1.2.3] - 2023-06-06

### Changed
- Dropped Python 3.7, added 3.11.

## [1.2.2] - 2022-05-06

### Changed
- Dropped Python 3.6, added 3.10.
- Linux wheels are built on manylinux2014.
- Simpler installation: Armadillo is no longer a prerequisite, and programs always link with `-lgmp -lmp++`.

### Fixed
- `mechanism.d_privacy`: method names.

## [1.2.1] - 2021-02-11

### Changed
- Upgraded Armadillo to 10.2.0.

### Fixed
- The `mechanism` package was not imported by `import qif`.

## [1.2.0] - 2020-12-28

### Added
- Type stubs (`.pyi`) for the Python library.

### Changed
- The Python library is organised in packages.

### Fixed
- `channel.hyper`: division by zero.

## [1.1.10] - 2020-12-27

### Changed
- Renamed `dirac` to `point` and `draw` to `sample`.
- The Python library is a `qif` package that imports the `_qif` extension module.
- Repository layout: `lib_cpp`, `lib_python` and `tests_cpp`.

## [1.1.9] - 2020-12-20

### Changed
- Upgraded Armadillo to 10.1.2, pybind11 to 2.6.1 and mp++ to 0.21.

## [1.1.8] - 2020-12-19

### Fixed
- Windows wheels: fixes to bundling `msvcp140.dll`.

## [1.1.7] - 2020-12-19

### Changed
- Python 3.6 or later is required.

### Fixed
- Windows wheels bundle `msvcp140.dll`.

## [1.1.6] - 2020-12-18

### Fixed
- `measure.bayes_vuln.strategy`: the size of the output.

## [1.1.5] - 2020-12-18

### Added
- Windows support (MSVC) and Windows wheels.

### Changed
- Armadillo is included as a submodule, and older Armadillo versions are no longer supported.

## [1.1.4] - 2020-12-14

### Changed
- Wheels are built with cibuildwheel, and the documentation is published to GitHub Pages.

## [1.1.3] - 2020-12-12

### Changed
- Dropped Python 3.5, added 3.9.
- Upgraded OR-Tools to 8.1.

## [1.1.2] - 2020-01-13

### Fixed
- Python: non-contiguous vectors are copied before being passed to the library.

## [1.1.1] - 2019-12-24

### Added
- `mechanism.d_privacy.min_loss_given_d` in Python.

### Changed
- Renamed `mechanism.geo` to `mechanism.geo_ind`.

## [1.1.0] - 2019-12-23

### Added
- `measure.pred_vuln` and `measure.pred_risk`.
- `metric.optimize.simplex_l1_min_enclosing_ball` and `metric.optimize.simplex_project`.

### Fixed
- `measure.g_vuln.g_to_bayes`.

## [1.0.0] - 2019-12-22

### Added
- Python: `qif.metric.optimize`, `qif.lp`, `channel.compose`.
- `measure.l_risk.l_zero_one`.
- A clear error message when numpy is not available.

## [0.8.0] - 2019-12-21

First tagged release: C++ library and Python bindings for quantitative information flow (channels, leakage
measures, metrics, mechanisms, refinement), distributed as wheels on PyPI.

[Unreleased]: https://github.com/chatziko/libqif/compare/v1.5.0...HEAD
[1.5.0]: https://github.com/chatziko/libqif/compare/v1.4.0...v1.5.0
[1.4.0]: https://github.com/chatziko/libqif/compare/v1.3.0...v1.4.0
[1.3.0]: https://github.com/chatziko/libqif/compare/v1.2.4...v1.3.0
[1.2.4]: https://github.com/chatziko/libqif/compare/v1.2.3...v1.2.4
[1.2.3]: https://github.com/chatziko/libqif/compare/v1.2.2...v1.2.3
[1.2.2]: https://github.com/chatziko/libqif/compare/v1.2.1...v1.2.2
[1.2.1]: https://github.com/chatziko/libqif/compare/v1.2.0...v1.2.1
[1.2.0]: https://github.com/chatziko/libqif/compare/v1.1.10...v1.2.0
[1.1.10]: https://github.com/chatziko/libqif/compare/v1.1.9...v1.1.10
[1.1.9]: https://github.com/chatziko/libqif/compare/v1.1.8...v1.1.9
[1.1.8]: https://github.com/chatziko/libqif/compare/v1.1.7...v1.1.8
[1.1.7]: https://github.com/chatziko/libqif/compare/v1.1.6...v1.1.7
[1.1.6]: https://github.com/chatziko/libqif/compare/v1.1.5...v1.1.6
[1.1.5]: https://github.com/chatziko/libqif/compare/v1.1.4...v1.1.5
[1.1.4]: https://github.com/chatziko/libqif/compare/v1.1.3...v1.1.4
[1.1.3]: https://github.com/chatziko/libqif/compare/v1.1.2...v1.1.3
[1.1.2]: https://github.com/chatziko/libqif/compare/v1.1.1...v1.1.2
[1.1.1]: https://github.com/chatziko/libqif/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/chatziko/libqif/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/chatziko/libqif/compare/v0.8.0...v1.0.0
[0.8.0]: https://github.com/chatziko/libqif/releases/tag/v0.8.0
