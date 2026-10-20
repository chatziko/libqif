# Installing libqif for C++

These instructions are for installing `libqif` to use with C++. For python, see the [README](README.md).

## Install via Homebrew

The easiest way to install libqif (on macOS or Linux) is via [Homebrew](https://brew.sh/), which also installs all
dependencies (including OR-Tools):
```bash
brew install --HEAD chatziko/tap/libqif
brew test --HEAD libqif
```

To upgrade to the latest version:
```bash
brew reinstall libqif
```

Homebrew installs everything under its prefix (`brew --prefix`), eg `/opt/homebrew` on Apple Silicon or
`/home/linuxbrew/.linuxbrew` on Linux, so you need to pass this prefix to the compiler (see below).

## Compile a sample program

Create a simple ```test.cpp``` file:
```c++
#include <qif>
using namespace qif;

int main() {
    chan C("1 0 0; 0 1 0; 0 0 1");
    prob pi = probab::uniform<double>(3);
    std::cout
        << "Bayes vulnerability of \n"
        << C << " under " << pi << " is "
        << measure::bayes_vuln::posterior(pi, C) << "\n";
}
```

Compile and run with:
```bash
g++ test.cpp -std=c++17 -lqif -larmadillo -lgmp -lmp++ -o test
./test

# or with clang
clang++ test.cpp -std=c++17 -lqif -larmadillo -lgmp -lmp++ -o test

# if installed via Homebrew
g++ test.cpp -std=c++17 -I$(brew --prefix)/include -L$(brew --prefix)/lib -lqif -larmadillo -lgmp -lmp++ -o test
```

You can find more sample programs in the [samples](https://github.com/chatziko/libqif/tree/master/samples) directory.

If libqif was built with OR-Tools you also need to link with `-lortools`.

## Build libqif from source

Prerequisites

* [CMake](https://cmake.org/) >= 3.24
* [GMP](https://gmplib.org/)
* [GSL](https://www.gnu.org/software/gsl/)
* BLAS/LAPACK (eg [OpenBLAS](https://www.openblas.net/))
* A C++17 compliant compiler (eg g++ or clang)

Optionally

* [OR-Tools](https://developers.google.com/optimization/) (recommended)
* [GLPK](https://www.gnu.org/software/glpk/)
* Python 3 with headers (only for building the python library)

On Ubuntu, these can be installed with:
```bash
sudo apt-get install g++ cmake libgmp-dev libgsl-dev libopenblas-dev python3-dev libglpk-dev
```

Get the code (note the `--recursive` to fetch the submodules).
```bash
git clone --recursive https://github.com/chatziko/libqif.git
```

To compile / install:
```bash
mkdir <path>/build && cd <path>/build
cmake ..
make
sudo make install
```

Notes:
* With CMake 4 you need to pass `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` (required by the included armadillo).
* On x86, libqif is compiled with `-march=haswell` by default (requires AVX2). Use eg `-DMARCH=native`
  to optimize for your CPU.

To run the tests
```bash
make tests_cpp
./tests_cpp/run
```

To build the samples:
```bash
make samples
./samples/<sample>
```

#### Use OR-Tools

Several libqif methods can use OR-Tools for linear optimization and network flow.
When installing via Homebrew, OR-Tools are installed automatically.

When installing manually, you should install OR-Tools (version 9) before compiling libqif.
The easiest way is to download the C++ `.tar.gz` (or `.zip`) file for your system from the
[OR-Tools releases](https://github.com/google/or-tools/releases), extract it, and pass its location to cmake:
```bash
sudo mkdir /opt/ortools
wget -O - <binary-distribution-url> | sudo tar -xzf - -C /opt/ortools --strip-components=1
cmake -DCMAKE_PREFIX_PATH=/opt/ortools ..
```
(when compiling your own programs, also pass `-I/opt/ortools/include -L/opt/ortools/lib -lortools`)
OR-Tools is found via the `ortoolsConfig.cmake` file it provides, so an OR-Tools installed in any standard location
(eg via a package manager, or built from source) also works.
