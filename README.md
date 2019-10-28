# libqif

[![Build Status](https://travis-ci.org/chatziko/libqif.svg?branch=master)](https://travis-ci.org/chatziko/libqif)

## Install via Homebrew

The easiest way to install libqif (especially on OSX) is via [Homebrew](http://brew.sh/):
```bash
brew tap chatziko/tap
brew install --HEAD libqif
brew test --HEAD libqif
```

To upgrade to the latest version:
```bash
brew reinstall --HEAD libqif
```

__On Linux__: the method above is available via [Linuxbrew](http://linuxbrew.sh/). Make sure you have ```csh``` installed (needed only for the installation). Also, Linuxbrew installs everything under ```$HOME/.linuxbrew```, so you need to either configure your system to use libraries from there, or symlink eveything under ```/usr/local```:
```bash
sudo ln -s $HOME/.linuxbrew/include/qif* /usr/local/include/
sudo ln -s $HOME/.linuxbrew/lib64/libqif* /usr/local/lib/
sudo ldconfig
```

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
        << bayes_vuln::posterior(pi, C) << "\n";
}
```

Compile and run with:
```bash
g++ test.cpp -std=c++17 -lqif -larmadillo -o test
./test

# or with clang
clang++ test.cpp -std=c++17 -lqif -larmadillo -o test
```

You can find more sample programs in the [samples](https://github.com/chatziko/libqif/tree/master/samples) directory.

If the `rat` type is used you also need to compile with `-lgmp -lgmpxx`.
If OR-Tools are used you also need to compile with `-lortools`.

## Build libqif from source

Prerequisites

* [CMake](http://www.cmake.org/)
* [Armadillo](http://arma.sourceforge.net/)
* [GMP](https://gmplib.org/)
* [GLPK](https://www.gnu.org/software/glpk/)
* [GSL](http://www.gnu.org/software/gsl/)
* A C++17 compliant compiler (eg g++ or clang)

Optionally

* [OR-Tools](https://developers.google.com/optimization/)

On Ubuntu, these can be installed with:
```bash
sudo apt-get install g++ cmake libarmadillo-dev libgmp-dev libglpk-dev libgsl0-dev
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

To run the tests
```bash
make tests
./tests/run
```

To build the samples:
```bash
make samples
./samples/<sample>
```

#### Use OR-Tools

If OR-Tools are installed in your system when libqif is compiled, they will
be used by some libqif methods for linear optimization and network flow.
The easiest way is to locate the `.tar.gz` file for your system in the
[OR-Tools binary distributions](https://developers.google.com/optimization/install/cpp/#binary-distributions),
then install it by simply copying the libraries and header files under `/usr/local` :
```
wget <binary-distribution-url> | tar -xzf -
sudo cp -r ortools*/{lib,include} /usr/local/
```

Alternatively, you can
[install OR-Tools from source via cmake](https://github.com/google/or-tools/blob/stable/cmake/README.md#building-or-tools-with-cmake).
The instructions in that link are dated, a method that has been tested is below:
```
git clone https://github.com/google/or-tools --depth 1 --branch=v7.2
mkdir or-tools/build && cd or-tools/build

cmake -DBUILD_DEPS:BOOL=ON ..
sudo cmake --build . --target install
sudo cp -r dependencies/install/* /usr/local
```

### Repository structure:

* `inc`: headers
* `src`: sources
* `tests`: test cases
* `samples`: sample programs
* `external`: external libs (googletest)
