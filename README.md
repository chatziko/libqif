# libqif

[![Build Status](https://travis-ci.org/chatziko/libqif.svg?branch=master)](https://travis-ci.org/chatziko/libqif)

## Install via Homebrew

The easiest way to install libqif (especially on OSX) is via [Homebrew](http://brew.sh/):

    brew tap chatziko/tap
    brew install --HEAD libqif
    brew test --HEAD libqif

__On Linux__: the method above is available via [Linuxbrew](http://linuxbrew.sh/). Make sure you have ```csh``` installed (needed only for the installation). Also, Linuxbrew installs everything under ```$HOME/.linuxbrew```, so you need to either configure your system to use libraries from there, or symlink eveything under ```/usr/local```:

    sudo ln -s $HOME/.linuxbrew/include/qif* /usr/local/include/
    sudo ln -s $HOME/.linuxbrew/lib64/libqif* /usr/local/lib/
    sudo ldconfig

## Compile a sample program

Create a simple ```test.cpp``` file:

    #include <qif>
    using namespace qif;

    int main() {
        chan C("1 0 0; 0 1 0; 0 0 1");
        prob pi = probab::uniform<double>(3);
        std::cout
            << "Bayes vulnerability of \n"
            << C << " under " << pi << " is "
            << bayes::post_vulnerability(pi, C) << "\n";
    }

Compile and run with:

    g++ test.cpp -std=c++11 -lqif -larmadillo -o test
    ./test
 
    # or with clang
    clang++ test.cpp -std=c++11 -lqif -larmadillo -o test

You can find more sample programs in the [samples](https://github.com/chatziko/libqif/tree/master/samples) directory.

If the `rat` type is used you also need to compile with `-lgmp -lgmpxx`.

## Build libqif from source

Prerequisites

* [CMake](http://www.cmake.org/)
* [Armadillo](http://arma.sourceforge.net/)
* [GMP](https://gmplib.org/)
* [GLPK](https://www.gnu.org/software/glpk/)
* [GSL](http://www.gnu.org/software/gsl/)
* A C++11 compliant compiler (eg g++ or clang)

On Ubuntu, these can be installed with:

    sudo apt-get install g++ cmake libarmadillo-dev libgmp-dev libglpk-dev libgsl0-dev

Get the code (note the `--recursive` to fetch the submodules).

    git clone --recursive https://github.com/chatziko/libqif.git

To compile / install:

    mkdir <path>/build && cd <path>/build
    cmake ..
    make
    sudo make install

To run the tests

    make tests
    ./tests/run

To build the samples:

    make samples
    ./samples/<sample>

### Repository structure:

* `inc`: headers
* `src`: sources
* `tests`: test cases
* `samples`: sample programs
* `external`: external libs (googletest)
