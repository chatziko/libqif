# libqif

## Installation

### Prerequisites

* [CMake](http://www.cmake.org/) (tested with version 2.8.12)
* [Armadillo](http://arma.sourceforge.net/) (tested with version 4.4)
* [GMP](https://gmplib.org/) (tested with version 6.0.0)
* [GLPK](https://www.gnu.org/software/glpk/) (tested with version 4.54)
* A C++11 compliant compiler (tested with g++ 4.9.1)

On Ubuntu, these can be installed with:

    sudo apt-get install g++ cmake libarmadillo-dev libgmp-dev libglpk-dev

### Get the code

    git clone --recursive https://github.com/chatziko/libqif.git

Don't forget the `--recursive` to fetch the submodules.

### Compiling

    mkdir <path>/build && cd <path>/build
    cmake ..
    make

To run the tests

    make tests
    ./tests/run

You can also build the samples and doc

    make samples
    ./samples/<sample>

    make doc
    firefox ./misc/doxygen/html/index.html

### Repository structure: 

* `inc`: headers
* `src`: sources
* `tests`: test cases
* `samples`: sample programs
* `external`: external libs (googletest)

