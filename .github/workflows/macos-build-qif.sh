#!/bin/bash

export HOMEBREW_NO_INSTALL_CLEANUP=1	# make homebrew
export HOMEBREW_NO_AUTO_UPDATE=1		# faster

cmake --version

rm -f /usr/local/include/c++	# brew install will fail if this exists

brew install gsl gmp armadillo

# or-tools
(brew unlink protobuf; true)	# unlink protobuf if present so that it doesn't conflict with the one installed by or-tools
brew tap chatziko/tap
brew install --HEAD or-tools

# build qif
mkdir build
cd build
cmake ..
make install -j 2
