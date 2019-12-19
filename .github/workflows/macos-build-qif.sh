#!/bin/bash

# stop on any error
set -e

export HOMEBREW_NO_INSTALL_CLEANUP=1	# make homebrew
export HOMEBREW_NO_AUTO_UPDATE=1		# faster

cmake --version

rm -f /usr/local/include/c++	# brew install will fail if this exists

brew install gsl gmp

# armadillo (manual install to use Accelerate, homebrew version uses openblas)
pushd /tmp
wget -nv http://sourceforge.net/projects/arma/files/armadillo-9.800.3.tar.xz
tar -xf armadillo*.tar.xz
cd armadillo-9.800.3
./configure
make install
popd

# or-tools
brew unlink protobuf || true	# unlink protobuf if present so that it doesn't conflict with the one installed by or-tools
brew tap chatziko/tap
brew install --HEAD or-tools

# build qif
mkdir build
cd build
cmake ..
make qif tests samples -j 2
./tests/run
make install -j 2
