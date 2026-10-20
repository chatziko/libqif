#!/bin/bash

# stop on any error
set -e

sysctl -a | grep machdep.cpu			# cpu info

export HOMEBREW_NO_INSTALL_CLEANUP=1	# make homebrew
export HOMEBREW_NO_AUTO_UPDATE=1		# faster

cmake --version

rm -rf /usr/local/include/c++	# brew install will fail if this exists

brew install gsl gmp

# or-tools, in /opt/ortools (CMAKE_PREFIX_PATH is set in build.yml)
brew unlink protobuf abseil || true	# unlink if present so that they don't conflict with the ones included in or-tools
curl -L https://github.com/google/or-tools/releases/download/v9.15/or-tools_$(uname -m)_macOS-26.2_cpp_v9.15.6755.tar.gz --output or-tools.tar.gz
export CMAKE_PREFIX_PATH=/opt/ortools
sudo mkdir -p /opt/ortools
sudo tar -xzf or-tools.tar.gz -C /opt/ortools --strip-components=1
rm or-tools.tar.gz
for lib in /opt/ortools/lib/*.dylib; do	# ortools' dylibs load each other via @rpath, but have no rpath (needed by delocate)
	[ -L "$lib" ] && continue
	sudo install_name_tool -add_rpath @loader_path "$lib" 2>/dev/null || true
	sudo codesign --force -s - "$lib"
done

# build qif
mkdir build
cd build
# the gh-actions vm crashes with march=native,sandybridge or haswell, so use x86-64 just for the tests
cmake $([ $(uname -m) = x86_64 ] && echo -DMARCH=x86-64) ..
make qif_cpp tests_cpp samples -j 2
./tests_cpp/run

rm -rf *
cmake ..						# then compile again with the default march
sudo make install -j 2			# /usr/local is not writable on arm64 runners (homebrew is in /opt/homebrew)
