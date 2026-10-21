#!/bin/bash

set -x		# echo commands
set -e		# stop on any error

cat /proc/cpuinfo

# docs dependencies, already in images >= v15 (see setup.sh), installed here for older ones
/opt/python/cp311-cp311/bin/pip install Sphinx myst-parser furo numpy

# build qif
# Python3_EXECUTABLE/PYTHON_EXECUTABLE are needed only for docs (python module / sphinx-build), they need to match the version used to install
# Sphinx in the libqif-manylinux container (see setup.sh)
mkdir build
cd build
cmake -DMARCH=x86-64 -DPython3_EXECUTABLE=/opt/python/cp311-cp311/bin/python -DPYTHON_EXECUTABLE=/opt/python/cp311-cp311/bin/python ..		# the gh-actions vm crashes with march=native,sandybridge or haswell, so use x86-64 just for the tests
make qif_cpp tests_cpp samples docs -j 2
./tests_cpp/run

mkdir -p /output																	# this is copied to host's 'wheelhouse' dir by cibuildwheel. At this moment it doesn't exist yet, so create
mv misc/docs/_build/html /output													# save in case we need to publish it
touch /output/html/.nojekyll														# disable jekyll processing, cause it hides folders starting with underscore!

rm -rf *
cmake -DPython3_EXECUTABLE=/opt/python/cp311-cp311/bin/python -DPYTHON_EXECUTABLE=/opt/python/cp311-cp311/bin/python ..						# then compile again with the default march
make install -j 2