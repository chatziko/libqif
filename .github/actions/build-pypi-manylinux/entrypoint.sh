#!/bin/bash

# stop on any error
set -e

cat /proc/cpuinfo

# Sphinx needed for docs, numpy to import the qif module while generating docs
/opt/python/cp36-cp36m/bin/pip install Sphinx numpy

# build qif
mkdir build
cd build
cmake -DMARCH=x86-64 -DPYTHON_EXECUTABLE=/opt/python/cp36-cp36m/bin/python ..		# the gh-actions vm crashes with march=native,sandybridge or haswell, so use x86-64 just for the tests
make qif tests samples docs -j 2
./tests/run

rm -rf *
cmake -DPYTHON_EXECUTABLE=/opt/python/cp36-cp36m/bin/python ..						# then compile again with the default march
make install -j 2


# build python
cd ../python_lib

# binary dists for each python version
for pyver in /opt/python/cp3{5,6,7,8}*
do
	rm -rf build
	$pyver/bin/python setup.py bdist_wheel
done

# repair wheels, then move under dist/
for wheel in dist/*.whl
do
	auditwheel repair $wheel
done

rm dist/*.whl
mv wheelhouse/*.whl dist/
