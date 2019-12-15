#!/bin/bash

# stop on any error
set -e

# build qif
mkdir build
cd build
cmake -DPYTHON_EXECUTABLE=/opt/python/cp36-cp36m/bin/python ..		# python needed just for the config, cause default is python2
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
