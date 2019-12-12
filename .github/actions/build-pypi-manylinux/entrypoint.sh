#!/bin/bash

# stop on any error
set -e

# build qif
mkdir build
cd build
cmake ..
make install


# build python
cd ../python

# source dist
/opt/python/cp36-cp36m/bin/python setup.py sdist

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
