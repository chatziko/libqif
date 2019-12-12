#!/bin/bash

# stop on any error
set -e

yum install -y wget openblas-devel lapack-devel dsl-devel gsl-devel xz


# cmake
yum remove -y cmake

cd /tmp
wget -nv https://github.com/Kitware/CMake/releases/download/v3.16.1/cmake-3.16.1-Linux-x86_64.tar.gz

tar -xf cmake-3.16.1-Linux-x86_64.tar.gz
cp -r cmake*/* /usr/local

# gmp
cd /tmp
wget -nv https://gmplib.org/download/gmp/gmp-6.1.2.tar.xz
tar -xf gmp*.tar.xz
cd gmp-6.1.2
./configure
make install

# armadillo
cd /tmp
wget -nv http://sourceforge.net/projects/arma/files/armadillo-9.800.3.tar.xz
tar -xf armadillo*.tar.xz
cd armadillo-9.800.3
./configure
make install

# ortools
# The C++ binary package of ortools is for CentOS 7, so we do a hack: we only copy the header files and static libs...
cd /tmp
wget -nv https://github.com/google/or-tools/releases/download/v7.4/or-tools_centos-7_v7.4.7247.tar.gz
tar -xf or-tools*.tar.gz
cp -r or-tools*/include /usr/local
cp or-tools*/lib/*.a /usr/local/lib

# ...and we extract the shared libraries from the python lib!
/opt/python/cp36-cp36m/bin/pip install ortools==7.4.7247
cp /opt/_internal/cpython-3.6.9/lib/python3.6/site-packages/ortools/.libs/* /usr/local/lib
ldconfig

# twine
/opt/python/cp36-cp36m/bin/pip install twine

# cleanup
rm -rf /tmp/*
