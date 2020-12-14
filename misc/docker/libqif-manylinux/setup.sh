#!/bin/bash

# stop on any error
set -e

yum install -y wget dsl-devel gsl-devel xz


# cmake
yum remove -y cmake

cd /tmp
wget -nv https://github.com/Kitware/CMake/releases/download/v3.16.1/cmake-3.16.1-Linux-x86_64.tar.gz

tar -xf cmake-3.16.1-Linux-x86_64.tar.gz
cp -r cmake*/* /usr/local

# openblas
wget -nv https://github.com/xianyi/OpenBLAS/archive/v0.3.7.zip
unzip v0.3.7
cd OpenBLAS-0.3.7
make TARGET=SANDYBRIDGE
make install PREFIX=/usr/local
ln -s /usr/local/lib/libopenblas.so /usr/local/lib/liblapack.so		# openblas actually contains a lapack
ln -s /usr/local/lib/libopenblas.a /usr/local/lib/liblapack.a		# implementation, we just need symlinks to find it

# gmp
cd /tmp
wget -nv https://gmplib.org/download/gmp/gmp-6.1.2.tar.xz
tar -xf gmp*.tar.xz
cd gmp-6.1.2
./configure --build=sandybridge-pc-linux-gnu
make install

# armadillo
cd /tmp
wget -nv http://sourceforge.net/projects/arma/files/armadillo-9.800.3.tar.xz
tar -xf armadillo*.tar.xz
cd armadillo-9.800.3
./configure
make install

# ortools
# The C++ binary package of ortools is for CentOS 8, so we do a hack: we only copy the header files and static libs...
cd /tmp
wget -nv https://github.com/google/or-tools/releases/download/v8.1/or-tools_centos-8_v8.1.8487.tar.gz
tar -xf or-tools*.tar.gz
cp -r or-tools*/include /usr/local
cp or-tools*/lib/*.a /usr/local/lib

# ...and we extract the shared libraries from the python lib!
/opt/python/cp36-cp36m/bin/pip install ortools==8.1.8487
cp /opt/_internal/cpython-3.6.12/lib/python3.6/site-packages/ortools/.libs/* /usr/local/lib
ldconfig

# Sphinx needed for docs, numpy to import the qif module while generating docs
/opt/python/cp36-cp36m/bin/pip install Sphinx numpy

# cleanup
rm -rf /tmp/*
