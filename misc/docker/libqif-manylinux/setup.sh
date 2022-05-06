#!/bin/bash

# stop on any error
set -e

yum install -y wget gsl-devel xz


# openblas
wget -nv https://github.com/xianyi/OpenBLAS/archive/v0.3.20.zip
unzip v0.3.20
cd OpenBLAS-0.3.20
make TARGET=SANDYBRIDGE
make install PREFIX=/usr/local
ln -s /usr/local/lib/libopenblas.so /usr/local/lib/liblapack.so		# openblas actually contains a lapack
ln -s /usr/local/lib/libopenblas.a /usr/local/lib/liblapack.a		# implementation, we just need symlinks to find it

# gmp
cd /tmp
wget -nv https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz
tar -xf gmp*.tar.xz
cd gmp-6.2.1
./configure --build=sandybridge-pc-linux-gnu
make install

# ortools
# The C++ binary package of ortools is for CentOS 8, so we do a hack: we only copy the header files and static libs...
cd /tmp
wget -nv https://github.com/google/or-tools/releases/download/v8.1/or-tools_centos-8_v8.1.8487.tar.gz
tar -xf or-tools*.tar.gz
cp -r or-tools*/include /usr/local
cp or-tools*/lib/*.a /usr/local/lib

# ...and we extract the shared libraries from the python lib!
/opt/python/cp37-cp37m/bin/pip install ortools==8.1.8487
cp /opt/_internal/cpython-3.7.13/lib/python3.7/site-packages/ortools/.libs/* /usr/local/lib
ldconfig

# Sphinx needed for docs, numpy to import the qif module while generating docs
/opt/python/cp37-cp37m/bin/pip install Sphinx numpy

# cleanup
rm -rf /tmp/*
