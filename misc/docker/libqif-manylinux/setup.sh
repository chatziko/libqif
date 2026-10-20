#!/bin/bash

# stop on any error
set -e

# tools and gcc 10 (the PATH is set in Dockerfile)
yum install -y \
    wget gsl-devel xz \
    \
    gcc-toolset-10-gcc gcc-toolset-10-gcc-c++ gcc-toolset-10-binutils \
    gcc-toolset-10-gcc-gfortran gcc-toolset-10-annobin \
    gcc-toolset-10-libstdc++-devel

# openblas
wget -nv https://github.com/OpenMathLib/OpenBLAS/archive/v0.3.34.zip
unzip v0.3.34
cd OpenBLAS-0.3.34
make -j$(nproc) DYNAMIC_ARCH=1									# kernels for all CPUs, selected at runtime
make install DYNAMIC_ARCH=1 PREFIX=/usr/local
ln -s /usr/local/lib/libopenblas.so /usr/local/lib/liblapack.so		# openblas actually contains a lapack
ln -s /usr/local/lib/libopenblas.a /usr/local/lib/liblapack.a		# implementation, we just need symlinks to find it

# gmp
cd /tmp
wget -nv https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz		# gmplib.org blocks CI runners
tar -xf gmp*.tar.xz
cd gmp-6.2.1
./configure --enable-fat										# code for all CPUs, selected at runtime
make -j$(nproc)
make install

# install an old cmake, cause the armadillo we use has an old CMakeLists.txt
curl -LO https://github.com/Kitware/CMake/releases/download/v3.18.6/cmake-3.18.6-Linux-x86_64.sh
chmod +x cmake-3.18.6-Linux-x86_64.sh
./cmake-3.18.6-Linux-x86_64.sh --skip-license --prefix=/usr/local

# ortools
cd /tmp
wget -nv 'https://github.com/google/or-tools/releases/download/v8.1/or-tools_centos-8_v8.1.8487.tar.gz'
tar -xf or-tools*.tar.gz
cp -r or-tools*/{lib,include} /usr/local

#### In the past we used this hack, not sure exactly why ####################
#
# The C++ binary package of ortools is for CentOS 8, so we do a hack: we only copy the header files and static libs...
# cd /tmp
# wget -nv 'https://github.com/google/or-tools/releases/download/v8.1/or-tools_centos-8_v8.1.8487.tar.gz'
# tar -xf or-tools*.tar.gz
# cp -r or-tools*/include /usr/local
# # cp or-tools*/lib/*.a /usr/local/lib

# # ...and we extract the shared libraries from the python lib!
# /opt/python/cp39-cp39/bin/pip install 'https://github.com/google/or-tools/releases/download/v8.1/ortools-8.1.8487-cp39-cp39-manylinux1_x86_64.whl'
# cp /opt/_internal/cpython-3.9.23/lib/python3.9/site-packages/ortools/.libs/* /usr/local/lib
#
##############################################################################


ldconfig

# Sphinx needed for docs, numpy to import the qif module while generating docs
/opt/python/cp39-cp39/bin/pip install Sphinx numpy

# cleanup
rm -rf /tmp/*
