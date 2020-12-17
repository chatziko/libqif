#!/bin/bash

# the source tree contains two symlinks to directories. In windows, directory symlinks are
# different than file symlinks, and git creates file symlinks. So we delete them, and create
# directory symlinks instead.
pushd lib/include/qif_bits
rm miniball osqp
cmd.exe /c 'mklink /D miniball ..\..\..\external\miniball\cpp\main'
cmd.exe /c 'mklink /D osqp ..\..\..\external\osqp\include'
popd

# Install gsl & mpir (GMP alternative for windows) via vcpkg
/c/vcpkg/vcpkg install gsl:x64-windows
/c/vcpkg/vcpkg install mpir:x64-windows
# /c/vcpkg/vcpkg integrate install

# Install openblas
mkdir -p /usr/local/lib /usr/local/include
pushd /tmp
curl https://github.com/xianyi/OpenBLAS/releases/download/v0.3.10/OpenBLAS-0.3.10-x64.zip -L --output openblas.zip
(mkdir openblas && cd openblas && 7z x ../openblas.zip)
cp openblas/bin/libopenblas.dll /usr/local/lib/openblas.dll
cp openblas/bin/libopenblas.dll /usr/local/lib/libopenblas.dll
cp openblas/lib/libopenblas.dll.a /usr/local/lib/openblas.lib

# Install ortols
curl https://github.com/google/or-tools/releases/download/v8.1/or-tools_VisualStudio2019-64bit_v8.1.8487.zip -L --output or-tools.zip
7z x or-tools.zip
cp -r or-tools*/include /usr/local
cp or-tools*/lib/*.lib /usr/local/lib
popd

echo Ready to build the modules