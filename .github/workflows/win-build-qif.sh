#!/bin/bash

# Install gsl & mpir (GMP alternative for windows) via vcpkg
/c/vcpkg/vcpkg install gsl:x64-windows
/c/vcpkg/vcpkg install mpir:x64-windows

# Install openblas
# Note: armadillo searches it as 'openblas' but in the resulting wheel it needs to be 'libopenblas', so we create both
mkdir -p /usr/local/lib /usr/local/include
pushd /tmp
curl https://github.com/xianyi/OpenBLAS/releases/download/v0.3.10/OpenBLAS-0.3.10-x64.zip -L --output openblas.zip
(mkdir openblas && cd openblas && 7z x ../openblas.zip)
cp openblas/bin/libopenblas.dll /usr/local/lib/openblas.dll
cp openblas/bin/libopenblas.dll /usr/local/lib/libopenblas.dll
cp openblas/lib/libopenblas.dll.a /usr/local/lib/openblas.lib
cp openblas/lib/libopenblas.dll.a /usr/local/lib/libopenblas.lib

# Install ortols
curl https://github.com/google/or-tools/releases/download/v8.1/or-tools_VisualStudio2019-64bit_v8.1.8487.zip -L --output or-tools.zip
7z x or-tools.zip
cp -r or-tools*/include /usr/local
cp or-tools*/lib/*.lib /usr/local/lib
popd

# Install depedencies of repair-wheel-windows.py
pip install pefile machomachomangler

echo Ready to build the modules