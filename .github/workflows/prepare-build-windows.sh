#!/bin/bash

# Install gsl & gmp via vcpkg
/c/vcpkg/vcpkg install gsl:x64-windows
/c/vcpkg/vcpkg install gmp:x64-windows

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

# Install ortools
curl https://github.com/google/or-tools/releases/download/v9.15/or-tools_x64_VisualStudio2022_cpp_v9.15.6755.zip -L --output or-tools.zip
7z x or-tools.zip > /dev/null
cp -r or-tools_*/* /usr/local
popd

echo Ready to build the modules