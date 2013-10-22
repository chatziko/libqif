#!/bin/sh

# creating the libraries needeed by libqif
# ----------------------------------------------------------------------------------
echo "Installing armadillo"
cd lib/armadillo-3.920.1/
cmake .
make
sudo make install
cd ..
cd ..
# ----------------------------------------------------------------------------------
echo "Installing gtest"
cd lib/gtest-1.7.0/make/
make
cp gtest_main.a ../../gtest_main.a
cd ..
cd ..
cd ..
# ----------------------------------------------------------------------------------
echo "Installing glpk"
cd lib/glpk-4.52/
./configure
make
sudo make install
cd ..
cd ..
# ----------------------------------------------------------------------------------