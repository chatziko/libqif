#!/bin/sh

clear
echo "Compiling LIBQIF Library"
#creating the "*.o" files ---------------------------------------------------------

#--- that needs armadillo
g++ -I./inc -c src/Channel.cpp
g++ -I./inc -c src/Graph.cpp
g++ -I./inc -c src/Gain.cpp
g++ -I./inc -c src/Prob.cpp
g++ -I./inc -c src/Mechanism.cpp

#--- that needs the basic things and scilab_call
g++ -I./inc -c src/EntropyModel.cpp
g++ -I./inc -c src/Shannon.cpp
g++ -I./inc -c src/MinEntropy.cpp
g++ -I./inc -c src/Guessing.cpp

#--- that needs glpk
g++ -I./inc -c src/LinearProgram.cpp -lglpk

#--- that needs LinearProgram
g++ -I./inc -c src/GLeakage.cpp -lglpk

#creating the library---------------------------------------------------------------

ar -rvs libqif.a Channel.o inc/Channel.h Graph.o inc/Graph.h Gain.o inc/Gain.h Prob.o inc/Prob.h Mechanism.o inc/Mechanism.h EntropyModel.o inc/EntropyModel.h Shannon.o inc/Shannon.h MinEntropy.o inc/MinEntropy.h Guessing.o inc/Guessing.h LinearProgram.o inc/LinearProgram.h GLeakage.o inc/GLeakage.h

mv libqif.a lib/

# Documentation --------------------------------------------------------------------
echo "Creating Library Documentation"

# preparing doxygen files
cd doxygen/
mv Doxyfile footer.html init stylesheet.css mainpage.dox header.html ../
cd ..
#requires doxygen
# doxygen Doxyfile

# restoring the doxygen files
mv Doxyfile footer.html init stylesheet.css mainpage.dox header.html doxygen/

# TESTING ---------------------------------------------------------------------------
echo "testing the LIBQIF library"

# creating .o files
g++ -I./inc -c tests/ChannelTest.cc
g++ -I./inc -c tests/GraphTest.cc
g++ -I./inc -c tests/GainTest.cc
g++ -I./inc -c tests/ProbTest.cc
g++ -I./inc -c tests/MechanismTest.cc
g++ -I./inc -c tests/EntropyModelTest.cc
g++ -I./inc -c tests/ShannonTest.cc
g++ -I./inc -c tests/MinEntropyTest.cc
g++ -I./inc -c tests/GuessingTest.cc
g++ -I./inc -c tests/LinearProgramTest.cc
g++ -I./inc -c tests/GLeakageTest.cc

# linking with gtest
g++ -I./inc Channel.o ChannelTest.o lib/gtest_main.a -lpthread -o tests/ChannelTest
g++ -I./inc Graph.o GraphTest.o lib/gtest_main.a -lpthread -o tests/GraphTest
g++ -I./inc Gain.o GainTest.o lib/gtest_main.a -lpthread -o tests/GainTest
g++ -I./inc Prob.o ProbTest.o lib/gtest_main.a -lpthread -o tests/ProbTest
g++ -I./inc Mechanism.o MechanismTest.o lib/gtest_main.a lib/libqif.a -lpthread -o tests/MechanismTest
g++ -I./inc EntropyModel.o EntropyModelTest.o lib/gtest_main.a -lpthread -o tests/EntropyModelTest
g++ -I./inc Shannon.o ShannonTest.o lib/gtest_main.a lib/libqif.a -lpthread -o tests/ShannonTest
g++ -I./inc MinEntropy.o MinEntropyTest.o lib/gtest_main.a lib/libqif.a -lpthread -o tests/MinEntropyTest
g++ -I./inc Guessing.o GuessingTest.o lib/gtest_main.a lib/libqif.a -lpthread -o tests/GuessingTest
g++ -I./inc LinearProgram.o LinearProgramTest.o lib/gtest_main.a -lglpk -lpthread -o tests/LinearProgramTest
g++ -I./inc GLeakage.o GLeakageTest.o lib/gtest_main.a lib/libqif.a -lpthread -o tests/GLeakageTest
rm *.o

# running test cases
cd tests/
 ./ChannelTest
 ./GraphTest
 ./GainTest
 ./ProbTest
 ./MechanismTest
 ./EntropyModelTest
 ./ShannonTest
 ./MinEntropyTest
 ./GuessingTest
 ./LinearProgramTest
 ./GLeakageTest
cd ..

# ----------------------------------------------------------------------------------
echo "LIBQIF Library is in bin/libqif.a"
echo "LIBQIF Library Documentation is in html/index.html"
echo "LIBQIF example is in samples/"
g++ -I./inc -L./ -o samples/programa1 samples/Main.cpp lib/libqif.a