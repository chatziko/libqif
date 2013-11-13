#!/bin/sh
clear
echo "testing the LIBQIF library"

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