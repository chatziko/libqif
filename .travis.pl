#!/usr/bin/perl -w
use strict;

my $linux = $ENV{TRAVIS_OS_NAME} eq 'linux';

sub run {
	my $cmd = shift;
	print "\n$cmd\n";
	system $cmd		and warn("command failed: $cmd\n\n"), exit 1;
}

# install dependencies
if($linux) {
	run qq{sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y};
	run qq{sudo apt-get -qq update};
	run qq{sudo apt-get install -y libgmp-dev libglpk-dev libgsl0-dev cmake g++-5 g++-6 g++-7 g++-8};

	run qq{git clone https://gitlab.com/conradsnicta/armadillo-code.git -b 7.800.x --depth 1};
	run qq{cd armadillo-code && ./configure && sudo make install};

	run qq{sudo apt-get -y install pkg-config build-essential autoconf libtool zlib1g-dev lsb-release}; # for ortools

} else {
	# on OSX we just install libqif via homebrew. This is useful to test by itself,
	# and also installs all dependencies needed for the actual build
	#
	run qq{rm -f /usr/local/include/c++};	# brew install will fail if this exists
	run qq{brew update};
	run qq{brew tap chatziko/tap};
	run qq{brew install --HEAD libqif};
	run qq{brew test --HEAD libqif};
}

# install ortools
run qq{git clone https://github.com/google/or-tools --depth 1};
run qq{cd or-tools && mkdir build && cd build && cmake .. && sudo cmake --build . --target install};

# build for each compiler
my @cxx = $linux ? qw/g++-5 g++-6 g++-7 g++-8/ : qw/clang++/;
for(@cxx) {
	run qq{mkdir -p build};
	chdir 'build';
	run qq{cmake -DCMAKE_CXX_COMPILER=$_ -DCMAKE_BUILD_TYPE=Release ..};
	run qq{make tests -j 2};
	run qq{./tests/run};
	chdir '..';
	run qq{rm -rf build};
}

exit 0;
