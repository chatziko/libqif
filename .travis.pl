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
	run qq{sudo apt-get install -y libgmp-dev libglpk-dev libgsl0-dev cmake g++-7 g++-8 g++-9 clang};

	run qq{git clone https://gitlab.com/conradsnicta/armadillo-code.git -b 8.400.x --depth 1};
	run qq{cd armadillo-code && ./configure && sudo make install};

} else {
	# on OSX we just install libqif via homebrew. This is useful to test by itself,
	# and also installs all dependencies needed for the actual build
	#
	$ENV{HOMEBREW_NO_INSTALL_CLEANUP} = 1;	# make homebrew
	$ENV{HOMEBREW_NO_AUTO_UPDATE} = 1;		# faster

	run qq{rm -f /usr/local/include/c++};	# brew install will fail if this exists
	# run qq{brew update};
	run qq{brew tap chatziko/tap};
	run qq{brew install --HEAD libqif};
	run qq{brew test --HEAD libqif};
}

# install ortools
my $ortools_url = $linux
	? 'https://github.com/google/or-tools/releases/download/v7.2/or-tools_ubuntu-16.04_v7.2.6977.tar.gz'
	: 'https://github.com/google/or-tools/releases/download/v7.2/or-tools_MacOsX-10.14.5_v7.2.6977.tar.gz';
run qq{wget -O - '$ortools_url' | tar -xzf -; mv or-tools* or-tools};
run qq{sudo cp -r or-tools/lib or-tools/include /usr/local/};

# build for each compiler
my @cxx = $linux ? qw/g++-7 g++-8 g++-9 clang++/ : qw/clang++/;
for(@cxx) {
	run qq{mkdir -p build};
	chdir 'build';
	run qq{cmake -DCMAKE_CXX_COMPILER=$_ -DCMAKE_BUILD_TYPE=Release ..};
	run qq{make allcode -j 2};
	run qq{./tests/run};
	chdir '..';
	run qq{rm -rf build};
}

exit 0;
