#!/usr/bin/perl -w
use strict;

use QIF::Channel;
use QIF::Compare;

# set this to 1 to use rational arithmetic (exlp is needed)
$QIF::Matrix::USE_RAT = 1;

my $C1 = QIF::Channel->new(matrix => q{
	2/3  1/3
    2/3   1/3
    1/4   3/4
});
my $C2 = QIF::Channel->new(matrix => q{
	1/2 1/2   0
	1/2   0 1/2
	  0 1/2 1/2
});

# sanity checks
die unless $C1->matrix->is_stochastic;
die unless $C2->matrix->is_stochastic;

# check whether C1 = C2 C3 for some C3
my $C3 = $C1->factorize($C2);
print "Factorizable:\n$C3"	if $C3;

# comparison. the various add_* methods allow to include different gain functions
# in the comparison
#
my $comp = QIF::Compare->new(c1 => $C1, c2 => $C2);
$comp->add_01_functions;
#$comp->add_2_block;
#$comp->add_identity;

my ($pi, $G) = $comp->compare;

if($pi) {
	print
		"C1 can leak more then C2\n\n" .
		"prior:\n$pi\n" .
		"gain function:\n$G\n" .
		"vulnerability of C1: " . $C1->g_vulnerability($pi, $G) .
		"\nvulnerability of C2: " . $C2->g_vulnerability($pi, $G) . "\n";
} else {
	print "C1 always leaks less than C2\n";
}


