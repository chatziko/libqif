#!/usr/bin/perl -w
use strict;

use QIF::Channel;
use QIF::Compare;

# set this to 1 to use rational arithmetic (exlp is needed)
$QIF::Matrix::USE_RAT = 1;

# Let's define some channels:

my $Ex4 = QIF::Channel->new(matrix => q{
	1/2 1/2
        1   0
        0   1
});
die unless $Ex4->matrix->is_stochastic;

my $Ex5 = QIF::Channel->new(matrix => q{
	.6  .4 
        0   1
        0   1
});
die unless $Ex5->matrix->is_stochastic;

# Let's define some gain functions:

my $G_d = QIF::Matrix->new(q{
	1   0       0
        0   1       .98
        0   .98     1
});

my $G_id = QIF::Matrix->new(q{
	1 0 0
	0 1 0
        0 0 1
});

# Let's define some priors:

my $pi = QIF::Matrix->uniform(3);
# I don't know how to create a prior like (.6, .2, .2)!
#$pi= QIF::Matrix->new('.6 .2 .2');
print $pi;

print "V_g_d(pi,Ex4)  = " .  $Ex4->g_vulnerability($pi, $G_d) . "\n";
print "V(pi,Ex5)      = " .  $Ex5->vulnerability($pi) . "\n";
print "V_g_id(pi,Ex5) = " .  $Ex5->g_vulnerability($pi, $G_id) . "\n";
print "V_g_d(pi,Ex5)  = " .  $Ex5->g_vulnerability($pi, $G_d) . "\n";

exit;

