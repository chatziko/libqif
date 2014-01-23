#!/usr/bin/perl -w
use strict;

use QIF::Channel;
use QIF::Compare;
use QIF::Compare::FindGain;
use QIF::Minvuln;

# set this to 1 to use rational arithmetic (exlp is needed)
$QIF::Matrix::USE_RAT = 1;

################################################
my $C = QIF::Channel->new(matrix => q{
    .9  .1  0
    .1  .7  .2 
    .1  .5  .4
});
die unless $C->matrix->is_stochastic;

# Note that gain functions need to be written as the
# *transposes* of the matrices shown in our CSF 2012 paper.
# (i.e. rows are secrets and columns are guesses)
my $g = QIF::Matrix->new(q{
    1 0 0 
    0 1 0
    0 0 1
});

# Using Minvuln, we can find that prior pi that minimizes 
# the posterior vulnerability

print "\nMinimizing vunerability ...\n\n";

my $comp = QIF::Minvuln->new(c1 => $C);
$comp->add_function($g);
my ($pi, $g) = $comp->compute;

if($pi) {
    my $additive = $C->g_vulnerability($pi, $g);
    print
        "The minimum g-vulnerability of\n\n" . $C .
        "\non gain function\n\n" . $g .
        "\nis " . $additive . ", realized on prior " . $pi . "\n"
} else {
     print "No prior found."
}
