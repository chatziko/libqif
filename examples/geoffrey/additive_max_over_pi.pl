# Examples of maximizing additive g-leakage using linear programming

#!/usr/bin/perl -w
use strict;

use QIF::Channel;
use QIF::Compare;
use QIF::Compare::FindGain;

# set this to 1 to use rational arithmetic (exlp is needed)
$QIF::Matrix::USE_RAT = 1;

################################################
my $C = QIF::Channel->new(matrix => q{
    .9  .1  0
    .1  .7  .2 
    .1  .5  .4
});
die unless $C->matrix->is_stochastic;

# A channel satisfying noninterference
my $NI = QIF::Channel->new(matrix => q{
    1 
    1
    1
});
die unless $NI->matrix->is_stochastic;

# Note that gain functions need to be written as the
# *transposes* of the matrices shown in our CSF 2012 paper.
# (i.e. rows are secrets and columns are guesses)
my $g = QIF::Matrix->new(q{
    .7  .1 
    .1  .4
    .4  .9
});

# Using Compare, we can find that prior pi that maximizes
#   V_g(pi,C) - V_g(pi,NI)
# which is the additive leakage

my $comp = QIF::Compare->new(c1 => $C, c2 => $NI);
$comp->add_function($g);
my ($pi, $g) = $comp->compare;

if($pi) {
    my $additive = $C->g_vulnerability($pi, $g) - $NI->g_vulnerability($pi, $g);
    print
        "The maximum additive leakage of \n\n" . $C .
        "\nunder gain function \n\n" . $g .
        "\nis " . $additive . ",\n" .
        "realized on prior " . $pi . "\n"
} else {
        print "The maximum additive leakage is 0.\n";
}

$comp = QIF::Compare->new(c1 => $C, c2 => $NI);
$comp->add_identity;
my ($pi, $g) = $comp->compare;

if($pi) {
    my $additive = $C->g_vulnerability($pi, $g) - $NI->g_vulnerability($pi, $g);
    print
        "The maximum additive leakage of \n\n" . $C .
        "\nunder gain function \n\n" . $g .
        "\nis " . $additive . ",\n" .
        "realized on prior " . $pi . "\n"
} else {
        print "The maximum additive leakage is 0.\n";
}

exit;

