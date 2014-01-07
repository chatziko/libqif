# Here are some examples of what can be done with Kostas's tools
# for analyzing g-leakage using linear programming techniques.

#!/usr/bin/perl -w
use strict;

use QIF::Channel;
use QIF::Compare;
use QIF::Compare::FindGain;

# set this to 1 to use rational arithmetic (exlp is needed)
$QIF::Matrix::USE_RAT = 1;

################################################
# Ex4 from our CSF 2012 paper:

my $Ex4 = QIF::Channel->new(matrix => q{
    1/2 1/2
    1   0
    0   1
});
die unless $Ex4->matrix->is_stochastic;

my $pi = QIF::Matrix->uniform(3);

# Note that gain functions are written as the *transposes* of the
# matrices shown in our CSF 2012 paper.
my $G_d = QIF::Matrix->new(q{
    1   0    0
    0   1    .98
    0   .98  1
});

print "V_g_d(pi,Ex4)  = " . $Ex4->g_vulnerability($pi, $G_d) . "\n\n";

################################################
# Ex5 from our CSF 2012 paper:

my $Ex5 = QIF::Channel->new(matrix => q{
    .6  .4 
    0   1
    0   1
});
die unless $Ex5->matrix->is_stochastic;

$pi = QIF::Matrix->new('.6 .2 .2');

print "V(pi,Ex5)      = " . $Ex5->vulnerability($pi) . "\n";
print "V_g_d(pi,Ex5)  = " . $Ex5->g_vulnerability($pi, $G_d) . "\n\n";

################################################
# The imperfect cancer test from my QEST 2011 paper:
my $Cancer = QIF::Channel->new(matrix => q{
    .90 .10 
    .07 .93
});
die unless $Cancer->matrix->is_stochastic;

$pi = QIF::Matrix->new('1/125 124/125');   # Equivalent to (.008, .992)

print "V(pi,Cancer)   = " . $Cancer->vulnerability($pi) . "\n\n";

# Suppose we want a gain function G that maximizes the G-vulnerability
# of Cancer under prior pi (which is not uniform).
# FindGain can do this! We compare Cancer with NI (which leaks nothing). 

my $NI = QIF::Channel->new(matrix => q{
    1 
    1
});
die unless $NI->matrix->is_stochastic;

my $comp = QIF::Compare::FindGain->new(c1 => $Cancer, c2 => $NI, prior => $pi);

my $G = $comp->compare(2);
# 2 is the size of the set W of possible guesses.
# (Increasing it seems to make no difference.)

if($G) {
    print
        "Cancer can outleak NI on prior " . $pi . "\n" .
        "gain function G:\n$G\n" .
        "V_G(pi,Cancer) = " . $Cancer->g_vulnerability($pi, $G) . "\n" .
        "V_G(pi,NI)     = " . $NI->g_vulnerability($pi, $G) . "\n\n";
} else {
        print "Cancer never outleaks NI on prior " . $pi . "\n\n";
}

##################################################
#The channels from section VI.C of our CSF 2012 paper:
my $Test1 = QIF::Channel->new(matrix => q{
        .2  .22 .58
        .2  .4  .4
        .35 .4  .25
});
die unless $Test1->matrix->is_stochastic;

my $Test2 = QIF::Channel->new(matrix => q{
        .1 .4 .1 .4
        .2 .2 .3 .3
        .5 .1 .1 .3
});
die unless $Test2->matrix->is_stochastic;

my $comp = QIF::Compare::FindGain->new(c1 => $Test1, c2 => $Test2);

my $G = $comp->compare(3);
# 3 is the size of the set W of possible guesses.
# (Increasing it seems to make no difference.)

my $pi = QIF::Matrix->uniform(3);

if($G) {
    print
        "Test1 can outleak Test2 on uniform prior and  gain function G:\n$G\n" .
        "V_G(pi,Test1) = " . $Test1->g_vulnerability($pi, $G) . "\n" .
        "V_G(pi,Test2) = " . $Test2->g_vulnerability($pi, $G) . "\n\n";
} else {
        print "Test1 never outleaks Test2 on uniform prior.\n\n";
}

exit;

