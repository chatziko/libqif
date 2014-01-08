# Here are some examples of what can be done with Kostas's tools
# for analyzing g-leakage using linear programming techniques.

#!/usr/bin/perl -w
use strict;

use QIF::Channel;
use QIF::Compare;
use QIF::Compare::FindGain;

# set this to 1 to use rational arithmetic (exlp is needed)
$QIF::Matrix::USE_RAT = 0;



my $n = 4;
my $Id = QIF::Channel->new(matrix => QIF::Matrix->identity($n));
die unless $Id->matrix->is_stochastic;

# Suppose we want a gain function G that maximizes the G-vulnerability
# of Id under prior pi (which is not uniform).
# FindGain can do this! We compare Id with NI (which leaks nothing). 

my $NI = QIF::Channel->new(matrix => QIF::Matrix->new(
	[([1]) x $n])
);
die unless $NI->matrix->is_stochastic;


print "aa:". max_leakage(QIF::Matrix->new('1/4 1/4 1/4 1/4'));
exit;

my ($max, $pi) = QIF::Channel::max_over_priors($n, 0.1, &max_leakage);
print "max: $max\n pi: $pi";


sub max_leakage {
	my ($pi) = @_;
	my $comp = QIF::Compare::FindGain->new(c1 => $Id, c2 => $NI, prior => $pi);

	my $G = $comp->compare($n)		or return 0;
	# 2 is the size of the set W of possible guesses.
	# (Increasing it seems to make no difference.)

	my $leak = ($Id->g_vulnerability($pi, $G) - $NI->g_vulnerability($pi, $G));
	print "-- $pi" if $leak > 0.599;
	return $leak;
}
