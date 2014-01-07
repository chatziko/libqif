#!/usr/bin/perl -w
use strict;

# Usage:
#   perl run.pl [--maximize=<ratio-sv|ratio-hr|prob-error>] [--fork]
#
# This script computes the set of corner points of the probability of error for a given
# channel. Using the corner points the script can compute 3 things (depending on the
# --maximize parameter)
#    - The maximum probability of error (--maximize=prob-error)
#    - The maximum ratio prob-error/Santhi-Vardy bound (--maximize=ratio-sv)
#    - The maximum ratio prob-error/Hellman-Raviv bound (--maximize=ratio-hr)
#
# NOTE: The script currently doesn't work for matrices with zero elements (this should be easy to fix)


# The channel matrix to use. Note that the complexity of the algorithm is exponential so
# the matrix has to be small!
#
# The following one comes from an instance of crowds
#
my $matrix = q{
	0.32961  0.32961  0.17039  0.17039
	0.32961  0.17039  0.32961  0.17039

	0.68482  0.07240  0.07240  0.17039
	0.07240  0.68482  0.17039  0.07240
	0.07240  0.17039  0.68482  0.07240
	0.17039  0.07240  0.07240  0.68482
};

# For any matrix with capacity 0 (all rows are the same)
# the maximum ratio must be 1 and the maximum
# probability of error must be n-1/n
#
$matrix = q{
	.5 .3 .1 .1
	.5 .3 .1 .1
	.5 .3 .1 .1
	.5 .3 .1 .1
};


use Getopt::Long;

use QIF::Channel;
use QIF::Corners;

my $rat = 1;

my ($maximize, $fork) = ("ratio-sv");
GetOptions(
	"fork"			=> \$fork,
	"maximize=s"	=> \$maximize,
)	or exit;

my %maximizeValues = (
	'ratio-sv' => 1,
	'ratio-hr' => 2,
	'prob-error' => 3,
);
exists $maximizeValues{$maximize}
	or die "Valid values for 'maximize' option are: ratio-sv, ratio-hr, prob-error";



my $chan = QIF::Channel->new(matrix => $matrix);
my $corn = QIF::Corners->new($chan);
$corn->{'fork'} = $fork;

my $res = $corn->best_ratio($maximizeValues{$maximize});

print "Maximum $maximize: $res->{ratio}\n\n";
print "Distribution giving the maximum: ", join(', ', @{$res->{x}}), "\n\n";
print "Solved systems: $res->{count}\n";

