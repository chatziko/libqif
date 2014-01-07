#!/usr/bin/perl -w
use strict;

use POSIX qw/ceil/;
use List::Util qw/max/;

use QIF::Channel;
use QIF::Graph;
use QIF::LinearProgram;
use QIF::Mechanism;


my $prog = QIF::LinearProgram->new(
	equalities => QIF::Matrix->new(''),
	inequalities => QIF::Matrix->new('
		0 -1 0
		0 -1 0
		120 210 -15000
		110 30 -4000
		1 1 -75
	'),
	objective => QIF::Matrix->new('
		143 60 0
	'),
);

use Data::Dumper;
print Dumper [$prog->solve];

exit;

my $C1 = QIF::Channel->new(matrix => q{
	1 0
	0 1
	0 1
});
my $C2 = QIF::Channel->new(matrix => q{
	0 1
	1 0
	0 1
});
my $C3 = QIF::Channel->new(matrix => q{
	0 1
	0 1
	1 0
});
my $prior = [0.6, 0.3, 0.1];

warn $C1->vulnerability($prior);
warn $C2->vulnerability($prior);
warn $C3->vulnerability($prior);
warn $C1->mutual_information($prior);
warn $C2->mutual_information($prior);
warn $C3->mutual_information($prior);

