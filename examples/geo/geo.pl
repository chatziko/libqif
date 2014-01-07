#!/usr/bin/perl -w
use strict;


use QIF::Channel;
use QIF::Compare;

$QIF::Matrix::USE_RAT = 0;


sub deterministic {
	my ($map) = @_;
	my @m;
	for my $val (@$map) {
		push @m, [ map { $val == $_ ? 1 : 0 } 0..@$map-1 ];
	}
	return QIF::Channel->new(matrix => QIF::Matrix->new(\@m));
}

sub euclid {
	my ($i, $j, $cols) = @_;
	my $xi = int($i / $cols);
	my $yi = $i % $cols;
	my $xj = int($j / $cols);
	my $yj = $j % $cols;
	return sqrt( ($xi - $xj)**2 + ($yi - $yj)**2 );
}

sub distance_gain {
	my ($rows, $cols) = @_;
	my $max = $rows*$cols-1;
	my $norm = euclid(0, $max, $cols);

	my @d;
	for my $i (0..$max) {
		for my $j (0..$max) {
			@d[$i][$j] = 1 - euclid($i, $j, $cols) / $norm;
		}
	}
	return QIF::Matrix->new(\@d);
}


my $g1 = deterministic([qw/
	0 1 2 1
/]);
my $d1 = distance_gain(2, 2, 1);



die $d1;
