#!/usr/bin/perl -w
use strict;


use QIF::Channel;

$QIF::Matrix::USE_RAT = 0;


while(1) {
	my $C = QIF::Channel->random(2,2)->matrix;
	my $p = QIF::Channel->random(1,2)->matrix->[0];
	my $G = QIF::Matrix->random(2,2);


	$C = QIF::Matrix->new(q{
		0.0181 0.9819
		0.5489 0.4511
	});
	$G = QIF::Matrix->new(q{
		0.2122 0.6462
		0.1664 0.6792
	});
	$p = 0.455513691093722;
	$p = [$p, 1-$p];



	sub prior_gain {
		my ($w) = @_;
		return
			$p->[0] * $G->[$w][0] +
			$p->[1] * $G->[$w][1];
	}
	sub post_gain {
		my ($y, $w) = @_;
			$p->[0] * $G->[$w][0] * $C->[0][$y] +
			$p->[1] * $G->[$w][1] * $C->[1][$y];
	}

	my $wp = prior_gain(0) >= prior_gain(1) ? 0 : 1;
	my $w0 = post_gain(0, 0) >= post_gain(0, 1) ? 0 : 1;
	my $w1 = post_gain(1, 0) >= post_gain(1, 1) ? 0 : 1;

	warn rand(), "--$w0-$w1";
	next if $w0 == $w1;
	warn "-";

	my $deriv_prior =
		$G->[$wp][0] - $G->[$wp][1];

	my $deriv_post =
		$C->[0][0] * $G->[$w0][0] - $C->[1][0] * $G->[$w0][1] +
		$C->[0][1] * $G->[$w1][0] - $C->[1][1] * $G->[$w1][1];


	if(abs($deriv_post) - abs($deriv_prior) > 1e-6) {
		warn "C:\n$C\nG:\n$G\np: $p->[0]\n deriv prior: $deriv_prior\n deriv_post: $deriv_post\n";
		warn "wp: $wp\nw0: $w0\nw1: $w1\n";
		last;
	}
}



#my $a1 = $G->[0][0] - $G->[0][1];
#my $a2 = $G->[1][0] - $G->[1][1];
#warn $G;
#warn $a1;
#warn $a2;
#
#my $b1 = $C->[0][0]*$G->[0][0] - $C->[1][0]*$G->[0][1] + $C->[0][1]*$G->[1][0] - $C->[1][1]*$G->[1][1];
#my $b2 = $C->[0][0]*$G->[1][0] - $C->[1][0]*$G->[1][1] + $C->[0][1]*$G->[0][0] - $C->[1][1]*$G->[0][1];
#warn $b1;
#warn $b2;
