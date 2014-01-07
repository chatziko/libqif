#!/usr/bin/perl -w
use strict;


use QIF::Channel;
use QIF::Compare;
use Data::Dumper;

$QIF::Matrix::USE_RAT = 0;


my $C1 = QIF::Channel->new(matrix => q{
	0.401512793319434 0.386437854618111 0.212049352062455 
	0.285301468829656 0.655709524810182 0.0589890063601622
	0.321834296212112 0.357309302613892 0.320856401173996
});
my $C2 = QIF::Channel->new(matrix => q{
	0.207644637734501 0.0323443768197586 0.489205183433045 0.270805802012696
	0.718485185673515 0.0326679945862032 0.074702089823173 0.174144729917109
	0.179920648024843 0.0693963166803889 0.249315288719217 0.50136774657555
});

print "factorizes? ".$C1->factorize($C2),"\n\n";

my $A = $C1->matrix;
my $B = $C2->matrix;
my $BR = $B->right_inverse;

print "BR * B:\n" . ($BR*$B)->sign . "\n";
print "BR * A:\n" . ($BR*$A)->sign . "\n";

print "BR\n" . $BR->sign;


my $comp = QIF::Compare->new(c1 => $C1, c2 => $C2);

my $Gw = QIF::Matrix->new(q{
	1 0
	1 0
	0 1
});

$comp->add_function($Gw);


print "Comparing:\nC1:\n", $comp->c1->matrix, "\nC2:\n", $comp->c2->matrix, "\n",
	"using " . @{$comp->Gs} . " gain function(s)\n\n";

my ($x, $G) = $comp->compare;
if($x) {
	print "C1 leaks more than C2 for the following prior:\n",
		$x, "and gain function:\n$G";
} else {
	print "C1 always leaks less than C2\n";
}


sub print_sign {
	my ($m) = @_;

	for my $row (@$m) {
		for (@$row) {
			print $_ < 0 ? '-' : '+', ' ';
		}
		print "\n";
	}
}



