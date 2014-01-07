#!/usr/bin/perl -w
use strict;


use QIF::Channel;
use QIF::Compare;
use QIF::Compare::FindGain;

$QIF::Matrix::USE_RAT = 1;

my $Ex1 = QIF::Channel->new(matrix => q{
	1/4  3/4
	1/4  3/4
	6/10 4/10
});
my $Ex2 = QIF::Channel->new(matrix => q{
	1/2   0 1/2
	  0 1/2 1/2
	1/2 1/2   0
});
my $Ex3 = QIF::Channel->new(matrix => q{
	1/3 1/3 1/3   0
	  0 1/3 1/3 1/3
	1/3 2/3   0   0
});
my $Ex4 = QIF::Channel->new(matrix => q{
	1/2   0   0 1/2
	  0 1/2 1/2   0
	1/2 1/2   0   0
});
my $Ex5 = QIF::Channel->new(matrix => q{
	2/3  1/3
    2/3   1/3
    1/4   3/4
});
my $Ex6 = QIF::Channel->new(matrix => q{
	1/2 1/2   0
	1/2   0 1/2
	  0 1/2 1/2
});

my $Cat1 = QIF::Channel->new(matrix => q{
	1/2 1/2 0
	0	1/2 1/2
	1/2	0	1/2
});
my $Cat2 = QIF::Channel->new(matrix => q{
	1/2	0	0	1/2
	0	1/2	0	1/2
	0	0	1/2	1/2
});

my $Gid_compl = QIF::Matrix->new(q{
	1 1 0
	1 0 1
	0 1 1
});

my $Gnonzero = QIF::Matrix->new(q{
	153/296 0       21/148
	0       289/296 1     
	1/2     63/296  0
});

my $Test1 = QIF::Channel->new(matrix => q{
	.2  .22 .58
	.2  .4  .4 
	.35 .4  .25
});
my $Test2 = QIF::Channel->new(matrix => q{
	.1 .4 .1 .4 
	.2 .2 .3 .3
	.5 .1 .1 .3 
});
#die unless $Test1->matrix->is_stochastic;
#die unless $Test2->matrix->is_stochastic;
#die "Factorizable:\n".$Test1->factorize($Test2) if $Test1->factorize($Test2);

die "Factorizable:\n".$Ex5->factorize($Ex6) if $Ex5->factorize($Ex6);

die unless $Ex5->matrix->is_stochastic;
die unless $Ex6->matrix->is_stochastic;

my $comp;
$comp = QIF::Compare->new(c1 => $Ex5, c2 => $Ex6);
$comp->add_01_functions;
#$comp->add_2_block;
#$comp->add_function($Gid_compl);
#$comp->add_function($Grand);
#$comp->add_identity;

#print "$_\n"	for @{ $comp->Gs };
print $comp->compare;
exit;
#die scalar @{ $comp->Gs };


my $unif = QIF::Matrix->uniform(3);
warn $Test1->g_vulnerability($unif, $Gnonzero);
warn $Test2->g_vulnerability($unif, $Gnonzero);
warn $Test1->g_vulnerability($unif, $Gnonzero)- $Test2->g_vulnerability($unif, $Gnonzero);
exit;

$comp = QIF::Compare::FindGain->new(c1 => $Test1, c2 => $Test2);
print $comp->compare(3);
exit;

while(1) {
	#$Test2 = QIF::Channel->random(3, 4);
	#$Test2 = $Cat2;

	#my $N = $Test2->matrix;
	#my $NR = $N->right_inverse;
	#!($NR * $N)->has_non_negative_row		or next;		# be sure there is no Cat-vector

	print "Random 3x4 matrix\n$Test2";

	my $c = 1;
	while($c <= 10) {
		#$Test1 = QIF::Channel->random(3,3);
		#$Test1 = $Cat1;

		# should not be factorizable
		next if $Test1->factorize($Test2);

		print "checking $c\n";
		$c++;

		#my $G = create_g($Test1, $Test2);

		$comp->c1($Test1);
		$comp->c2($Test2);

		my ($x, $G) = $comp->compare;
		next if $x;		# C1 leaks more for this $x, as expected

		die "This C1:\n$Test1\nleaks less for all " . scalar(@{$comp->Gs}) . " gain function(s)\n";
	}
}
exit;


print "Comparing:\nC1:\n", $comp->c1->matrix, "\nC2:\n", $comp->c2->matrix, "\n",
	"using " . @{$comp->Gs} . " gain function(s)\n\n";

my ($x, $G) = $comp->compare;
if($x) {
	print "C1 leaks more than C2 for the following prior:\n",
		$x, "and gain function:\n$G";
} else {
	print "C1 always leaks less than C2\n";
}





## failed methods to create a 01 G for a channel
##
#sub create_g {
#	my ($C1, $C2) = @_;
#	my $A = $C1->matrix;
#	my $B = $C2->matrix;
#	my $BR = $B->right_inverse;
#
#	print "create_g\nBR:\n".$BR->sign . "\n";
#	print "BR * A:\n" . ($BR*$A)->sign . "\n";
#	print "BR * B:\n" . ($BR*$B)->sign . "\n";
#
#	my $X = $B->rows;
#	my $vecs = [];
#
#	for my $i (1..$BR->rows) {
#		my $row = $BR->row($i);
#		!($row * $A)->is_non_negative	or next;
#
#		my @pos = grep { $row->[0][$_] > -1e-6 } 0..$X-1;
#		#@pos < @$row	or next; # need both positive and negatice
#
#		my $vec1 = [(0)x$X];
#		$vec1->[$_] = 1	for @pos;
#
#		my $vec2 = [(1)x$X];
#		$vec2->[$_] = 0	for @pos;
#
#		push @$vecs, $vec1, $vec2;
#	}
#
#	# remove subsets (or duplicates)
#	my $vecs2 = [];
#	while(my $vec = shift @$vecs) {
#		push @$vecs2, $vec
#			unless grep { subset($vec, $_) } @$vecs, @$vecs2;
#	}
#
#	return QIF::Matrix->new($vecs2)->transpose;
#}
#
## $a, $b are 0-1 vectors
## true if $a \subseteq $b
#sub subset {
#	my ($a, $b) = @_;
#	for(0..@$a-1) {
#		return 0 if $a->[$_] && !$b->[$_];
#	}
#	return 1;
#}
#
#
