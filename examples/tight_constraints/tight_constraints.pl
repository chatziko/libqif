#!/usr/bin/perl -w
use strict;

use POSIX qw/ceil/;
use List::Util qw/max/;
use Parallel::Iterator qw/iterate_as_array/;

use QIF::Channel;
use QIF::Graph;
use QIF::Mechanism;
use QIF::LocationSet;

$QIF::Matrix::USE_RAT = 0;
$QIF::Matrix::CONVERT_FRACTIONS = 0;

my $u = 150;
my $v = 5;

my $g;
my $epsilon;
#$g = QIF::Graph->hamming($u, $v);
#$g = QIF::Graph->line(2);
#$g = QIF::Graph->cycle($n);
#$g = QIF::Graph->cycle($n);
#$g = QIF::Graph->complete(4);
#$g = QIF::Graph->random($n, 0.7);
#my $g2 = QIF::Graph->line($n);

#$g = $g->strong_product(QIF::Graph->line(3));
#die $g->to_png('test.png');
#die $g;




# build prior
#
#my @value_prob = (.3, .27, .23, .2);
#my $prior = [[]];
#for (@{$g->labels}) {
#	my $p = 1;
#	$p *= $value_prob[$_]	for split ',';
#	push @{$prior->[0]}, $p;
#}
#$prior = QIF::Matrix->new($prior);
#die 'bad prior' unless $prior->is_stochastic;


my @c = `cat /proc/cpuinfo | grep processor`;
my $cpus = @c;


sub experiment_geo {
	my $size = 5;
	my $step = 1;
	my $ls = QIF::LocationSet->grid($size, $size, $step);
	my $prior = QIF::Matrix->uniform(scalar @{$ls->points});
	my $gain = $ls->gain_function;

	print "size: $size x $size, step: $step\n\n";

	my @epsilons;
	for($epsilon = 0.5; $epsilon <= 1.3; $epsilon += 0.2) {
		push @epsilons, $epsilon;
	}

	$cpus = 1;
	iterate_as_array
		{ workers => $cpus },
		sub {
			my ($id, $epsilon) = @_;

			#my $mech = eval { QIF::Mechanism->tight_constraint_grid($size, $step, $epsilon) }
			my $mech = eval { QIF::Mechanism->tight_constraint($epsilon, metric => $ls->metric) }
				or print("$epsilon, not in range\n"), return;
			#my $pl = QIF::Mechanism->plannar_laplace_grid($size, $step, $epsilon);
			my $pl = QIF::Mechanism->plannar_laplace($ls, $epsilon);

			print join ', ',
				$epsilon,
				$mech->utility($prior),
				$pl->utility($prior);
				#$mech->utility($prior, $gain),
				#$pl->utility($prior, $gain);
				#$mech->smallest_epsilon // 'infty',
				#$pl->smallest_epsilon // 'infty';
			print "\n";
		},
		\@epsilons;
}


sub experiment_sum {
	$u = 150;
	$v = 5;

	$g = QIF::Graph->radius($u*$v+1, $v);

	my $prior = QIF::Matrix->uniform($g->nodes);

	warn "building distance matrix";
	my $dist = radius_metric($g, $v);
	#my $dist = $g->metric;
	warn "got distance matrix";


	for($epsilon = 0.79; $epsilon <= 1.31; $epsilon += .01) {
		print "$epsilon, ";
		#print new_bound($prior, $dist, $epsilon), ", ";

		my $mech = eval { QIF::Mechanism->tight_constraint($epsilon, graph => $g, metric => $dist) }
			or print("not in range\n"), next;
		print $mech->utility($prior), ", ";

		my $geo = QIF::Mechanism->geometric($g->nodes, $epsilon / $v);
		print $geo->utility($prior), "\n";
	}
}

sub experiment_2_count {
	$u = 30;

	$g = QIF::Graph->line($u+1);
	$g = $g->strong_product($g);

	#die $g->to_png('test.png');

	my $prior = QIF::Matrix->uniform($g->nodes);

	warn "building distance matrix";
	my $dist = n_count_metric($g);
	#my $dist = $g->metric;
	warn "got distance matrix";

	for($epsilon = 0.89; $epsilon <= 1.301; $epsilon += .01) {
		print "$epsilon, ";
#		my $nb = eval { new_bound($prior, $dist, $epsilon) }
#			or print("not in range\n"), next;
#		print $nb, ", ";

		my $mech = eval { QIF::Mechanism->tight_constraint($epsilon, graph => $g, metric => $dist) }
			or print("not in range\n"), next;
		print $mech->utility($prior), ", ";
		#print "---\n$mech";

		my $geo = QIF::Mechanism->geometric($u+1, $epsilon / 2);
		my $geo2 = $geo->comp_parallel($geo);
		$geo2->matrix->rows == $g->nodes	or die 'bad size';
		#print "\n\n$geo2";

		$geo2->graph($g);
		$geo2->is_diff_private($epsilon) or die "not dp";

		print $geo2->utility($prior), "\n";
	}
}

#experiment_2_count();
#experiment_sum();
experiment_geo();


exit;


#	print $epsilon . ": " .
#		"old bound: " . old_bound($u, $v, $epsilon) . ", " .
#		"new bound: " . new_bound($prior, $dist, $epsilon) . "\n";
#}

exit;

sub new_bound {
	my ($prior, $dist, $epsilon) = @_;

	my $phi = $dist->scalar_powers(exp(-$epsilon));
	my $phi_inv = $phi->inverse;

	my $y = $prior * $phi_inv;
	die "not in range"		unless $y->is_non_negative;

	#return QIF::Channel::log2( $y->sum / $prior->max );
	return $y->sum;
}


sub old_bound {
	my ($u, $v, $epsilon) = @_;
	my $exp_eps = exp($epsilon);
	return $u * QIF::Channel::log2(($v * $exp_eps) / ($v -1 + $exp_eps));
}

# distance matrix for hamming graphs, much faster than $g->metric
sub hamming_metric {
	my ($g) = @_;
	my $n = $g->nodes;
	my $lab = $g->labels;

	my @dist;
	for my $i (0..$n-1) {
		$dist[$i][$i] = 0;

		my @val1 = split ',', $lab->[$i];
		for my $j ($i+1..$n-1) {
			my @val2 = split ',', $lab->[$j];
			$dist[$i][$j] =
			$dist[$j][$i] =
				grep { $val1[$_] != $val2[$_] } 0..@val1-1;
		}
	}
	return QIF::Matrix->new(\@dist);
}

sub n_count_metric {
	my ($g) = @_;
	my $n = $g->nodes;
	my $lab = $g->labels;

	my @dist;
	for my $i (0..$n-1) {
		$dist[$i][$i] = 0;

		my @val1 = split ',', $lab->[$i];
		for my $j ($i+1..$n-1) {
			my @val2 = split ',', $lab->[$j];
			$dist[$i][$j] =
			$dist[$j][$i] =
				max map { abs($val1[$_] - $val2[$_]) } 0..@val1-1;
		}
	}
	return QIF::Matrix->new(\@dist);
}

sub radius_metric {
	my ($g, $v) = @_;
	my $n = $g->nodes;

	my @dist;
	for my $i (0..$n-1) {
		$dist[$i][$i] = 0;

		for my $j ($i+1..$n-1) {
			$dist[$i][$j] =
			$dist[$j][$i] =
				ceil(abs($i - $j) / $v);
		}
	}
	return QIF::Matrix->new(\@dist);
}



my $nodes = $g->nodes;
my ($tight, $z) = QIF::Mechanism->tight_constraint($epsilon, graph => $g);

my $a = 0.25;
my $b = (1-$a)/($nodes-1);
my $p = QIF::Matrix->new([[$a, ($b)x($nodes-1)]]);

#$p = QIF::Matrix->new([[(2/$nodes, 0)x($nodes/2)]]);

#my $p = QIF::Matrix->uniform(5);
#$p = QIF::Matrix->new('0.2 0 0.2 0 0.2 0 0.2 0 0.2 0');

print "p: $p\n";

my $optimal = QIF::Mechanism->optimal($g, $epsilon, $p);
print "tight:\n$tight\n";
print "optimal:\n$optimal\n";
warn $tight->matrix->is_stochastic;

my $tu = $tight->utility($p);
my $ou = $optimal->utility($p);
warn "tight: $tu\noptimal: $ou\ndiff: " . ($ou - $tu) . "\n"; 
warn "ratio: " . $tu/$ou . "\n";
exit;

for(my $alpha = 0.1; $alpha < 0.99; $alpha += 0.1) {
	eval {
		my ($m, $z) = QIF::Mechanism->tight_constraint($alpha, graph => $g);
		printf "%.1f: %.2f - %.2f\n", $alpha , @$z;
	};
}
