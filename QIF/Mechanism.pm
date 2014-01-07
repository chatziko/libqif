package QIF::Mechanism;
use Moose;
use feature 'state';

use Storable qw/freeze thaw/;
use IPC::ShareLite;

use QIF::Matrix qw/_equal _less_than_or_eq/;
use QIF::LinearProgram;

extends 'QIF::Channel';

has graph	=> (is => 'rw', isa => 'Maybe[QIF::Graph]', required => 0, );
has metric	=> (is => 'rw', isa => 'Maybe[QIF::Matrix]', required => 0, );

use constant PI	=> 4 * atan2(1, 1);
use constant E	=> exp(1);


sub optimal {
	my ($class, $graph, $epsilon, $p, $G) = @_;

	my $X = $graph->nodes;
	my $A = $graph->adjacency;
	$p ||= QIF::Matrix->uniform($X);
	$G ||= QIF::Matrix->identity($X);

	$p->rows == 1 && $p->cols == $X		or die 'invalid p';
	$G->rows == $X && $G->cols == $X	or die 'invalid G';

	# Variables are C_x,y
	# 3 things to build
	my $eq   = QIF::Matrix->new([]);
	my $ineq = QIF::Matrix->new([]);
	my $obj  = QIF::Matrix->new([]);

	# objective function (assume direct mechanism)
	#   V(C) = sum_x,y p_x C_x,y G_x,y
	#
	my @coeff = (0)x($X*$X);
	for my $x (0..$X-1) {
		for my $y (0..$X-1) {
			$coeff[ $x*$X+$y ] += $p->[0][$x] * $G->[$x][$y];
		}
	}
	push @coeff, 0;		# constant coeff 0
	push @$obj, \@coeff;

	# equality to express that the sum of values is 1
	#
	for my $x (0..$X-1) {
		my @coeff = (0)x($X*$X);
		for my $y (0..$X-1) {
			$coeff[ $x*$X+$y ] = 1;
		}
		push @coeff, -1;		# constant coeff 1
		push @$eq, \@coeff;
	}

	# inequalities for each edge
	#
	for my $x1 (0..$X-1) {
		for my $x2 (0..$X-1) {
			$A->[$x1][$x2]	or next;	# only connected nodes

			for my $y (0..$X-1) {
				# alpha C[x1,y] <= C[x2,y]
				my @coeff = (0)x($X*$X);
				$coeff[ $x1*$X+$y ] = exp(-$epsilon);
				$coeff[ $x2*$X+$y ] = -1;
				push @coeff, 0;		# constant coeff 0
				push @$ineq, \@coeff;
			}
		}
	}

	my $program = QIF::LinearProgram->new(objective => $obj, equalities => $eq, inequalities => $ineq);
	my ($z, $sol) = $program->solve;

	# build $C from the solution
	my $C = QIF::Matrix->empty;	
	for my $x (0..$X-1) {
		for my $y (0..$X-1) {
			$C->[$x][$y] = $sol->[0][$x*$X+$y];
		}
	}

	return $class->new(matrix => $C, graph => $graph);
}

sub tight_constraint {
	my ($class, $epsilon, %opt) = @_;

	my $graph = $opt{graph};
	my $metric = $opt{metric};
	$graph || $metric			or die "must give graph or metric (or both)";

	$metric ||= $graph->metric;
	my $nodes = $graph ? $graph->nodes : $metric->rows;

	# NOTE: using $uniform instead of $ones somehow improves stability, making $diag non-negative when it should be
	#
	my $phi = $metric->scalar_powers(exp(-$epsilon));
	my $phi_inv = $phi->inverse;
	#my $ones = QIF::Matrix->ones(1, $nodes);
	my $uniform = QIF::Matrix->uniform($nodes);

	my $diag = $uniform * $phi_inv;
	#$diag->cleanup_zeroes;
	$diag->is_non_negative		or die "negative diag";

	my @c;
	for my $i (0..$nodes-1) {
		for my $j (0..$nodes-1) {
			$c[$i][$j] = $nodes * $diag->[0][$j] * $phi->[$i][$j];
			#warn("$c[$i][$j] = $nodes * $diag->[0][$j] * $phi->[$i][$j]\n")	if $c[$i][$j] < 0;
		}
	}
	my $matrix = QIF::Matrix->new(\@c);
	#$matrix->is_stochastic	or die 'non-stochastic';	# sanity check

	# find optimality range
#	my ($high, $low);
#	my $step = 0.01;
#	for ($high = 0; $high <= 1; $high += $step) {
#		my $a = 1/$nodes * (1+$high);
#		my $p = QIF::Matrix->new([[$a, ((1-$a)/($nodes-1))x($nodes-1)]])->transpose;
#		my $Y = $Ginv * $p;
#		$Y->is_non_negative or last;
#	}
#	for ($low = 0; $low <= 1; $low += $step) {
#		my $a = 1/$nodes * (1-$low-$step);
#		my $p = QIF::Matrix->new([[$a, ((1-$a)/($nodes-1))x($nodes-1)]])->transpose;
#		my $Y = $Ginv * $p;
#		$Y->is_non_negative or last;
#	}

	return $class->new(matrix => $matrix, graph => $graph, metric => $metric);#, [$low, $high];
}

sub tight_constraint_grid {
	my ($class, $size, $step, $epsilon) = @_;

	my $s = `./tight_grid/tight_grid $size $size $step $epsilon`;
	if($? != 0) {
		warn "out of memory"	if $! =~ /allocate memory/;
		die "cannot create tight_constraints $s $!";
	}

	$s =~ s/.*\n.*\n//;		# remove 2 lines
	my $matrix = QIF::Matrix->new($s);

	my $ls = QIF::LocationSet->grid($size, $size, $step);
	return $class->new(matrix => $matrix, metric => $ls->metric);
}

sub geometric {
	my ($class, $n, $epsilon) = @_;

	my $g = QIF::Graph->line($n);

	# compute metricance matrix fast
	my @metric;
	for my $i (0..$n-1) {
		$metric[$i][$i] = 0;

		for my $j ($i+1..$n-1) {
			$metric[$i][$j] =
			$metric[$j][$i] = abs($i - $j);
		}
	}
	my $metric = QIF::Matrix->new(\@metric);

	return $class->tight_constraint($epsilon, graph => $g, metric => $metric);
}

sub comp_parallel {
	my ($self, $mech) = @_;
	my $M1 = $self->matrix;
	my $M2 = $mech->matrix;
	
	my $n1 = $M1->rows;
	my $n2 = $M2->rows;
	my $m1 = $M1->cols;
	my $m2 = $M2->cols;

	my $new = QIF::Matrix->new($n1*$n2, $m1*$m2);
	for my $i1 (0..$n1-1) {
	for my $i2 (0..$n2-1) {
		for my $j1 (0..$m1-1) {
		for my $j2 (0..$m2-1) {
				$new->[$i1*$n2 + $i2][$j1*$m2 + $j2] = $M1->[$i1][$j1] * $M2->[$i2][$j2];
		}}
	}}
	return QIF::Mechanism->new(matrix => $new);	
}

sub is_private {
	my ($self, $eps) = @_;
	my $e_eps = exp($eps);

	my $adj = $self->graph->adjacency	if $self->graph;
	my $d = $self->metric;

	my $M = $self->matrix;
	my $n = $M->rows;
	my $m = $M->cols;

	for my $i1 (0..$n-1) {
		for my $i2 ($i1+1..$n-1) {
			my $a;
			if($adj) {
				$adj->[$i1][$i2]	or next;
				$a = $e_eps;
			} else {
				$a = exp($eps * $d->[$i1][$i2]);
			}

			for my $j (0..$m-1) {
				$M->[$i1][$j] <= $M->[$i2][$j]
					? _less_than_or_eq($M->[$i2][$j], $a * $M->[$i1][$j])
					: _less_than_or_eq($M->[$i1][$j], $a * $M->[$i2][$j])
					or return 0;
			}
		}
	}
	return 1;
}

sub smallest_epsilon {
	my ($self, $only_first_row) = @_;
	my $max = 0;

	my $adj = $self->graph->adjacency	if $self->graph;
	my $d = $self->metric;

	my $M = $self->matrix;
	my $n = $M->rows;
	my $m = $M->cols;

	for my $i1 (0..$n-1) {
		for my $i2 ($i1+1..$n-1) {
			my $dist;
			if($adj) {
				$adj->[$i1][$i2]	or next;
				$dist = 1;
			} else {
				$dist = $d->[$i1][$i2];
			}

			for my $j (0..$m-1) {
				my ($small, $large) = sort { $a <=> $b } $M->[$i1][$j], $M->[$i2][$j];

				next if $small < 0.0001;
				if(_equal($small, 0)) {
					_equal($large, 0)	or return undef;		# zero vs non-zero, no eps exists
					next;													# both zeros, no test needed
				}

				my $eps = log($large/$small)/$dist;
				$max = $eps	if $eps > $max;			# the smallest eps satisfying privacy is the _max of the $eps
			}
		}
		last if $only_first_row;
	}
	return $max;
}

sub utility {
	my ($self, $p, $G) = @_;
	$p ||= QIF::Matrix->uniform($self->input_no);
	#$G ||= QIF::Matrix->identity($self->input_no);
	return $G ? $self->g_vulnerability($p, $G) : $self->vulnerability($p);
}

sub plannar_laplace {
	my ($class, $loc_set, $epsilon) = @_;
	my $draws = 10000;
	my $rows_per_run = 50;

	my $points = $loc_set->points;
	my $n = @$points;
	my @m;

	while (@m < $n) {
		# Algorithm::KNN::XS has a memory leak. As a workaroud
		# we do a batch of $rows_per_run rows in a child, so that memory is freed
		#
		my $new_rows = exec_in_child(sub {
			my @r;
			for my $i (0..$rows_per_run-1) {
				#warn @m+$i;
				my $real = $points->[@m + $i]	or last;
				$r[$i] = [(0) x $n];

				for (1..$draws) {
					my $z = addPolarNoise($epsilon, $real);
					my $closest = $loc_set->closest($z);
					$r[$i][$closest]++;
				}
				$_ /= $draws	for @{ $r[$i] };
			}
			return \@r;
		});
		push @m, @$new_rows;
	}

	my $matrix = QIF::Matrix->new(\@m);
	$matrix->is_stochastic	or die "not stochastic";

	return $class->new(matrix => $matrix, metric => $loc_set->metric);
}

sub plannar_laplace_grid {
	my ($class, $size, $step, $epsilon) = @_;

	my $s = `./plannar_laplace/laplace $size $size $step $epsilon`;
	$s =~ s/.*\n.*\n//;		# remove 2 lines
	my $matrix = QIF::Matrix->new($s);

	my $ls = QIF::LocationSet->grid($size, $size, $step);
	return $class->new(matrix => $matrix, metric => $ls->metric);
}

sub lambertW {
	my ($x) = @_;

	#min_diff decides when the while loop should stop
	my $min_diff = 1e-10;
	if (_equal($x, -1/E)) {
		return -1;

	} elsif ($x < 0 && $x > -1/E) {
		my $q = log(-$x);
		my $p = 1;
		while (abs($p-$q) > $min_diff) {
			$p = ($q * $q + $x / exp($q)) / ($q+1);
			$q = ($p * $p + $x / exp($p)) / ($p+1);
		}
		#This line decides the precision of the float number that would be returned
		return int(1000000 * $q - .5) / 1000000;

	} elsif (_equal($x, 0)) {
		return 0;

	} else {
		#TODO why do you need this if branch? 
		return 0;
	}
}

sub inverseCumulativeGamma {
	my ($epsilon, $z) = @_;
	my $x = ($z-1) / E;
	return - (lambertW($x) + 1) / $epsilon;
}

sub addPolarNoise {
	my ($epsilon, $pos) = @_;

	# random number in [0, 2*PI)
	my $theta = rand() * PI * 2;
	# random variable in [0,1)
	my $z = rand();
	my $r = inverseCumulativeGamma($epsilon, $z);

	return [
		$pos->[0] + $r * cos($theta),
		$pos->[1] + $r * sin($theta)
	];
}

# runs code after forking, waits and returns result
sub exec_in_child {
	my $code = shift;

	state $share = IPC::ShareLite->new(
		-key     => 1971,
		-create  => 'yes',
		-destroy => 'no'
	) or die $!;

	if(my $pid = fork) {
		waitpid $pid, 0;
		return ${ thaw $share->fetch };

	} else {
		my $res = $code->();
		$share->store(freeze \$res);
		exit;
	}
}

1;
