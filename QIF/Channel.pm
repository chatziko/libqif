package QIF::Channel;
use Moose;
use Moose::Util::TypeConstraints;

use List::Util qw/sum max min/;

use QIF::Matrix;

use overload fallback => 1,
	'""' => \&stringify,
	'*' => \&multiply,
	'==' => \&is_equal_to;

has matrix	=> (is => 'rw', isa => 'QIF::Matrix', required => 1, coerce => 1);

# Channel.pm
#
# Class for channel-related operations
#

our $log2 = log(2);

sub new_symmetric {
	my ($class, $size, $bigval) = @_;
	my $lowval = (1 - $bigval) / ($size-1);

	my $c = [];
	for my $i (0..$size-1) {
		push @$c, [];
		for my $j (0..$size-1) {
			$c->[$i][$j] = $i == $j ? $bigval : $lowval;
		}
	}

	my $self = { matrix => $c };
	return bless $self;
}

sub random {
	my ($class, $m, $n) = @_;
	return $class->new(
		matrix => QIF::Matrix->random($m, $n)->map_rows(sub {
			my ($row) = @_;
			my $s = sum(@$row);
			return [ map { $_/$s } @$row ];
		})
	);
}

# input_no ()
#   Returns the number of input values of the channel (the number of rows of the matrix)
#
sub input_no {
	my ($self) = @_;
	return $self->matrix->rows;
}

# output_no ()
#   Returns the number of output values of the channel (the number of columns of the matrix)
#
sub output_no {
	my ($self) = @_;
	return $self->matrix->cols;
}

# capacity ($e)
#   Computes the capacity of the channel using the Arimoto-Blahut algorithm
#   $e is the accuracy to which the capacity should be computed
#
sub capacity {
	my ($self, $e) = @_;
	$e ||= 0.001;		# default value
	my $matrix = $self->matrix;

	my $m = scalar @{ $matrix };		# no of X's values
	my $n = scalar @{ $matrix->[0] };	# no of Y's values

	my @F = map 0, 1..$m;
	my @Px = map 1/$m, 1..$m;
	my @Py = map 0, 1..$n;

	while(1) {
		# update Py
		for my $k (0..$n-1) {
			$Py[$k] = 0;
			$Py[$k] += $matrix->[$_][$k] * $Px[$_]
				for 0..$m-1;
		}

		# update F
		for my $j (0..$m-1) {
			my $s = 0;
			# NOTE: this is e-base log, not 2-base!
			$s += $matrix->[$j][$_] > 0
				? $matrix->[$j][$_] * log( $matrix->[$j][$_] / $Py[$_])
				: 0
				for 0..$n-1;
			$F[$j] = exp($s);
		}

		# check stop condition
		my $x = 0;
		$x += $F[$_] * $Px[$_]
			for 0..$m-1;
		my $IL = log2($x);
		my $IU = log2(max(@F));

		return ($IL, \@Px)		if $IU - $IL < $e;

		# update Px
		$Px[$_] *= $F[$_] / $x
			for 0..$m-1;
	}
}

# capacity_sym ()
#   Computes the capacity of the channel assuming it is symmetric (in this case C = log|O| - H(row)
#   NOTE: symmetry is assumed, we don't test it
#
sub capacity_sym {
	my ($self) = (@_);
	my $row = $self->matrix->[0];

	my $s = 0;
	$s += $_ == 0 ? 0 : - $_ * log2($_)
		for @$row;

	return log2(scalar @$row) - $s;
}

# capacity_partsym ($indexes)
#   Computes the capacity of the channel assuming that it is partially symmetric (some columns
#   are constant, the rest are symmetric). It uses the formula from the channels paper.
#   $indexes is an array with the number of the columns
#   that form the symmetric part. No testing of symmetry is performed
#
sub capacity_partsym {
	my ($self, $indexes) = @_;
	my $row = $self->matrix->[0];

	# calculate the sum and the entropy of the symmetric part of a row
	my ($entr, $sum) = (0, 0);
	for(@$indexes) {
		my $val = $row->[$_];
		$sum += $val;
		$entr += $val == 0 ? 0 : - $val * log2($val);
	}

	return $sum * log2(@$indexes / $sum) - $entr;
}

sub probable_innocence {
	my ($self) = @_;
	my $matrix = $self->matrix;

	my $n = scalar @{ $matrix };		# no of rows
	my $m = scalar @{ $matrix->[0] };	# no of cols

	for my $z (0..$m-1) {
		for my $i (0..$n-1) {
			my $test = ($n-1) * $matrix->[$i][$z];
			for my $j (0..$n-1) {
				# check that (n-1)p(o_z | a_i) >= p(o_z | a_j)
				return undef, $z, $i, $j
					if $test < $matrix->[$j][$z];
			}
		}
	}
	return 1;
}

sub vulnerability {
	my ($self, $x) = @_;
	$x = $x->[0]	if ref $x->[0] eq 'ARRAY';	# support both vectors and matrices

	my $m = $self->{matrix};
	my $s = 0;
	for my $i(0.. @{$m->[0]}-1) {
		my $max = 0;
		for my $j(0..@$m-1) {
			$max = $m->[$j][$i] * $x->[$j]
				if $max < $m->[$j][$i] * $x->[$j];
		}
		$s += $max;
	}
	return $s;
}

sub g_vulnerability {
	my ($self, $p, $G) = @_;
	$p = $p->[0]	if ref $p->[0] eq 'ARRAY';	# support both vectors and matrices

	my $M = $self->matrix;

	my $X = @$M;
	$X == @$p		or die "invalid p";
	my $Y = @{$M->[0]};
	my $K = @{$G->[0]};

	# Compute sum_y max_k sum_x p_x M_x,y G_x,k
	#
	my $res = 0;
	for my $y (0..$Y-1) {
		my $max = 0;
		for my $k(0..$K-1) {
			my $g = 0;
			for my $x (0..$X-1) {
				$g += $p->[$x] * $M->[$x][$y] * $G->[$x][$k];
			}
			$max = $g	if $max < $g;
		}
		$res += $max;
	}
	return $res;
}

# prob_error ($x)
#
# Returns the bayesian probability of error, using $x as the input distribution
#
sub prob_error {
	my ($self, $x) = @_;
	return 1 - $self->vulnerability($x);
}

# 1 - 2^-H(A|O)
#
sub santhi_vardy {
	my ($self, $x) = @_;
	return 1 - 2**(- $self->cond_entropy($x));
}

# 1/2 H(A|O)
#
sub hellman_raviv {
	my ($self, $x) = @_;
	return $self->cond_entropy($x) / 2;
}

# computes H(A|O) for given input distribution
#
sub cond_entropy {
	my ($self, $x) = @_;

	my $c = $self->matrix;
	my $n = @$c;
	my $m = @{$c->[0]};

	my $s = 0;
	for my $j (0..$m-1) {
		#my $xo = post($c, $x, $j);
		my $xo = [ map $x->[$_] * $c->[$_][$j], 0..@$c-1 ];		# $x .* the j-th column
		my $sum = 0; $sum += $_ for @$xo;
		if($sum > 0) {
			$_ /= $sum for @$xo;
		}

		for my $i (0..$n-1) {
			next unless $xo->[$i];
			$s -= $x->[$i] * $c->[$i][$j] * log2($xo->[$i]);
		}
	}

	return $s;
}

sub mutual_information {
	my ($self, $x) = @_;
	return entropy($x) - $self->cond_entropy($x);
}

sub has_zero {
	my ($self) = @_;
	for my $row(@{$self->matrix}) {
		return 1 if grep $_ == 0, @$row;
	}
	return 0;
}

# Returns a channel X such that $self = $C * $X
#
sub factorize {
	my ($self, $C) = @_;

	# A: our matrix, M x N
	# B: C's matrix, M x R
	# X: unknowns    R x N
	#
	my $A = $self->matrix;
	my $B = $C->matrix;
	
	my $M = $A->rows;
	my $N = $A->cols;
	my $R = $B->cols;
	$B->rows == $M	or return undef;
	my $vars = $R * $N;

	# Build equations for A = B X
	# We have R x N variables, that will be unfolded in a vector.
	# The varialbe X[r,n] will have variable number rN+n.
	# For each element m,n of A we have an equation, specifying that the inner product of
	# the m=th row of B and n-th column of X is A[m,n]
	#
	my $eq = [];
	for my $m (0..$M-1) {
		for my $n (0..$N-1) {
			my $row = [(0) x $vars, -$A->[$m][$n]];		# start with 0 coeff. The result (negated last element) is A[m,n]
			for my $r (0..$R-1) {
				$row->[ $r*$N+$n ] = $B->[$m][$r];		# coeff B[m,r] for variable X[r,n]
			}
			push @$eq, $row;
		}
	}

	# equalities for summing up to 1
	#
	for my $r (0..$R-1) {
		my $row = [(0) x $vars, -1];		# start with 0 coeff. The result (negated last element) is 1
		for my $n (0..$N-1) {
			$row->[ $r*$N+$n ] = 1;	# coeff 1 for variable X[r,n]
		}
		push @$eq, $row;
	}
	$eq = QIF::Matrix->new($eq);	# upgrade

	# we don't really care to maximize anything, so objective function is f(x) = 0
	#
	my $obj = QIF::Matrix->zeroes(1, $vars+1);

	# no inequalities
	#
	my $ineq = QIF::Matrix->empty;

	# solve program
	#
	my $program = QIF::LinearProgram->new(equalities => $eq, inequalities => $ineq, objective => $obj);
	my ($z, $sol) = $program->solve;
	$sol	or return undef;		# no solution

	# reconstrict channel from solution
	#
	my $X = [];
	for my $r (0..$R-1) {
		for my $n (0..$N-1) {
			$X->[$r][$n] = $sol->[0][ $r*$N+$n ];
		}
	}
	$X = QIF::Channel->new(matrix => QIF::Matrix->new($X));

	# verify solution
	$self == $C * $X		or die "factorize solution is invalid";

	return $X;
}

sub bounded_entropy_distance {
	my ($class, $x, $y) = @_;
	my $n = $x->cols;
	$n == $y->cols	or die 'size mismatch';

	my $max = 0;
	for(0..$n-1) {
		my $m = max($x->[0][$_], $y->[0][$_]);
		next if $m < 1e-7;

		my $l = abs($x->[0][$_] - $y->[0][$_]) / $m;
		$max = $l	if $max < $l;
	}
	return $max;
}

# to_string
#    Returns a string representation of the channel (useful for printing)
#
sub to_string {
	my ($self) = @_;

	my $s = "";
	for my $row (@{$self->matrix}) {
		$s .= join('  ', map {defined $_ ? sprintf("%.5f", $_) : "undef  "} @$row) . "\n";
	}

	return $s;
}


# Auxiliary functions
#
sub log2 {
	return log($_[0]) / $log2;
}

sub entropy {
	my ($row) = (@_);

	my $s = 0;
	$s += $_ == 0 ? 0 : - $_ * log2($_)
		for @$row;

	return $s;
}

sub stringify {
	my ($self) = @_;
	return "".$self->matrix;
}

sub multiply {
	my ($self, $C) = @_;
	return QIF::Channel->new(
		matrix => $self->matrix * $C->matrix
	);
}

sub is_equal_to {
	my ($self, $C) = @_;
	return $self->matrix == $C->matrix;
}

sub for_all_priors {
	my ($n, $step, $sub, $so_far, $k, $remaining) = @_;

	return for_all_priors($n, $step, $sub, QIF::Matrix->new([[]]), 0, 1)
		unless $so_far;

	if($k == $n-1) {
		$so_far->[0][$k] = $remaining;
		$sub->($so_far);

	} else {
		for(my $p = 0; $p <= $remaining; $p += $step) {
			$so_far->[0][$k] = $p;
			for_all_priors($n, $step, $sub, $so_far, $k+1, $remaining-$p);
		}
	}
}

sub max_over_priors {
	my ($n, $step, $sub) = @_;
	my $best_pi;
	my $max = -1e10;

	for_all_priors($n, $step, sub {
		my ($pi) = @_;
		my $z = $sub->($pi);
		return unless $z > $max;

		$max = $z;
		$best_pi = $pi->clone;
	});

	return wantarray ? ($max, $best_pi) : $max;
}

# max_prob_error
#
# Returns the maximum probability of error over all priors, using the subgradient descent method.
# Kostas: if I remember correctly this is based on the Java code to be included in Prism.
# Not sure whether I have tested the perl code.

sub max_prob_error {
	my ($self, $e) = @_;
	$e ||= 1E-6;		# default value

	my $c = $self->matrix;
	my $n = @$c;

	my $alpha = 1;

	my $x = [(1/$n)x$n];
	my $max = $self->prob_error($x);
	my $max_x = $x;

	my $k;
	LOOP: for($k = 1; $k < 10000; $k++) {
		my $new_x = [@$x];
		my $grad = _prob_error_subgradient($c, $x);
		$new_x->[$_] -= $alpha/$k * $grad->[$_]
			for 0..$n-1;

		# projection
		_project_to_simplex($new_x);

		my $new_val = $self->prob_error($new_x);
		if($new_val > $max) {
			last LOOP if $new_val - $max < $e;
			$max = $new_val;
			$max_x = $new_x;
		}

		$x = $new_x;
	}

	return $max, $max_x, $k;
}

sub _project_to_simplex {
	my ($x) = @_;
	my %done;
	my $n = @$x;

	while(1) {
		my $t = (sum(@$x)-1)/$n;
		my $n1 = 0;

		for my $i (0..@$x-1) {
			next if $done{$i};

			$x->[$i] -= $t;
			if($x->[$i] < 0) {
				$x->[$i] = 0;
				$done{$i} = 1;
				$n1++;
			}
		}

		last if $n1 == 0;		# no negative elements
		$n -= $n1;
	}
}

sub _prob_error_subgradient {
	my ($c, $x) = @_;

	my $n = @$c;
	my $m = @{$c->[0]};

	my $subgrad = [(0)x$n];

	for my $j (0..$m-1) {
		my $max = $c->[0][$j] * $x->[0];
		my $max_i = 0;
		for my $i (1..$n-1) {
			my $temp = $c->[$i][$j] * $x->[$i];
			if($temp > $max) {
				$max = $temp;
				$max_i = $i;
			}
		}

		$subgrad->[$max_i] += $c->[$max_i][$j];
	}

	return $subgrad;
}

1;
