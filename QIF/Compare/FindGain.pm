package QIF::Compare::FindGain;
use Moose;

use QIF::Channel;
use QIF::LinearProgram;


has c1	=> (is => 'rw', isa => 'QIF::Channel', required => 1);
has c2	=> (is => 'rw', isa => 'QIF::Channel', required => 1);
has prior => (is => 'rw', isa => 'QIF::Matrix', required => 0);
has solved => (is => 'rw', isa => 'Num', required => 0, default => 0);


my $epsilon = 1e-6;

# $K is the number of guesses
sub compare {
	my ($self, $K) = @_;
	$K	or die "must give the number of guesses";

	my $M1 = $self->c1->matrix;
	my $M2 = $self->c2->matrix;

	my $X = $M1->rows;
	my $Y1 = $M1->cols;
	my $Y2 = $M2->cols;

	my $p = $self->prior || $self->prior(QIF::Matrix->uniform($X));

	# we generate all possible combinations of "best guesses" for each column
	# for each combination we generate the corresponding linear program, if it's feasible and the optimal is negative
	# then C1 has greater vulnerability for that solution
	#
	my $combs = _create_comb($K, $Y1+$Y2);
	my $G;
	my $n = 0;
	for my $comb (@$combs) {
		$n++;
		my $k_per_y1 = [splice @$comb, 0, $Y1];		# first $Y1 elements are for C1
		my $k_per_y2 = $comb;						# remaining $Y2 elements are for C2

		my $program = $self->_build_program($K, $k_per_y1, $k_per_y2);
		my ($z, $sol) = $program->solve;

		# continue in case of no solution or 0 solution
		defined $z && !_equal($z, 0)		or next;

		# build $G from the solution
		$G = QIF::Matrix->empty;	
		for my $x (0..$X-1) {
			for my $k (0..$K-1) {
				$G->[$x][$k] = $sol->[0][$x*$K+$k];
			}
		}

		# precaution check: test that the answer is the same as the vulnerability difference
		my $v1 = $self->c1->g_vulnerability($p, $G);
		my $v2 = $self->c2->g_vulnerability($p, $G);
		_equal($v1 - $v2, $z) 		or die "$z is not the same as the leakage difference (".($v1-$v2).")";

		last;
	}
	$self->solved($n);

	return $G;
}

sub _build_program {
	my ($self, $K, $k_per_y1, $k_per_y2) = @_;

	my $M1 = $self->c1->matrix;
	my $M2 = $self->c2->matrix;
	my $p = $self->prior;

	my $X = $M1->rows;
	my $Y1 = $M1->cols;
	my $Y2 = $M2->cols;

	# 3 things to build
	my $eq   = QIF::Matrix->new([]);
	my $ineq = QIF::Matrix->new([]);
	my $obj  = QIF::Matrix->new([]);

	# inequalities
	#
	push @$ineq, @{ $self->_build_ineq_for_c($K, $self->c1, $k_per_y1) };
	push @$ineq, @{ $self->_build_ineq_for_c($K, $self->c2, $k_per_y2) };

	# inequalities to express that G's elements are <= 1
	#
	for my $i (0..$X*$K-1) {
		# build x_i - 1 <= 0
		my $vec = [(0)x($X*$K), -1];
		$vec->[$i] = 1;
		push @$ineq, $vec;
	}

	# TEST: cols of G add up to 1
#	for my $k (0..$K-1) {
#		my @coeff = ((0)x($X*$K), -1);
#		for my $x (0..$X-1) {
#			$coeff[ $x*$K+$k ] = 1;
#		}
#		push @$eq, \@coeff;
#	}

	# objective function
	#   V(C1) - V(C2) = sum_x p_x ( sum_y1 M1_x,y1 G_x,kb - sum_y2 M2_x,y2 G_x,kb )
	# if positive the C1 has greater leakage
	#
	my @coeff = (0)x($X*$K);
	for my $x (0..$X-1) {
		for my $y (0..$Y1-1) {
			my $k = $k_per_y1->[$y];
			$coeff[ $x*$K+$k ] += $p->[0][$x] * $M1->[$x][$y];
		}
		for my $y (0..$Y2-1) {
			my $k = $k_per_y2->[$y];
			$coeff[ $x*$K+$k ] -= $p->[0][$x] * $M2->[$x][$y];
		}
	}
	push @coeff, 0;		# constant coeff 0
	push @$obj, \@coeff;

	# add one more inequality asking that the objective function is positive
	#
	push @$ineq, [map {-$_} @coeff];		# negate to turn <= into >=

	# done
	return QIF::LinearProgram->new(objective => $obj, equalities => $eq, inequalities => $ineq);
}


sub _build_ineq_for_c {
	my ($self, $K, $c, $k_per_y) = @_;

	my $M = $c->matrix;
	my $p = $self->prior;

	my $X = $M->rows;
	my $Y = $M->cols;

	@$k_per_y == $Y		or die "wrong k_per_y";

	# for each y let kb be the best guess, stored in $k_per_y.
	# for each k != kb we create one inequality of the form
	#   sum_x p_x M_x,y (G_x,k - G_x,kb) <= 0
	#
	my $ineq = QIF::Matrix->new([]);

	for my $y(0..$Y-1) {
		my $kb = $k_per_y->[$y];

		for my $k (0..$K-1) {		# for each k != kb
			next if $k == $kb;

			my @coeff = (0)x($X*$K);		# vector containing the coefficients of the inequality
			for my $x (0..$X-1) {
				$coeff[ $x*$K+$k  ] += $p->[0][$x] * $M->[$x][$y];
				$coeff[ $x*$K+$kb ] -= $p->[0][$x] * $M->[$x][$y];
			}
			push @coeff, 0;		# constant coeff 0
			push @$ineq, \@coeff;
		}
	}

	return $ineq;
}

# create_comb ($n, $l)
# returns all ordered combinations with repetition of numbers 0..$n-1, $l times
#
# eg create_comb(2, 2) will give [[0,0], [0,1], [1,0], [1,1]]
#
sub _create_comb {
	my ($n, $l) = @_;
	my $perm = [[]];

	for my $i (1..$l) {					# add one element each time
		my $temp = [];
		for(@$perm) {					# for all existing permutation
			for my $j (0..$n-1) {			# for all numbers
				push @$temp, [@$_, $j];	# extend existing permutation, put back in $temp
			}
		}
		$perm = $temp;
	}
	return $perm;
}

sub _equal {
	my ($a, $b) = @_;
	return $QIF::Matrix::USE_RAT
		? $a == $b
		: $a > $b - $epsilon && $a < $b + $epsilon;
}

1;
