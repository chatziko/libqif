# This module reports the prior that minimizes the posterior
# vulnerability of a channel
# Try it with examples/minvuln.pl
# It is just a slight modification of the module Compare, in future
# we'll find a better place to put it.
# I changed the sign of the coefficients of the objective functions,
# so that the linear problem will minimize it. The value of the
# objective function is negated as well, that's why we start with a
# max of -2. The prior obtained is just fine.

package QIF::Minvuln;
use Moose;

use QIF::Matrix qw/:all/;
use QIF::Channel;
use QIF::LinearProgram;


has c1	=> (is => 'rw', isa => 'QIF::Channel', required => 1);
# has c2	=> (is => 'rw', isa => 'QIF::Channel', required => 1);
has Gs	=> (is => 'rw', isa => 'ArrayRef[QIF::Matrix]', default => sub { [] });
has solved => (is => 'rw', isa => 'Num', required => 0, default => 0);


sub compute {
	my ($self) = @_;

	for my $G (@{$self->Gs}) {
		my $res = $self->_minimize_vuln($G);
		return $res, $G if $res;
	}
	return undef;
}

sub _minimize_vuln {
	my ($self, $G) = @_;

	my $M1 = $self->c1->matrix;
#	my $M2 = $self->c2->matrix;

	my $Y1 = $M1->cols;
#	my $Y2 = $M2->cols;
	my $K = $G->cols; # |W| in the paper

	# we generate all possible combinations of "best guesses" for each column
	# for each combination we generate the corresponding linear program, if it's feasible and the optimal is negative
	# then C1 has greater vulnerability for that solution
	#
	my $combs = _create_comb($K, $Y1);
	my $n = 0;

	my $max = -2;		# record the biggest vulnerability difference seen so far, and the prior that causes it.
	my $bestprior;

	for my $comb (@$combs) {
		$n++;
		my $k_per_y1 = [splice @$comb, 0, $Y1];		# first $Y1 elements are for C1
#		my $k_per_y2 = $comb;						# remaining $Y2 elements are for C1

		my $program = $self->_build_program($G, $k_per_y1);
		my ($z, $x) = $program->solve;
                
		# continue in case of no solution or a non-record solution
#                if (defined $z){ print "[found : " . $z . " with prior " . $x . "]\n\n"; }
		defined $z && !_less_than_or_eq($z, $max)		or next;

#		print "[found a min: " . $z . " with prior " . $x . "]\n\n";
		$max = $z;
		$bestprior = $x;

		# precaution check: test that the answer is the same as the vulnerability difference
		my $v1 = $self->c1->g_vulnerability($x, $G);
#		my $v2 = $self->c2->g_vulnerability($x, $G);
		_equal($v1, -$z) 		or die "$z is not the same as the leakage difference (".($v1).")";

		#last;		# commented this out so that we try *all* combinations
	}
	$self->solved($n);

	return $bestprior;
}


# Same as above but returns all the results, so that we can check if
# there are more than one minimum and if the corresponding priors
# provides better prior vulnerability. 
# This method should just be called by OCaml.

sub compute_all {
	my ($self) = @_;

	for my $G (@{$self->Gs}) {
		my (@priors, @vulns) = $self->_minimize_vuln_all($G);
		return @priors, @vulns
	}
	return undef;
}

sub _minimize_vuln_all {
	my ($self, $G) = @_;

	my $M1 = $self->c1->matrix;
#	my $M2 = $self->c2->matrix;

	my $Y1 = $M1->cols;
#	my $Y2 = $M2->cols;
	my $K = $G->cols; # |W| in the paper

	# we generate all possible combinations of "best guesses" for each column
	# for each combination we generate the corresponding linear program, if it's feasible and the optimal is negative
	# then C1 has greater vulnerability for that solution
	#
	my $combs = _create_comb($K, $Y1);
	my $n = 0;

	my @vulns;		# record the biggest vulnerability difference seen so far, and the prior that causes it.
	my @priors;

	for my $comb (@$combs) {
		$n++;

		my $k_per_y1 = $comb;

		my $program = $self->_build_program($G, $k_per_y1);
		my ($z, $x) = $program->solve;
                
		# continue in case of no solution or a non-record solution
#                if (defined $z){ print "[found : " . $z . " with prior " . $x . "]\n\n"; }
#		defined $z && !_less_than_or_eq($z, $max)		or next;
		defined $z 		or next;
#		print "[found a min: " . $z . " with prior " . $x . "]\n\n";
		push (@vulns, $z);
		push (@priors, $x);

		# precaution check: test that the answer is the same as the vulnerability difference
		my $v1 = $self->c1->g_vulnerability($x, $G);
		_equal($v1, -$z) 		or die "$z is not the same as the leakage difference (".($v1).")";

		#last;		# commented this out so that we try *all* combinations
	}
	$self->solved($n);

	return (\@priors,\@vulns);
}




sub add_function {
	my ($self, $G) = @_;
	push @{$self->Gs}, $G;
	return $G;
}

sub add_identity {
	my ($self) = @_;
	push @{$self->Gs}, QIF::Matrix->identity($self->c1->matrix->rows);
}

# generates all 2-block partition gain functions
#
sub add_2_block {
	my ($self) = @_;

	my $M1 = $self->c1->matrix;
	my $X = $M1->rows;
	my $Gs = $self->Gs;

	# create all possible 2-block partitions of inputs. To break the symmetry, we always put the first
	# input in the first partition
	#
	my $partitions = _create_comb(2, $X-1);	# combinations of the form [0,1,0,1], 1 means belongs to the first partition

	for my $part (@$partitions) {
		grep { !$_ } @$part		or next;	# ignore [1,1,...1]
		unshift @$part, 1;					# always take first input

		push @$Gs, QIF::Matrix->new([
			map({ $_ == 1 ? [1, 0] : [0, 1] } @$part)
		]);
	}
	return 1;
}

# generates all 0-1 gain functions
#
sub add_01_functions {
	my ($self) = @_;

	my $M1 = $self->c1->matrix;
	my $X = $M1->rows;
	my $Gs = $self->Gs;

	# create all possible subsets of inputs, except the empty and full subset
	#
	my $subsets = _create_comb(2, $X);	# combinations of the form [0,1,0,1], 1 means belongs to the subset
	@$subsets = grep { grep({ $_ > 0 } @$_) && grep({ $_ == 0 } @$_) } @$subsets;

	# we take all combinations of subsets with the following rules
	#  - none should be a subset of any other
	#  - we should have at least 2 subsets
	#
	my $k = @$subsets;
	my $combs = _create_comb(2, $k);

	for my $comb (@$combs) {
		my $subs = [ map { $subsets->[$_] } grep { $comb->[$_] } 0..$k-1 ];
		$subs = _remove_subsets($subs);
		@$subs >= 2		or next;

		my $G = QIF::Matrix->new($subs)->transpose;
		#!grep { $G->row($_)->is_zero } 1..$G->rows	or next;	# no zero rows (WE SHOULD PROVE THAT THIS IS SAFE
		!grep { $G == $_ } @$Gs						or next;	# might be already there because we remove the subsets

		push @$Gs, $G;
	}

	# sort based on the number of guesses, cause it's faster to check the ones with fewer guesses
	@$Gs = sort { $a->cols <=> $b->cols } @$Gs;
}


sub _build_program {
	my ($self, $G, $k_per_y1) = @_;

	my $M1 = $self->c1->matrix;
#	my $M2 = $self->c2->matrix;

	my $X = $M1->rows;
	my $Y1 = $M1->cols;
#	my $Y2 = $M2->cols;

	# 3 things to build
	my $eq   = QIF::Matrix->new([]);
	my $ineq = QIF::Matrix->new([]);
	my $obj  = QIF::Matrix->new([]);

	# inequalities
	#
	push @$ineq, @{ $self->_build_ineq_for_c($G, $self->c1, $k_per_y1) };
#	push @$ineq, @{ $self->_build_ineq_for_c($G, $self->c2, $k_per_y2) };

	# equality to express that the sum of values is 1
	#
	push @$eq, [(1)x($X), -1];		# x_1 + ... + x_n - 1 = 0

	# objective function
	#   V(C1) - V(C2) = sum_x p_x ( sum_y1 M1_x,y1 G_x,kb - sum_y2 M2_x,y2 G_x,kb )
	#   V(C1) = sum_x p_x ( sum_y1 M1_x,y1 G_x,kb)
	# if positive the C1 has greater leakage
	#
	my @coeff;
	for my $x (0..$X-1) {
		my $c = 0;
		$c += $M1->[$x][$_] * $G->[$x][ $k_per_y1->[$_] ]	for 0..$Y1-1;
#		$c -= $M2->[$x][$_] * $G->[$x][ $k_per_y2->[$_] ]	for 0..$Y2-1;

		push @coeff, $c;
	}
	push @coeff, 0;		# constant coeff 0
        @coeff = map {-$_} @coeff; # changing sign we _minimize_ the objective function
	push @$obj, \@coeff;

	# add one more inequality asking that the objective function is positive
	#
	#push @$ineq, [map {-$_} @coeff];		# negate to turn <= into >=

	# done
	return QIF::LinearProgram->new(objective => $obj, equalities => $eq, inequalities => $ineq);
}


sub _build_ineq_for_c {
	my ($self, $G, $c, $k_per_y) = @_;

	my $M = $c->matrix;

	my $X = $M->rows;
	my $Y = $M->cols;
	my $K = $G->cols;

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

			my @coeff;		# vector containing the coefficients of the inequality
			for my $x (0..$X-1) {
				push @coeff, $M->[$x][$y] * ($G->[$x][$k] - $G->[$x][$kb]);
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

# $a, $b are 0-1 arrays seen as sets
# returns true if $a \subseteq $b
sub _subseteq {
	my ($a, $b) = @_;
	for(0..@$a-1) {
		return 0 if $a->[$_] && !$b->[$_];
	}
	return 1;
}

# $sets is array of 0-1 arrays (seen as sets)
# removes the sets that are subseteq of some other set
# (by consequence it also removes duplicates)
#
sub _remove_subsets {
	my ($sets) = @_;

	# remove each set to avoid comparing it with itself
	my @res;
	while(my $set = shift @$sets) {
		push @res, $set		unless grep { _subseteq($set, $_) } @$sets, @res;
	}
	return \@res;
}

1;
