package QIF::Matrix;
use Moose::Util::TypeConstraints;

# A bit of a hack: this is not a proper Moose class, cause we want to use a
# blessed arrayref directly. But we define a 'QIF::Matrix' type, and a coercion from Str
#

use overload fallback => 1,
	'""' => \&_stringify,
	'*' => \&multiply, 
	'==' => \&is_equal_to;

use Math::BigRat lib => 'GMP';
use List::Util ();

our $USE_RAT = 0;
our $CONVERT_FRACTIONS = 1;

class_type 'QIF::Matrix';
coerce 'QIF::Matrix'
	=> from 'Str'
	=> via { parse($_) };

my $epsilon = 1e-4;
our $display_precision = 4;

use Exporter 'import';
our @EXPORT_OK = qw/_equal _less_than _less_than_or_eq/;
our %EXPORT_TAGS = ('all' => \@EXPORT_OK);


######### CONSTRUCTORS ###########
#
sub new {
	my ($class, $m) = @_;

	# if string, parse
	return parse($m)	unless ref $m;

	# convert elements of the form 'n/m', convert to rats if needed
	if($USE_RAT || $CONVERT_FRACTIONS) {
		for(@$m) {
			for(@$_) {
				if($USE_RAT) {
					$_ = Math::BigRat->new($_)	unless ref $_;
				} elsif($CONVERT_FRACTIONS) {
					$_ = $1/$2	if /(\d+)\/(\d+)/;
				}
			}
		}
	}

	return bless $m;
}

sub empty {
	my ($class) = @_;
	return $class->new([]);
}

sub identity {
	my ($class, $n) = @_;
	
	my @m = map { [(0)x$_, 1, (0)x($n-$_-1)] } 0..$n-1;

	return $class->new(\@m);
}

sub random {
	my ($class, $m, $n) = @_;
	return $class->new([
		map { [ map { rand() } 1..$n ] }
			1..$m
	]);
}

sub ones {
	my ($class, $m, $n) = @_;
	return $class->constant($m, $n, 1);
}

sub zeroes {
	my ($class, $m, $n) = @_;
	return $class->constant($m, $n, 0);
}

sub constant {
	my ($class, $m, $n, $c) = @_;
	return $class->new([
		map { [ map { $c } 1..$n ] }
			1..$m
	]);
}

sub uniform {
	my ($class, $n) = @_;
	my $p = $USE_RAT ? Math::BigRat->new("1/$n") : 1/$n;
	return $class->new([[($p) x $n]]);
}


######## NEW MATRIX FROM EXISTING ONE ####
#
#
sub map {
	my ($self, $f) = @_;
	my $m = [];
	for my $x (0..$self->rows-1) {
		for my $y (0..$self->cols-1) {
			$m->[$x][$y] = $f->($self->[$x][$y]);
		}
	}
	return QIF::Matrix->new($m);
}

sub map_rows {
	my ($self, $f) = @_;
	my $m = [];
	for my $x (0..$self->rows-1) {
		$m->[$x] = $f->($self->[$x]);
	}
	return QIF::Matrix->new($m);
}

sub clone {
	my ($self) = @_;
	return $self->map(sub {
		return $USE_RAT ? $_[0]->copy : $_[0];
	});
}

sub sign {
	my ($self) = @_;
	return $self->map(sub {
		return $_[0] < 0 ? -1 : 1
	});
};

sub scalar_powers {
	my ($self, $p) = @_;
	return $self->map(sub { $p ** $_[0] });
}



######## ROW/COL MANIPULATION #############
#

sub rows {
	return scalar @{$_[0]};
}

sub cols {
	return @{$_[0]} ? scalar @{ $_[0]->[0] } : 0;
}

sub row {
	my ($self, $i) = @_;
	my $m = QIF::Matrix->new([]);
	$m->[0] = $self->[$i-1];
	return $m;
}

sub column {
	my ($self, $j) = @_;
	my $m = QIF::Matrix->new([]);
	for (0..$self->rows-1) {
		$m->[$_][0] = $self->[$_][$j-1];
	}
	return $m;
}

sub as_list {
	my ($self) = @_;
	return map { @$_ } @$self;
}

sub sum {
	my ($self) = @_;
	return List::Util::sum($self->as_list);
}

sub max {
	my ($self) = @_;
	return List::Util::max($self->as_list);
}



######## MATRIX PROPERTIES #########
#
sub has_non_negative_row {
	my ($self) = @_;
	for my $row (@$self) {
		grep { _less_than($_, 0) } @$row
			or return 1;
	}
	return 0;
}

sub all {
	my ($self, $f) = @_;
	for ($self->as_list) {
		$f->($_)	or return 0;	
	}
	return 1;
}

sub any {
	my ($self, $f) = @_;
	for ($self->as_list) {
		$f->($_)	and return 1;
	}
	return 0;
}

sub is_zero {
	return shift->all(sub { _equal($_[0], 0) });
}

sub is_negative {
	return shift->all(sub { _less_than($_[0], 0) });
}

sub is_non_negative {
	return shift->all(sub { !_less_than($_[0], 0) });
}

sub is_stochastic {
	my ($self) = @_;
	$self->is_non_negative	or return 0;

	for (1..$self->rows) {
		#warn $_;
		_equal($self->row($_)->sum, 1)	or return 0;
	}
	return 1;
}

sub is_equal_to {
	my ($self, $M) = @_;

	my $I = $self->rows;
	my $J = $self->cols;
	$M->rows == $I && $M->cols == $J	or return 0;

	for my $i (0..$I-1) {
		for my $j (0..$J-1) {
			_equal($self->[$i][$j], $M->[$i][$j])	or return 0;
		}
	}
	return 1;
}



######## MATRIX OPERATIONS #########
#

sub transpose {
	my ($self) = @_;
	my $m = [];

	for my $x (0..$self->rows-1) {
		for my $y (0..$self->cols-1) {
			$m->[$y][$x] = $self->[$x][$y];
		}
	}
	return QIF::Matrix->new($m);
}

sub multiply {
	my ($self, $m) = @_;
	my $res = QIF::Matrix->new([]);
	for my $x (0..$self->rows-1) {
		for my $z (0..$m->cols-1) {
			my $w = 0;
			for my $y (0..$self->cols-1) {
				$w += $self->[$x][$y] * $m->[$y][$z];
			}
			$res->[$x][$z] = $w;
		}
	}
	return $res;
}

# generate the inverse.
#
sub inverse {
	my ($self) = @_;
	my $n = $self->rows;
	$self->cols == $n		or die 'not square';

	if($USE_RAT) {
		return $self->_inverse_lp;
	} else {
		return $self->_inverse_armadillo;
		#return $self->_inverse_scilab;

		# old method, use Math::MatrixReal
		#my $inv = $self->_to_math_matrixreal->inverse	or die "not invertible";
		#return QIF::Matrix->_from_math_matrixreal($inv);
	}
}

sub _inverse_lp {
	my ($self) = @_;
	my $n = $self->rows;

	# when using rats, we employ a linear program (but without actually optimizing anything)
	#
	my $eq   = QIF::Matrix->empty;
	my $obj  = QIF::Matrix->empty;

	push @$obj, [(0)x($n*$n), 0];		# objective function is 0, we don't really optimize

	for my $i (0..$n-1) {
		for my $j (0..$n-1) {
			my @coeff = (0)x($n*$n);
			for my $k (0..$n-1) {
				$coeff[ $k*$n+$j ] = $self->[$i][$k];
			}
			push @coeff, ($i == $j ? -1 : 0);		# constant coeff 0 or 1
			push @$eq, \@coeff;
		}
	}

	my $program = QIF::LinearProgram->new(objective => $obj, equalities => $eq, inequalities => QIF::Matrix->empty, non_negative => 0);
	my ($z, $sol) = $program->solve;
	$sol		or die "not invertible";

	# build inverse from the solution
	my $inverse = QIF::Matrix->empty;	
	for my $i (0..$n-1) {
		for my $j (0..$n-1) {
			$inverse->[$i][$j] = $sol->[0][$i*$n+$j];
		}
	}
	return $inverse;
}

sub _inverse_scilab {
	my ($self) = @_;
	my $n = $self->rows;

	my $tempfile = '/tmp/matrix';
	open F, "> $tempfile";
	print F join(' ', @$_), "\n"
		for @{$self->_to_float_array};
	close F;
	system 'cp /tmp/matrix /tmp/a';

	my $command = qq{scilab-cli -e 'a = read("$tempfile", $n, $n); mdelete("$tempfile"); write("$tempfile", inv(a)); exit;'};
	my $output = `$command 2>&1`;
	$? == 0		or die "scilab failed: $output $!";

	open F, "< $tempfile";
	$/ = undef;
	my $s = <F>;
	close F;
	unlink $tempfile;

	return QIF::Matrix->new($s);
}

sub _inverse_armadillo {
	my ($self) = @_;
	my $n = $self->rows;

	my $r = int(rand()*1000000);
	my $tempfile_out = "/tmp/matrix_${r}.mat";
	my $tempfile_in = "/tmp/matrix_${r}_inv.mat";

	open F, "> $tempfile_out";
	print F "ARMA_MAT_TXT_FN008\n$n $n\n";
	print F join(' ', @$_), "\n"
		for @{$self->_to_float_array};
	close F;

	my $command = qq{./invert $tempfile_out $tempfile_in};
	my $output = `$command 2>&1` // '';
	if($? != 0) {
		warn "out of memory"	if $! =~ /allocate memory/;
		die "invert failed: $output $!";
	}

	open F, "< $tempfile_in";
	<F>; <F>;
	local $/ = undef;
	my $s = <F>;
	close F;

	unlink $tempfile_out, $tempfile_in;

	return QIF::Matrix->new($s);
}

sub right_inverse {
	my ($self) = @_;
	my $t = $self->transpose;
	return $t * ($self * $t)->inverse;
}

sub left_inverse {
	my ($self) = @_;
	my $t = $self->transpose;
	return ($t * $self)->inverse * $t;
}


# Multiply each row by the least common multiple of all denominators, to remove fractions
#
sub remove_fractions {
	my ($self) = @_;

	for my $row (@$self) {
		my $n = 1;
		for(@$row) {
			ref $_ or next;		# not a fraction
			my $den = $_->denominator;
			$n = $den * $n / _gcd($den, $n);		# LCM of $den and $n
		}
		$_ *= $n	for @$row;
	}
}

# replace small values by a real 0
sub cleanup_zeroes {
	my ($self) = @_;
	for my $row (@$self) {
		for(@$row) {
			$_ = 0	if _equal($_, 0);
		}
	}
	return $self;
}


##### AUXILIARY OBJECT METHODS
#
sub _to_float_array {
	my ($self) = @_;
	return $self	unless $USE_RAT;

	my $m = [];
	for my $x (0..$self->rows-1) {
		for my $y (0..$self->cols-1) {
			$m->[$x][$y] = $self->[$x][$y]->numify;
		}
	}
	return $m;
}

sub _to_math_matrixreal {
	my ($self) = @_;

	my $m = $self->_to_float_array;

	require Math::MatrixReal;
	return Math::MatrixReal->new_from_rows($m);
}

sub _format_number {
	$USE_RAT ? $_[0] : sprintf("%.${display_precision}f", $_[0])
}

sub _stringify {
	my ($self) = @_;

	my @width =
		map { List::Util::max(map { length(_format_number($_)) }
								  $self->column($_)->as_list)
			}
			1..$self->cols;

	my $s = '';
	for my $x (0..$self->rows-1) {
		for my $y (0..$self->cols-1) {
			my $a = _format_number($self->[$x][$y]);
			$s .= ($y ? ' ' : '') . $a . ' 'x($width[$y]-length($a));
		}
		$s .= "\n";
	}
	return $s;
}


######## AUXILIARY CLASS METHODS #####
#
sub _from_math_matrixreal {
	my ($class, $A) = @_;
	my ($M, $N) = $A->dim;
	my $B = [];
	for my $m (0..$M-1) {
		for my $n (0..$N-1) {
			$B->[$m][$n] = $A->element($m+1, $n+1);
		}
	}
	return $class->new($B);
}



######### AUXILIARY BARE METHODS #######
#

sub parse {
	my ($s) = @_;
	my $m = [];
	for my $line (split /[\r\n]/, $s) {
		$line =~ s/^\s+//g;
		$line =~ s/\s+$//g;
		$line ne ''		or next;

		push @$m, [split /\s+/, $line];
	}
	return QIF::Matrix->new($m);
}

sub _gcd {
	my ($a, $b) = @_;
	($a,$b) = ($b,$a) if $a > $b;
	while ($a) {
		($a, $b) = ($b % $a, $a);
	}
	return $b;
}

sub _equal {
	my ($a, $b) = @_;
	return $USE_RAT
		? $a == $b
		: $a > $b - $epsilon && $a < $b + $epsilon;
}

sub _less_than {
	my ($a, $b) = @_;
	return $USE_RAT
		? $a < $b
		: $a < $b - $epsilon;
}

sub _less_than_or_eq {
	my ($a, $b) = @_;
	return $USE_RAT
		? $a <=$b
		: $a < $b + $epsilon;
}


# Returns a matrix compatible with Algorithm::Simplex, containing either Num's or Math::Cephes::Fraction's
#
#sub to_simplex_tableau {
#	my ($self) = @_;
#
#	my $m = [];
#	for my $x (0..$self->rows-1) {
#		for my $y (0..$self->cols-1) {
#			my $w = $self->[$x][$y];
#			$m->[$x][$y] = $USE_RAT
#				? Math::Cephes::Fraction->new(ref $w ? $w->parts : ($w, 1))
#				: $w;
#		}
#	}
#	return $m;
#}

1;
