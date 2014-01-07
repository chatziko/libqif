package QIF::LinearProgram;
use Moose;

use QIF::Matrix;

use overload fallback => 1,
	'""' => \&_stringify;

# maximization problem
# equalities are = 0    (use negated last col for constant)
# inequalities are <= 0 (use negated last col for constant)

has equalities		=> (is => 'rw', isa => 'QIF::Matrix', required => 1, coerce => 1);
has inequalities	=> (is => 'rw', isa => 'QIF::Matrix', required => 1, coerce => 1);
has objective		=> (is => 'rw', isa => 'QIF::Matrix', required => 1, coerce => 1);
has non_negative	=> (is => 'rw', isa => 'Bool', default => 1);


sub solve {
	my ($self) = @_;

	my $model = $self->generate_mps_model;

	my $file = '/tmp/' . int(rand()*1e8);
	open F, "> $file";
	print F $model;
	close F;

	my $output = $QIF::Matrix::USE_RAT
		? `exlp --print-solution --max $file`
		: `glpsol -m $file --freemps --max -o /dev/stdout`;
	unlink $file;

	#print $model, $output;
	return $QIF::Matrix::USE_RAT
		? $self->parse_exlp_output($output)
		: $self->parse_glpsol_output($output);
}

sub parse_exlp_output {
	my ($self, $output) = @_;

	die 'unbounded' if $output =~ /unbounded/;
	return undef if $output =~ /infeasible/;

	# feasible, parse solution
	my $z;
	my $x = [];
	for(split /[\r\n]/, $output) {
		if(/^optimum\s*:\s*([-\d\/]+)/) {
			$z = Math::BigRat->new($1);
		}
		if(/^X(\d+):\s+([-\d\/]+)/) {
			$x->[$1] = Math::BigRat->new($2);
		}
	}
	defined $z && grep({ defined } @$x) == $self->objective->cols - 1
		or die "cannot parse exlp output:\n\n$output";

	return $z, QIF::Matrix->new([$x]);
}

sub parse_glpsol_output {
	my ($self, $output) = @_;

	die 'unbounded' if $output =~ /UNBOUNDED/;
	return undef if $output =~ /INFEASIBLE|PROBLEM HAS NO PRIMAL FEASIBLE/;

	# feasible, parse solution
	my $z;
	my $x = [];
	for(split /[\r\n]/, $output) {
		if(/Objective:\s*OBJ = ([-+\d.e]+)/) {
			$z = $1;
		}
		if(/\s+\d+\s+X(\d+)\s+\w+\s+([-+\d.e]+)/) {
			$x->[$1] = $2;
		}
	}
	defined $z && grep({ defined } @$x) == $self->objective->cols - 1
		or die "cannot parse glpsol output:\n\n$output";

	return $z, QIF::Matrix->new([$x]);
}

sub generate_mps_model {
	my ($self) = @_;
	my $obj = $self->objective;
	my $eq = $self->equalities;
	my $ineq = $self->inequalities;
	my $X = $obj->cols;

	$obj->[0][$X-1] == 0		or die "constant not supported in objective function";

	my $model = "NAME  QIF\n";

	# ROWS
	$model .= "ROWS\n N  OBJ\n";
	$model .= " E  EQ$_\n"	for 0..@$eq-1;
	$model .= " L  INEQ$_\n"	for 0..@$ineq-1;

	# COLUMS
	$model .= "COLUMNS\n";
	for my $col (0..$X-2) {
		$model .= " X$col  OBJ  $obj->[0][$col]\n";

		$model .= " X$col  EQ$_  $eq->[$_][$col]\n"
			for 0..@$eq-1;

		$model .= " X$col  INEQ$_  $ineq->[$_][$col]\n"
			for 0..@$ineq-1;
	}

	# RHS
	$model .= "RHS\n";

	$model .= " RHS  EQ$_  " . -$eq->[$_][$X-1] . "\n"
		for 0..@$eq-1;

	$model .= " RHS  INEQ$_  " . -$ineq->[$_][$X-1] . "\n"
		for 0..@$ineq-1;

	unless($self->non_negative) {
		# No BOUNDS assumes >= 0
		$model .= "BOUNDS\n";
		$model .= " FR  BND  X$_\n"
			for 0..$X-2;
	}

	$model .= "ENDATA\n";

	return $model;
}

sub _stringify {
	my ($self) = @_;
	return
		"objective:\n" . $self->objective . "\n" .
		"equalities:\n" . $self->equalities . "\n" .
		"inequalities:\n" . $self->inequalities . "\n";
}

## not used anymore, glpsol supports mps
##
#sub generate_glpsol_model {
#	my ($self) = @_;
#	my $obj = $self->objective;
#	my $X = $obj->cols;
#
#	my $model = '';
#
#	# variable declaration (note: last column is constant, on variable)
#	#
#	for (0..$X-2) {
#		$model .= "var X$_ >= 0;\n";
#	}
#	$model .= "\n";
#
#	# objective function
#	#
#	$model .= "maximize OBJ: " . vector_to_expr($obj->[0]). ";\n\n";
#
#	# equalities
#	#
#	my $c = 0;
#	for(@{$self->equalities}) {
#		$model .= 'c' . $c++ . ': ' . vector_to_expr($_) . " = 0;\n";
#	}
#	$model .= "\n\n";
#
#	# inequalities
#	#
#	for(@{$self->inequalities}) {
#		$model .= 'c' . $c++ . ': ' . vector_to_expr($_) . " <= 0;\n";
#	}
#	$model .= "\n\n";
#
#	return $model;
#}
#
## turns a vector into a linear expression of variables
#sub vector_to_expr {
#	my ($row) = @_;
#	my $N = @$row;
#	return
#		join ' + ',
#			 map({ $row->[$_].'*X'.$_ } grep { $row->[$_] } 0..$N-2),
#			 $row->[$N-1];
#}


1;
