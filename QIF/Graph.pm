package QIF::Graph;
use Moose;
use Moose::Util::TypeConstraints;
use List::Util qw/min sum/;

use QIF::Matrix;

use overload fallback => 1,
	'""' => \&stringify;

has adjacency	=> (is => 'rw', isa => 'QIF::Matrix', required => 1, coerce => 1);
has layout		=> (is => 'rw', isa => 'Str', required => 0, );
has labels		=> (is => 'rw', isa => 'ArrayRef', required => 0, );
has metric	=> (is => 'rw', isa => 'QIF::Matrix', required => 0, lazy_build => 1 );


sub line {
	my ($class, $n) = @_;
	return $class->_undirected_from_edge_test($n, sub {
		$_[1] - $_[0] == 1
	});
}

sub radius {
	my ($class, $n, $r) = @_;
	return $class->_undirected_from_edge_test($n, sub {
		$_[1] - $_[0] <= $r
	});
}

sub cycle {
	my ($class, $n) = @_;
	my $g = $class->_undirected_from_edge_test($n, sub {
		$_[1] - $_[0] == 1 ||
		$_[1] - $_[0] == $n - 1
	});
	$g->layout('circo');
	return $g;
}

sub grid {
	my ($class, $n, $m) = @_;
	return $class->line($n)->cartesian_product( $class->line($m) );
}

sub hamming {
	my ($class, $d, $q) = @_;
	# cartesian product of $d complete graphs of $q elements
	my $c = $class->complete($q);
	my $h = $c;
	$h = $h->cartesian_product($c)	for 1..$d-1;
	return $h;
}

sub hypercube {
	my ($class, $d) = @_;
	return $class->hamming($d, 2);
}

sub complete {
	my ($class, $n) = @_;
	my $g = $class->_undirected_from_edge_test($n, sub { 1 });
	$g->layout('circo');
	return $g;
}

sub star {
	my ($class, $n) = @_;
	return $class->_undirected_from_edge_test($n, sub {
		$_[0] == 0
	});
}

sub random {
	my ($class, $n, $p) = @_;
	$p ||= 0.5;
	my $g = $class->_undirected_from_edge_test($n, sub {
		rand() < $p
	});
	$g->layout('circo');
	return $g;
}

# not really efficient
sub histogram {
	my ($class, $n, $bins) = @_;

	# nodes are tuples with the counts summing to $n
	my @nodes = grep { sum(@$_) == $n } @{ _create_comb($n+1, $bins) };
	my @labels = map { join(',', @$_) } @nodes;

	my $g = $class->_undirected_from_edge_test(\@nodes, sub {
		my ($x, $y) = @_;

		# must differ in exactly 2 elements, the one should be +1 and the other -1
		my @diff = grep { $x->[$_] != $y->[$_] } 0..@$x-1;
		@diff == 2	or return 0;
		my ($i, $j) = @diff;

		return
			($x->[$i] == $y->[$i] + 1 && $x->[$j] == $y->[$j] - 1) ||
			($x->[$i] == $y->[$i] - 1 && $x->[$j] == $y->[$j] + 1);
	}, labels => \@labels);
	$g->layout('circo');
	return $g;
}

sub _undirected_from_edge_test {
	my ($class, $nodes, $test, %opt) = @_;

	# $nodes can either be an array of arbitrary scalars representing nodes, or a number
	# in which case the nodes list will be 0..$nodes-1
	#
	$nodes = [0..$nodes-1]	unless ref $nodes;

	my @adj;
	for my $i (0..@$nodes-1) {
		for my $j ($i..@$nodes-1) {
			$adj[$i][$j] =
			$adj[$j][$i] =
				$i == $j || !$test->($nodes->[$i], $nodes->[$j]) ? 0 : 1;
		}
	}
	return $class->new(adjacency => QIF::Matrix->new(\@adj), %opt);
}

sub nodes {
	my ($self) = @_;
	return $self->adjacency->rows;
}

sub _build_metric {
	my ($self) = @_;
	my $n = $self->nodes;

	# Floyd-Warshall algorithm
	# $d is initialiazed to the cost to go from i to j. 0 if i=j, 1 if there's a path,
	# maximum ($n+1) otherwise
	#
	my $d = $self->adjacency->map(sub { $_[0] ? 1 : $n+1 });
	$d->[$_][$_] = 0	for 0..$n-1;

	for my $k (0..$n-1) {
	for my $i (0..$n-1) {
	for my $j (0..$n-1) {
		my $z = $d->[$i][$k] + $d->[$k][$j];
		$d->[$i][$j] = $z	if $d->[$i][$j] > $z;
	}}}

	!$d->any(sub { $_[0] == $n+1 })		or die "not connected";

	return $d;
}


######## GRAPH OPERATIONS ####
#
#
sub strong_product {
	my ($self, $g) = @_;

	my $A1 = $self->adjacency;
	my $n1 = $self->nodes;
	my $A2 = $g->adjacency;
	my $n2 = $g->nodes;

	my $lab1 = $self->labels || [0..$n1-1];
	my $lab2 = $g->labels    || [0..$n2-1];
	my @labels = map { "$_->[0],$_->[1]" } @{_cartesian($lab1, $lab2)};

	return ref($self)->_undirected_from_edge_test(_cartesian($n1, $n2), sub {
		my ($x1, $x2) = @{$_[0]};
		my ($y1, $y2) = @{$_[1]};
		return ($x1 == $y1 || $A1->[$x1][$y1]) &&
			   ($x2 == $y2 || $A2->[$x2][$y2]);
	}, labels => \@labels);
}

sub cartesian_product {
	my ($self, $g) = @_;

	my $A1 = $self->adjacency;
	my $n1 = $self->nodes;
	my $A2 = $g->adjacency;
	my $n2 = $g->nodes;

	my $lab1 = $self->labels || [0..$n1-1];
	my $lab2 = $g->labels    || [0..$n2-1];
	my @labels = map { "$_->[0],$_->[1]" } @{_cartesian($lab1, $lab2)};

	return ref($self)->_undirected_from_edge_test(_cartesian($n1, $n2), sub {
		my ($x1, $x2) = @{$_[0]};
		my ($y1, $y2) = @{$_[1]};
		return ($x1 == $y1 && $A2->[$x2][$y2]) ||
			   ($x2 == $y2 && $A1->[$x1][$y1]);
	}, labels => \@labels);
}


# For the moment this assumes an undirected graph
sub to_png {
	my ($self, $filename) = @_;

	require GraphViz;

	my $n = $self->nodes;
	my $A = $self->adjacency;

	my $graph;
	my $labels = $self->labels || [0..$n-1];
	my $layout = $self->layout || 'dot';
	#$layout = 'grid';

	if($layout eq 'line') {
		# force nodes on a line
		$graph = GraphViz->new(directed => 0, layout => 'neato', ovarlap => 0, splines => 1);
		my $z = int(sqrt($n));
		$graph->add_node($_, label => $labels->[$_], pos=> "$_,0!", width => 0.1)
			for 0..$n-1;

	} elsif($layout eq 'grid') {
		# for products. assumes square grid
		$graph = GraphViz->new(directed => 0, layout => 'neato');
		my $z = int(sqrt($n));
		$graph->add_node($_, label => $labels->[$_], pos=> int($_ / $z).','.($_ % $z).'!', width => 0.1)
			for 0..$n-1;

	} else {
		#create a new graph
		$graph = GraphViz->new(directed => 0, layout => $layout);
		$graph->add_node($_, label => $labels->[$_], width => 0.1)
			for 0..$n-1;
	}

	# create edges as defined by your adjacecncy matrix
	for my $i (0..$n-1) {
		for my $j ($i+1..$n-1) {
			$graph->add_edge($i => $j)
				if $A->[$i][$j];
		}
	}

	#render the graph to a png-file
	$graph->as_png($filename);
}

sub stringify {
	my ($self) = @_;
	return "".$self->adjacency;
}

# $n, $m can be arrayrefs, or numbers
sub _cartesian {
	my ($n, $m) = @_;
	$n = [0..$n-1]	unless ref $n;
	$m = [0..$m-1]	unless ref $m;
	my @res;
	for my $i (@$n) {
		for my $j (@$m) {
			push @res, [$i, $j];
		}
	}
	return \@res;
}

# _create_comb ($n, $l)
# returns all ordered combinations with repetition of numbers 0..$n-1, $l times
#
# eg _create_comb(2, 2) will give [[0,0], [0,1], [1,0], [1,1]]
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

1;
