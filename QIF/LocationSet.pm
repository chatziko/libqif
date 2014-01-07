package QIF::LocationSet;
use Moose;

use List::Util qw/max/;

my $USE_KNN;
eval {
	require Algorithm::KNN::XS;
	$USE_KNN = 1;
};

has points	=> (is => 'rw', isa => 'ArrayRef[Any]', required => 1, );
has metric	=> (is => 'rw', isa => 'QIF::Matrix', required => 0, lazy_build => 1 );
has knn		=> (is => 'rw', isa => 'Algorithm::KNN::XS', required => 0, lazy_build => 1 );
has hash	=> (is => 'rw', isa => 'HashRef', required => 0, lazy_build => 1 );


sub grid {
	my ($class, $width, $height, $step_x, $step_y) = @_;
	$step_y //= $step_x;

	my @points;
	for my $i (0..$width-1) {
		for my $j (0..$height-1) {
			push @points, [$i * $step_x, $j * $step_y];
		}
	}

	return $class->new(points => \@points);
}

sub _build_metric {
	my ($self) = @_;

	my $d;
	my $points = $self->points;

	for my $i (0..@$points-1) {
		$d->[$i][$i] = 0;

		for my $j ($i+1..@$points-1) {
			$d->[$i][$j] = $d->[$j][$i] = euclidean($points->[$i], $points->[$j]);
		}
	}
	return QIF::Matrix->new($d);
}

sub _build_knn {
	my ($self) = @_;
	return Algorithm::KNN::XS->new(points => $self->points);
}

sub _build_hash {
	my ($self) = @_;
	my $points = $self->points;
	return {
		map { $points->[$_][0].",".$points->[$_][1] => $_ } 0..@$points-1
	};
}

sub euclidean {
	my ($a, $b) = @_;
	my $d1 = $a->[0] - $b->[0];
	my $d2 = $a->[1] - $b->[1];
	return sqrt($d1*$d1 + $d2*$d2);
}

sub closest {
	my ($self, $p) = @_;
	my $closest = 0;

	if($USE_KNN) {
		my $res = $self->knn->annkSearch(
			query_point		=> $p,
			limit_neighbors	=> 1,
			epsilon			=> 0,
		)->[0]{point};
		$closest = $self->hash->{$res->[0].",".$res->[1]};

	} else {
		my $points = $self->points;
		my $mind = euclidean($p, $points->[0]);

		for my $i (1..@$points-1) {
			my $d = euclidean($p, $points->[$i]);
			if($d < $mind) {
				$closest = $i;
				$mind = $d;
			}
		}
	}
	defined $closest	or die;
	return $closest;
}

sub gain_function {
	my ($self) = @_;
	my $m = $self->metric;
	my $max = max map { max(@$_) } @$m;
	return $m->map(sub { 1 - $_[0]/$max });
}


1;
