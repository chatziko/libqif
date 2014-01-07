package QIF::Corners;
use strict;

use QIF::Channel;


# Constructor ($channel)
#   $channel: an object of Channel type
#
sub new {
	my ($class, $channel) = @_;
	my $self = { channel => $channel };
	return bless $self;
}


sub best_ratio {
	my ($self, $mode, $returnPoints) = @_;
	my $c = $self->{channel}{matrix};

	my $subsets = create_comb(2, scalar @$c);		# create all possible subsets
	$subsets = [grep sum(@$_) >= 2, @$subsets];		# keep only those with at least 2 elements

	my $ratio = 0;
	my $xr;				# the x that gives the ratio

	my $count = 0;
	my $countok = 0;

	if($returnPoints) {
		!$self->{'fork'}		or die 'cannot use $returnPoints in fork mode';
		$self->{points} = [];
	}

	# loop for all subsets if users (the rest are zero)
	for my $sub(@$subsets) {
		# create a new matrix with only the selected users
		my $cs = QIF::Matrix->new([]);
		for my $i(0..@$c-1) {
			push @$cs, $c->[$i]
				if $sub->[$i];
		}
		my $chans = QIF::Channel->new(matrix => $cs);

		my $n = @$cs;				# number of users
		my $m = @{ $cs->[0] };		# number of observables

		my ($lcount, $lcountok) = (0, 0);	# local counters
		my $lxr;							# local $xr (only for the subset of users)

		# rebuild the point by filling 0 for the removed users
		my $rebuild_point = sub {
			my ($lx) = @_;
			my $x = [];
			my $i = 0;
			push @$x, ($_ ? $lx->[$i++] : 0)
				for @$sub;
			return $x;
		};

		my $alist = create_perm($n);
		my $olist = create_comb($m, $n-1);

		my $pid;
		if($self->{'fork'}) {
			$pid = open(KID, "-|");
			die "can't fork" unless defined $pid;
			if($pid) {	# parent
				splice @$alist, @$alist/2;		# remove last half
			} else {
				splice @$alist, 0, @$alist/2;	# remove first half
			}
		}

		for my $al (@$alist) {
		for my $ol (@$olist) {
			$lcount++;

			my $xs = solve_system($cs, $al, $ol)	or next;
			check_solution($cs, $al, $ol, $xs)		or next;

			$self->store_point($rebuild_point->($xs))	if $returnPoints;
			$lcountok++;

			my $pe = $chans->prob_error($xs);
			my $ob = $mode == 1 ? $chans->santhi_vardy($xs) :
					 $mode == 2 ? $chans->hellman_raviv($xs) :
					 1;

			if($pe > $ratio * $ob) {
				$ratio = $pe/$ob;
				$lxr = $xs;
			}
		}}

		$count += $lcount;
		$countok += $lcountok;

		# rebuild $xr from $lxr
		$xr = $rebuild_point->($lxr)
			if $lxr;

		if($self->{'fork'}) {
			if($pid) {
				# get results from child process
				my @kid = <KID>; chomp for @kid;
				close(KID);

				$count += $kid[0];
				$countok += $kid[1];

				if($kid[2] > $ratio) {
					$ratio = $kid[2];
					@$xr = @kid[3..@kid-1];
				}

			} else {
				print $lcount, "\n";
				print $lcountok, "\n";
				print $ratio, "\n";
				print "$_\n" for @$xr;
				CORE::exit;
			}
		}
	}

	my $res = {
		ratio =>	$ratio,
		x =>		$xr,
		count =>	$count,
		countok =>	$countok,
		points =>	$self->{points}
	};
	return $res;
}

# TODO: restructure the code, for now we get the corner points through best_ratio
#
sub corner_points {
	my ($self) = @_;
	return $self->best_ratio(3, 1)->{points};
}

sub compare_vulnerabilities {
	my ($class, $c1, $c2) = @_;

	# vuln(c1) is convex, we check the inequality on the points of v(c2)
	#
	my $corn = Corners->new($c2);
	my $points = $corn->corner_points;

	for my $p (@$points) {
		my $v1 = $c1->vulnerability($p);
		my $v2 = $c2->vulnerability($p);
		return $p, $v1, $v2
			if $v1 > $v2;
	}
}

sub best_ratio_sym {
	my ($self, $returnPoints) = @_;

	my $ratio = 0;
	my $xr;				# the x that gives the ratio
	my $count = 0;

	$self->{points} = []
		if $returnPoints;

	my $chan = $self->{channel};
	my $m = $chan->{matrix};
	my $phi = $m->[0][0] / $m->[0][1];
	my $n = @$m;	# number of rows

	# number of non zero elements
	for my $i (2..$n) {
		# number of 'low' value elements. There must be at least one low and one high
		for my $j (0..$i) {
			$count++;

			my $low = 1 / ($j + ($i-$j)*$phi);
			my $high = $low * $phi;

			# construct point
			my $x = [];
			push @$x, ($high)x($i-$j);
			push @$x, ($low)x$j;
			push @$x, (0)x($n-$i);
			
			$self->store_point($x)		if $returnPoints;

			# check if it gives better ration
			my $pe = $chan->prob_error($x);
			my $ob = $chan->santhi_vardy($x);

			if($pe > $ratio * $ob) {
				$ratio = $pe/$ob;
				$xr = $x;
			}
		}
	}

	my $res = {
		ratio =>	$ratio,
		x =>		$xr,
		count =>	$count,
		countok =>	$count,
		points =>	$self->{points},
	};
	return $res;
}

sub best_ratio_exhaustive {
	my ($self) = @_;

	my $n = $self->{channel}->input_no;

	my $curx = [(.5) x $n];
	my ($r, $x);
	for(my $wid = .5; $wid > 0.01; $wid /= 10) {
		($r, $x) = $self->best_ratio_exhaustive_rec([], $curx, $wid);
		$curx = $x;
	}

	my $res = {
		ratio =>	$r,
		x =>		$x,
		count =>	0,
		countok =>	0,
	};
	return $res;
}

sub best_ratio_exhaustive_rec {
	my ($self, $px, $curx, $wid) = @_;

	# take a copy
	my $x = [@$px];
	my $n = @$x;

	# if solution is complete
	if($n == @$curx-1) {
		push @$x, 1-sum(@$x);

		my $pe = $self->{channel}->prob_error($x);
		my $ob = $self->{channel}->santhi_vardy($x);
		my $r = $ob > .0001 ? $pe/$ob : 0;

		return $r, $x;

	} else {
		my $r = 0;
		my $xr;

		my $from = max($curx->[$n] - $wid, 0);
		my $to = min($curx->[$n] + $wid, 1-sum(@$x));

		for(my $v = $from; $v <= $to; $v += $wid/40) {
			$x->[$n] = $v;
			my ($nr, $nx) = $self->best_ratio_exhaustive_rec($x, $curx, $wid);
			if($nr > $r) {
				$r = $nr;
				$xr = $nx;
			}
		}

		return $r, $xr;
	}
}

sub store_point {
	my ($self, $xs) = @_;

	my $copy = [@$xs];
	my $points = $self->{points};

	vector_eq($_, $copy) and return
		for @$points;

	push @$points, $copy;
}

sub vector_eq {
	my ($a, $b) = @_;
	@$a == @$b	or return 0;
	for(0..@$a-1) {
		# if numbers are Math::BigRat then use bcmp
		#
		return 0
			if ref($a->[$_])
				? $a->[$_]->bcmp($b->[$_])
				: abs($a->[$_] - $b->[$_]) > 0;
	}
	return 1;
}



# --- Auxiliary functions ---

sub vec_eq {
	my ($x, $y) = @_;
	for my $i (0..@$x-1) {
		return 0 if abs($x->[$i] - $y->[$i]) > .0001;
	}
	return 1;
}

sub check_solution {
	my ($c, $l, $k, $x) = @_;

	my $f;# = vec_eq($x, [0.166666666666667, 0.5, 0.166666666666667, 0.166666666666667]);
	if($f) {
		print "found\n";
		print "\t" . join(', ', @$l) . "\n";
		print "\t" . join(', ', @$k) . "\n";
	}

	my $failed = 0;

	for my $i(0..@$l-2) {
		my ($li, $ki) = ($l->[$i], $k->[$i]);
		my $t = $x->[$li] * $c->[$li][$ki] + .0001;

		for my $i2 (0..@$c-1) {
			unless($t > $x->[$i2] * $c->[$i2][$ki]) {
				printf "failed rows $li $i2 col $ki %f %f\n", $t, $x->[$i2] * $c->[$i2][$ki] if $f;
				$failed++;
				last;
			}
		}
	}
	print "succ\n" if $f;
	return $failed <= 0;
}

sub solve_system {
	my ($c, $l, $k) = @_;
	my $n = @$l;

	my @coeff = (1);		# coeff_i is the ratio of x_{l_i} wrt to x_{l_0}
	for my $i (1..$n-1) {
		my $a = $c->[ $l->[$i-1] ][ $k->[$i-1] ];
		my $b = $c->[ $l->[$i]   ][ $k->[$i-1] ];
		$b > 0	or return undef;		# no solution

		$coeff[$i] = $coeff[$i-1] *  $a/$b; 
	}

	my $x = [];
	$x->[$l->[0]] = 1/sum(@coeff);		# compute x_{l_0}

	for my $i (1..$n-1) {				# compute the rest based on x_{l_0}
		$x->[$l->[$i]] = $x->[$l->[0]] * $coeff[$i];
	}

	return $x;
}


# create_comb ($n, $l)
# returns all ordered combinations with repetition of numbers 0..$n-1, $l times
#
# eg create_comb(2, 2) will give [[0,0], [0,1], [1,0], [1,1]]
#
sub create_comb {
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

# create_perm ($n)
# returns all permutations of numbers 0..$n-1 where 0 is always first
#
sub create_perm {
	my ($n) = @_;

	my $perm = create_perm_rec([1..$n-1]);
	unshift @$_, 0	for @$perm;

	return $perm;
}
sub create_perm_rec {
	my ($array) = @_;
	my $n = @$array		or return [[]];
	my $perm = [];

	for my $i (0..$n-1) {
		my $temp = [ @$array ];			# copy array
		my $el = splice @$temp, $i, 1;	# remove element in pos $i

		my $rest = create_perm_rec($temp);	# permute the rest

		for(@$rest) {					# for all permutations of the remaining elements
			unshift @$_, $el;			# put back element
			push @$perm, $_;			# add permutation to the array
		}
	}
	return $perm;
}

sub sum {
	my $s = 0;
	$s += $_ for @_;
	return $s;
}

sub max {
	my $m = $_[0];
	for my $i (1..@_-1) {
		$m = $_[$i] if $_[$i] > $m;
	}
	return $m;
}

sub min {
	my $m = $_[0];
	for my $i (1..@_-1) {
		$m = $_[$i] if $_[$i] < $m;
	}
	return $m;
}

# return p(.|o_j)
sub post {
	my ($c, $x, $j) = @_;
	my $xo = [ map $x->[$_] * $c->[$_][$j], 0..@$c-1 ];		# $x .* the j-th column

	my $s = sum(@$xo) or return [map 0, @$c];				# if the sum is 0 then o_j has 0 prob
	$_ /= $s for @$xo;										# normalize to get probabilities

	return $xo;
}

1;
