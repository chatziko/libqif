use QIF::Matrix;
use QIF::Channel;

use List::Util qw/min max/;

while(1) {

my $p1 = rand();
my $p2 = rand();
#$p1 =0.00241107807519469;
#$p2 =0.0024195291276996;

$a = QIF::Matrix->new([[$p1, 1-$p1]]);
$b = QIF::Matrix->new([[$p2, 1-$p2]]);

my $r1 = abs($p1 - $p2) / min( max($p1, $p2), max(1-$p1, 1-$p2) );

my $r2 = QIF::Channel->bounded_entropy_distance($a, $b);

die "$p1\n$p2\n$r1\n$r2\n"	if abs($r1-$r2) > 0.0001;

}
