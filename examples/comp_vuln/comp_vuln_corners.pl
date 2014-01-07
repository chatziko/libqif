#!/usr/bin/perl -w
use strict;

use Channel;
use Corners;

use CGI qw/:standard/;



if(param("action")) {

	my $rat = 1;

	my $m1 = param("m1")		or error("no matrix");
	my $m2 = param("m2")		or error("no matrix");

	my $c1 = Channel->new($m1, rat => $rat);
	my $c2 = Channel->new($m2, rat => $rat);

	my ($p, $v1, $v2) = Corners->compare_vulnerabilities($c1, $c2);

	print header, "<html><body>";
	if($p) {
		print "C1 has greater vulnerability ($v1) than C2 ($v2) for the prior:<br>" .
				join(" ", @$p), "\n";
	} else {
		print "C1 has always less vulnerability than C2";
	}

	print "<br><br><a href='javascript:history.back()'>Back</a></body></html>";

} else {

	print header, q{
		<html><body>
		Give two matrices to compare their vulnerability (and by consequence their leakage) over all
		priors. Enter one row per line, and elements separated by whitespace. You can use fractions
		(eg 1/2) for elements. If you select the fractional arithmetic option, the computation becomes slower but
		the results are exact.
		<br><br>
		<form method="POST" action="comp_vuln.pl">
		C1:<br>
		<textarea name="m1" rows=5 cols=20>
1/3 1/3 1/3   0
  0 1/3 1/3 1/3
1/3 2/3   0   0
		</textarea>
		<br><br>
		C2:<br>
		<textarea name="m2" rows=5 cols=20>
1/2   0   0 1/2
  0 1/2 1/2   0
1/2 1/2   0   0
		</textarea>
		<br><bR>
		<input type="submit" name="action" value="Go"/>
		</form>
		</bodu></html>
	};
}



sub error {
	my ($s) = @_;
	print header, "<html><body>$s</body></html>";
	exit;
}
