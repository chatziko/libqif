#!/usr/bin/perl -w 
use strict; 

use Tk; 
use Tk::LabFrame;
use Tk::DialogBox;

use QIF::Channel;
use QIF::Corners;


my $defMatrix = 
"0.32961  0.32961  0.17039  0.17039
0.32961  0.17039  0.32961  0.17039
0.68482  0.07240  0.07240  0.17039
0.07240  0.68482  0.17039  0.07240
0.07240  0.17039  0.68482  0.07240
0.17039  0.07240  0.07240  0.68482";

my @maxValues = (3, 1, 2);
my %maxLabels = (
	3 => 'Probability of error (PE)',
	1 => 'Ratio PE/Santhi-Vardy bound',
	2 => 'Ratio PE/Hellman-Raviv bound',
);
my $maxSelected = 3;
my $fork = 0;


my $mw = MainWindow->new; 

$mw->Label(-text => "Channel matrix:")->pack(-anchor => 'w');

my $txt = $mw->Text(-width => 50, -height => 10);
$txt->Contents($defMatrix);
$txt->pack(-fill => 'both', -expand => 1);

my $frame = $mw->LabFrame(-label => "Compute the maximum:", -labelside => 'acrosstop',
                           )->pack(qw/-side left -fill none/);



$frame->Radiobutton(-text => $maxLabels{$_}, -variable => \$maxSelected, -value => $_,
	)->pack(qw/-anchor w -fill none -expand 0 -side top/)
    for @maxValues;

$mw->Checkbutton(-text => 'Use multiple processors', -variable => \$fork,
	)->pack(qw/-anchor sw -fill none -expand 0 -side bottom/);

my $button = $mw->Button(-text => 'GO', -command => \&go); 
$button->pack(-side => 'top', -fill => 'x', -expand => 1); 

my $button2 = $mw->Button(-text => 'Quit', -command => \&exit); 
$button2->pack(-side => 'top', -fill => 'x', -expand => 1); 

MainLoop; 


sub go {
	my $matrix = QIF::Matrix->new($txt->Contents);
	my $chan = QIF::Channel->new(matrix => $matrix);

	alert("Error", "Currently matrices with zero elements are not supported"), return
		if $chan->has_zero;

	my $corn = new QIF::Corners $chan;
	$corn->{'fork'} = $fork;

	my $start = time;
	my $res = $corn->best_ratio($maxSelected);

	my $s =
		"Maximum $maxLabels{$maxSelected}: $res->{ratio}\n\n" .
		"Distribution giving the maximum: " . join(', ', @{$res->{x}}) . "\n\n" .
		"Solved systems: $res->{count}\n\n" .
		"Time: " . (time-$start) . " seconds\n\n" .
		"Matrix:\n" . $chan->to_string;

	alert("Result", $s);
}

sub alert {
	my ($title, $msg) = @_;
    my $d = $mw->DialogBox(-title => $title, -buttons => ["OK"]);
	my $l = $d->add("Label", -text => $msg);
	$l->pack(-side => "left", -anchor => "w");
    my $button = $d->Show;
}


