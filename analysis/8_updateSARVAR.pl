#! /usr/bin/env perl

use strict;
use warnings;

use DateTime qw( );
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=  "Using sequence demix proportion files located in RUNDIR/7 DASHBOARD/, update the WaTCH-WV SARS variant database.\n";
$usage   .=   "\n";

my $rundir;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } else {
		$rundir = $arg;
	}
}

die "FATAL: $progname requires a valid run directory.\n$usage\n" unless defined $rundir and -d "$rundir";

#
# Get time stamp for now
#
my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');



my @files = glob(qq("$rundir/7 DASHBOARD/*.txt"));
#print Dumper(\@files);
#die;

my %metadata = ();
my %props    = ();


my $WATCHDB_MAIN  = "updates/watchdb.LATEST.txt";
my $SARVARDB_MAIN = "updates/sarvardb.LATEST.txt";
my $SARVARDB_INCR = "updates/sarvardb/sarvardb.$NOW.txt";


`cp $SARVARDB_MAIN $SARVARDB_MAIN.OLD`;

open my $fh, "<", "$SARVARDB_MAIN" or die "Unable to open $SARVARDB_MAIN: $!\n";
while (my $line = <$fh>) {
	next if "$line" =~ /^sample_id/;
	chomp $line;
	my ($at, $runid, $facility, $county, $start, $end, $varname, $prop) = split /\t/, "$line", -1;
	next unless defined $runid and $runid ne "";
	$props{"$at"} = {} unless defined $props{"$at"};
	$props{"$at"}->{"$runid"} = {} unless defined $props{"$at"}->{"$runid"};
	$props{"$at"}->{"$runid"}->{"$varname"} = $prop;
}
close $fh;


foreach my $f (@files) {
	open my $fh, "<", "$f" or die "Unable to open $f: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		my ($runid, $at, $varname, $prop) = split /\t/, "$line", -1;
		next unless defined $runid and $runid ne "";
		$props{"$at"} = {} unless defined $props{"$at"};
		$props{"$at"}->{"$runid"} = {} unless defined $props{"$at"}->{"$runid"};
		$props{"$at"}->{"$runid"}->{"$varname"} = $prop;
		$metadata{"$at"} = {} unless defined $metadata{"$at"};
	}
	close $fh;
}

#print Dumper(\%props);
#die;

open my $fhw, "<", "$WATCHDB_MAIN" or die "Unable to open $WATCHDB_MAIN: $!\n";
while (my $line = <$fhw>) {
	chomp $line;
	my @cols = split /\t/, "$line", -1;
	my ($at, $county, $wwtp, $start, $end) = ($cols[0], $cols[5], $cols[6], $cols[11], $cols[12]);
	if (defined $props{$at}) {
		$metadata{$at} = {
			"at" => "$at",
			"county" => "$county",
			"start_datetime" => "$start",
			"end_datetime" => "$end",
			"facility" => "$wwtp"
		};
	}
}
close $fhw;

#print Dumper(\%metadata);
#die;

open my $fho, ">", "$SARVARDB_INCR" or die "Unable to open $SARVARDB_INCR for writing: $!\n";
print $fho "sample_id\tseqrun_id\tfacility\tcounty\tstart_datetime\tend_datetime\tvariant\tproportion\n";
foreach my $at (keys %props) {
	foreach my $rid (keys %{$props{"$at"}}) {
		foreach my $varname (keys %{$props{"$at"}->{"$rid"}}) {
			print $fho "$at\t";
			print $fho "$rid\t";
			print $fho "$metadata{$at}->{facility}\t";
			print $fho "$metadata{$at}->{county}\t";
			print $fho "$metadata{$at}->{start_datetime}\t";
			print $fho "$metadata{$at}->{end_datetime}\t";
			print $fho "$varname\t";
			print $fho "$props{$at}->{$rid}->{$varname}\n";
		}
	}
}
close $fho;

`cp $SARVARDB_INCR $SARVARDB_MAIN`;

exit;

