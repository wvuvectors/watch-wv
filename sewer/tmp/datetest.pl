#! /usr/bin/env perl

use strict;
use warnings;

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use DateTime;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] [DAYS]\n";
$usage   .=   "Build update table for AssetTiger database. Pass an integer DAYS to determine ";
$usage   .=   "how many days of samples to include in the update (default = 90). Writes to STDOUT.\n";
$usage   .=   "\n";

my $dbdir  = "data/latest";

my $status = 0;
my $days = 90;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } else {
		$days = $arg;
		$days =~ s/\D//gi;
	}
}


#
# Get time stamp for now
#
my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');

my $THRESHOLD_DATE = DateTime->now->clone->subtract( days => $days );

my %table = ();

my $linenum = 0;
my @headers = ();
open (my $tFH, "<", "$dbdir/watchdb.sample.txt") or die "Unable to open $dbdir/watchdb.sample.txt for reading: $!\n";
while (my $line = <$tFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;

	my @cols = split "\t", "$line", -1;

	if ($linenum == 0) {
		$linenum++;
		for (@cols) {
			push @headers, "$_";
		}
	} else {
		# first field is always the uid regardless of table
		my $uid = $cols[0];

		$table{"$uid"} = {};
		for (my $i=0; $i < scalar(@cols); $i++) {
			$table{"$uid"}->{"$headers[$i]"} = trim("$cols[$i]");
		}
	}
}
close $tFH;

#print Dumper(\%table);
#die;



foreach my $sample_id (keys %table) {
	my $resultRef = $table{"$sample_id"};
	
	# AT exports times at AM/PM!
	
	# Skip unless sample_collection_datetime >= TODAY - days
	my $dt_string = $resultRef->{"sample_collection_start_datetime"};	# MM-DD-YYYY HH:MM:SS
	next if !defined $dt_string or "$dt_string" eq "";
	
	my ($d_str, $t_str) = split /\s/, "$resultRef->{sample_collection_start_datetime}", -1;
	my ($mo, $dd, $yyyy) = split /\//, "$d_str", -1;
	#my ($hh, $mi, $ss) = split /:/, "$t_str", -1;
	unless (defined $mo and defined $dd and defined $yyyy) {
		print "$sample_id\t$dt_string\n";
	}

	my $dt = 
		 DateTime->new(
				time_zone => 'local',
				year   => $yyyy,
				month  => $mo,
				day    => $dd,
				hour   => 0,
				minute => 0,
				second => 0
	);
				
	unless (DateTime->compare($dt, $THRESHOLD_DATE) == -1) {
		print "$sample_id\t$dt_string\n";
	}
	

}


sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


