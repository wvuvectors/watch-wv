#! /usr/bin/env perl

use strict;
use warnings;

use DateTime qw( );
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options]\n";
$usage   .=  "Print summary info about the WaTCH-WV SARS variant database to STDOUT.\n";
$usage   .=  "   [-c FILE] Print county data in a tabular format amenable to R.\n";
$usage   .=  "   [-f FILE] Print facility data in a tabular format amenable to R.\n";
$usage   .=  "   [-r FILE] Print run data in a tabular format amenable to R.\n";
$usage   .=  "   [-v FILE] Print variant data in a tabular format amenable to R.\n";
$usage   .=  "   [-d FILE] Print date data in a tabular format amenable to R.\n";
$usage   .=   "\n";

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
	}
}

#
# Get time stamp for now
#
my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%m/%d/%Y');



my $WATCHDB_MAIN  = "updates/watchdb.LATEST.txt";
my $SARVARDB_MAIN = "updates/sarvardb.LATEST.txt";

=cut
Columns in the sarvarDB:
	0	sample_id
	1	seqrun_id
	2	facility
	3	county
	4	start_datetime
	5	end_datetime
	6	variant
	7	proportion
=cut

my %sarvar      = ("by_county" => {}, "by_facility" => {}, "by_run" => {}, "by_variant" => {}, "by_month" => {});
my %sample_ids  = ();
my %last_sample = ("yr" => 2020, "mo" => 1, "day" => 1, "run" => "");

open my $fh, "<", "$SARVARDB_MAIN" or die "Unable to open $SARVARDB_MAIN: $!\n";
while (my $line = <$fh>) {
	next if "$line" =~ /^sample_id/;
	chomp $line;
	my ($sample_id, $runid, $facility, $county, $start, $end, $varname, $prop) = split /\t/, "$line", -1;
	
	my $lineage = "$varname";
	if ($varname =~ /^XBB/i) {
		$lineage =~ s/XBB\.(.+?)\.(.+?)\..*$/XBB.$1.$2/i;
	} else {
		$lineage =~ s/(.+?)\.(.+?)\..*$/$1.$2/i;
	}
	
	unless (defined $sample_ids{"$sample_id"}) {
		
		$sample_ids{"$sample_id"} = 1;
		my ($mo, $day, $yr) = (0,0,0);
		if ($end =~ m/(\d+?)\/(\d+?)\/(\d+?)\s/) {
			($mo, $day, $yr) = ($1, $2, $3);
		}
		
		if ($yr >= $last_sample{"yr"} and $mo >= $last_sample{"mo"} and $day > $last_sample{"day"}) {
			$last_sample{"yr"} = $yr;
			$last_sample{"mo"} = $mo;
			$last_sample{"day"} = $day;
			$last_sample{"run"} = "$runid";
		}
		
		
		$sarvar{"by_month"}->{"$yr-$mo"} = 0 unless defined $sarvar{"by_month"}->{"$yr-$mo"};
		$sarvar{"by_month"}->{"$yr-$mo"} = $sarvar{"by_month"}->{"$yr-$mo"}+1;
		
		$sarvar{"by_county"}->{"$county"} = 0 unless defined $sarvar{"by_county"}->{"$county"};
		$sarvar{"by_county"}->{"$county"} = $sarvar{"by_county"}->{"$county"}+1;

		$sarvar{"by_facility"}->{"$facility"} = 0 unless defined $sarvar{"by_facility"}->{"$facility"};
		$sarvar{"by_facility"}->{"$facility"} = $sarvar{"by_facility"}->{"$facility"}+1;
	
		$sarvar{"by_run"}->{"$runid"} = 0 unless defined $sarvar{"by_run"}->{"$runid"};
		$sarvar{"by_run"}->{"$runid"} = $sarvar{"by_run"}->{"$runid"}+1;
	}

	$sarvar{"by_variant"}->{"$varname"} = 0 unless defined $sarvar{"by_variant"}->{"$varname"};
	$sarvar{"by_variant"}->{"$varname"} = $sarvar{"by_variant"}->{"$varname"}+1;

	$sarvar{"by_lineage"}->{"$lineage"} = 0 unless defined $sarvar{"by_lineage"}->{"$lineage"};
	$sarvar{"by_lineage"}->{"$lineage"} = $sarvar{"by_lineage"}->{"$lineage"}+1;

}
close $fh;

#print Dumper(\%sarvar);
#die;

# summary info to STDOUT
print "#\n# Summary of sarsvarDB from $NOW.\n#\n\n";
print "Total num sequenced : " . scalar(keys %sample_ids) . ".\n";
print "Date of last sample : " . $last_sample{"mo"} . "/" . $last_sample{"day"} . "/" . $last_sample{"yr"} . ".\n";
print "Last sequencing run : " . $last_sample{"run"} . ".\n";
print "\n";
print "Top 20 most abundant lineages across all samples:\n\n";

my @lineages_sorted = sort {$sarvar{"by_lineage"}->{"$b"} <=> $sarvar{"by_lineage"}->{"$a"}} keys %{$sarvar{"by_lineage"}};
my $i = 0;
foreach my $lineage (@lineages_sorted) {
	$i++;
	print "$sarvar{by_lineage}->{$lineage} samples contain $lineage\n";
	last if $i == 20;
}
print "\n";
print "Overall there are " . scalar(keys %{$sarvar{"by_variant"}}) . " variants comprising " . scalar(keys %{$sarvar{"by_lineage"}}) . " lineages in the database.\n";
print "\n# EOF.\n";

=cut
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
=cut


exit;

