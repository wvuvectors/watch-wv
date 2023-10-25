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


my $warn_level = 2;

my $WATCHDB_MAIN  = "updates/watchdb.LATEST.txt";
my $SARVARDB_MAIN = "updates/sarvardb.LATEST.txt";
my $SARVAR_VTW    = "resources_seq/VTW.txt";

=cut
Columns in the sarvarDB:
	0	sample_id
	1	seqrun_id
	2	facility
	3	county
	4	start_datetime
	5	end_datetime
	6 lineage
	7	variant
	8	proportion
=cut

my %sarvar      = ("by_county" => {}, "by_facility" => {}, "by_run" => {}, "by_variant" => {}, "by_month" => {});
my %sample_ids  = ();
my %last_sample = ("yr" => 2020, "mo" => 1, "day" => 1);
my %last_run    = ("yr" => 2020, "mo" => 1, "day" => 1);

my %vtw = ();
open my $wfh, "<", "$SARVAR_VTW" or die "Unable to open $SARVAR_VTW: $!\n";
while (my $line = <$wfh>) {
	next if "$line" =~ /^#/ or "$line" =~ /^WHO_label/;
	chomp $line;
	my ($who_label, $lineage, $cdc_status, $watch_status) = split /\t/, "$line", -1;
	$vtw{"$lineage"} = {"cdc_status" => "$cdc_status", "label" => "$who_label", "watch_status" => "$watch_status"};
}
close $wfh;


open my $fh, "<", "$SARVARDB_MAIN" or die "Unable to open $SARVARDB_MAIN: $!\n";
while (my $line = <$fh>) {
	next if "$line" =~ /^sample_id/;
	chomp $line;
	my ($sample_id, $runid, $facility, $county, $start, $end, $parental, $lineage, $varname, $prop) = split /\t/, "$line", -1;
	
	my ($mo, $day, $yr) = (0,0,0);
	if ($end =~ m/(\d+?)\/(\d+?)\/(\d+?)\s/) {
		($mo, $day, $yr) = ($1, $2, $3);
	}
	$yr = "20" . "$yr" unless length "$yr" > 2;
	my $datecomp = "$yr" . "$mo" . "$day";
	
	unless (defined $sample_ids{"$sample_id"}) {
		
		$sample_ids{"$sample_id"} = 1;

		if ($yr >= $last_sample{"yr"} and $mo >= $last_sample{"mo"} and $day > $last_sample{"day"}) {
			$last_sample{"yr"} = $yr;
			$last_sample{"mo"} = $mo;
			$last_sample{"day"} = $day;
		}
		
		my ($rmo, $rday, $ryr) = (0,0,0);
		if ($runid =~ m/SARS_(\d{4})(\d{2})(\d{2})/) {
			($ryr, $rmo, $rday) = ($1, $2, $3);
		}
		if ($ryr >= $last_run{"yr"} and $rmo >= $last_run{"mo"} and $rday > $last_run{"day"}) {
			$last_run{"yr"} = $ryr;
			$last_run{"mo"} = $rmo;
			$last_run{"day"} = $rday;
		}
		
		$sarvar{"by_sample_month"}->{"$yr-$mo"} = 0 unless defined $sarvar{"by_sample_month"}->{"$yr-$mo"};
		$sarvar{"by_sample_month"}->{"$yr-$mo"} = $sarvar{"by_sample_month"}->{"$yr-$mo"}+1;
		
		$sarvar{"by_county"}->{"$county"} = 0 unless defined $sarvar{"by_county"}->{"$county"};
		$sarvar{"by_county"}->{"$county"} = $sarvar{"by_county"}->{"$county"}+1;

		$sarvar{"by_facility"}->{"$facility"} = 0 unless defined $sarvar{"by_facility"}->{"$facility"};
		$sarvar{"by_facility"}->{"$facility"} = $sarvar{"by_facility"}->{"$facility"}+1;
	
		$sarvar{"by_run"}->{"$runid"} = 0 unless defined $sarvar{"by_run"}->{"$runid"};
		$sarvar{"by_run"}->{"$runid"} = $sarvar{"by_run"}->{"$runid"}+1;
	}
	
	if (defined $vtw{"$varname"}) {
		$vtw{"$varname"}->{"by_county"} = {} unless defined $vtw{"$varname"}->{"by_county"};
		$vtw{"$varname"}->{"by_county"}->{"$county"} = $datecomp unless defined $vtw{"$varname"}->{"by_county"}->{"$county"};
		$vtw{"$varname"}->{"by_county"}->{"$county"} = $datecomp if $datecomp > $vtw{"$varname"}->{"by_county"}->{"$county"};
	} elsif (defined $vtw{"$lineage"}) {
		$vtw{"$lineage"}->{"by_county"} = {} unless defined $vtw{"$lineage"}->{"by_county"};
		$vtw{"$lineage"}->{"by_county"}->{"$county"} = $datecomp unless defined $vtw{"$lineage"}->{"by_county"}->{"$county"};
		$vtw{"$lineage"}->{"by_county"}->{"$county"} = $datecomp if $datecomp > $vtw{"$lineage"}->{"by_county"}->{"$county"};
	} elsif (defined $vtw{"$parental"}) {
		$vtw{"$parental"}->{"by_county"} = {} unless defined $vtw{"$parental"}->{"by_county"};
		$vtw{"$parental"}->{"by_county"}->{"$county"} = $datecomp unless defined $vtw{"$parental"}->{"by_county"}->{"$county"};
		$vtw{"$parental"}->{"by_county"}->{"$county"} = $datecomp if $datecomp > $vtw{"$parental"}->{"by_county"}->{"$county"};
	}
	
	$sarvar{"by_variant"}->{"$varname"} = 0 unless defined $sarvar{"by_variant"}->{"$varname"};
	$sarvar{"by_variant"}->{"$varname"} = $sarvar{"by_variant"}->{"$varname"}+1;

	$sarvar{"by_lineage"}->{"$lineage"} = 0 unless defined $sarvar{"by_lineage"}->{"$lineage"};
	$sarvar{"by_lineage"}->{"$lineage"} = $sarvar{"by_lineage"}->{"$lineage"}+1;

	$sarvar{"by_parent"}->{"$parental"} = 0 unless defined $sarvar{"by_parent"}->{"$parental"};
	$sarvar{"by_parent"}->{"$parental"} = $sarvar{"by_parent"}->{"$parental"}+1;

}
close $fh;

#print Dumper(\%vtw);
#die;


# Print summary info to STDOUT
print "#\n# Summary of sarvarDB from $NOW.\n#\n\n";
print "Total num sequenced : " . scalar(keys %sample_ids) . ".\n";
print "Date of last sample : " . $last_sample{"mo"} . "/" . $last_sample{"day"} . "/" . $last_sample{"yr"} . ".\n";
print "Last sequencing run : " . $last_run{"mo"} . "/" . $last_run{"day"} . "/" . $last_run{"yr"} . ".\n";

print "\n";
print "Distribution of variants to watch across all samples (WaTCH level <= $warn_level, sorted by most recent detection):\n\n";
my @vtw_sorted = sort {$vtw{"$a"} cmp $vtw{"$b"}} keys %vtw;
foreach my $lineage (@vtw_sorted) {
	next if $vtw{"$lineage"}->{"watch_status"} > $warn_level;
	if (defined $vtw{"$lineage"}->{"by_county"}) {
		print "$vtw{$lineage}->{label} variant $lineage found in " . scalar(keys %{$vtw{"$lineage"}->{"by_county"}}) . " count";
		if (scalar(keys %{$vtw{"$lineage"}->{"by_county"}}) == 1) {
			print "y:\n";
		} else {
			print "ies:\n";
		}
		foreach my $county (sort({$vtw{"$lineage"}->{"by_county"}->{"$b"} <=> $vtw{"$lineage"}->{"by_county"}->{"$a"}} keys %{$vtw{"$lineage"}->{"by_county"}})) {
			my $datestr = $vtw{"$lineage"}->{"by_county"}->{"$county"};
			$datestr =~ s/(\d{4})(\d{2})(\d{2})/$2\/$3\/$1/i;
			print "   $county, most recent detection on $datestr.\n";
		}
		print "\n";
	}
}

=cut
print "\n";
print "Top 20 most abundant parental groups across all samples:\n\n";
my @parentals_sorted = sort {$sarvar{"by_parent"}->{"$b"} <=> $sarvar{"by_parent"}->{"$a"}} keys %{$sarvar{"by_parent"}};
my $i = 0;
foreach my $parental (@parentals_sorted) {
	$i++;
	print "$sarvar{by_parent}->{$parental} samples contain $parental\n";
	last if $i == 20;
}

print "\n";
print "Top 20 most abundant lineages across all samples:\n\n";
my @lineages_sorted = sort {$sarvar{"by_lineage"}->{"$b"} <=> $sarvar{"by_lineage"}->{"$a"}} keys %{$sarvar{"by_lineage"}};
$i = 0;
foreach my $lineage (@lineages_sorted) {
	$i++;
	print "$sarvar{by_lineage}->{$lineage} samples contain $lineage\n";
	last if $i == 20;
}
=cut

print "\n";
print "Top 20 most abundant variants across all samples:\n\n";
my @variants_sorted = sort {$sarvar{"by_variant"}->{"$b"} <=> $sarvar{"by_variant"}->{"$a"}} keys %{$sarvar{"by_variant"}};
my $i = 0;
foreach my $variant (@variants_sorted) {
	$i++;
	print "$sarvar{by_variant}->{$variant} samples contain $variant\n";
	last if $i == 20;
}

print "\n";
print "Overall there are " . scalar(keys %{$sarvar{"by_variant"}}) . " variants comprising " . scalar(keys %{$sarvar{"by_lineage"}}) . " lineages in the database.\n";
print "\n# EOF.\n";


exit;

