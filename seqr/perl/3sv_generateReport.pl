#! /usr/bin/env perl

use strict;
use warnings;

use DateTime qw( );
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options]\n";
$usage   .=  "Print summary report on the WaTCH-WV seqr database to STDOUT.\n";
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

my $WATCHDB_MAIN  = "../patchr/data/latest/watchdb.sample.txt";
my $SARVARDB_MAIN = "data/latest/seqrdb.txt";
my $SARVAR_VTW    = "resources/vtws.txt";

=cut
Columns in the seqrDB table:
	0	sample_id
	1	sbatch_id
	2	location
	3	county
	4	sample_collection_start_datetime
	5	sample_collection_end_datetime
	6 lineage_group
	7	variant
	8	variant_proportion
=cut

my %sarvar      = ("by_county" => {}, "by_location" => {}, "by_run" => {}, "by_variant" => {}, "by_month" => {});
my %sample_ids  = ();
my $last_sample = 0;
#my %last_sample = ("yr" => 2020, "mo" => 1, "day" => 1);
my %last_run_h  = ("yr" => 2020, "mo" => 1, "day" => 1);

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
	my ($sample_id, $sbatch_id, $location, $county, $start, $end, $lineage, $varname, $prop) = split /\t/, "$line", -1;
	
	my ($mo, $day, $yr) = (0,0,0);
	if ($end =~ m/(\d+?)\/(\d+?)\/(\d+?)\s/) {
		($mo, $day, $yr) = ($1, $2, $3);
	}
	$yr = "20" . "$yr" unless length "$yr" > 2;
	my $datecomp = "$yr" . "$mo" . "$day";
	
	unless (defined $sample_ids{"$sample_id"}) {
		
		$sample_ids{"$sample_id"} = {"sbatch_id" => "$sbatch_id", "print_date" => "$mo/$day/$yr", "sort_date" => "$datecomp", "county" => "$county", "location" => "$location"};
		
		$last_sample = $datecomp if $datecomp > $last_sample;
		
		my ($rmo, $rday, $ryr) = (0,0,0);
		if ($sbatch_id =~ m/SARS_(\d+)/) {
			my $rdate = "$1";
			($ryr, $rmo, $rday) = (substr($rdate, 0, 4), substr($rdate, 4, 2), substr($rdate, 6, 2));
		}

		if ($ryr > $last_run_h{"yr"}) {
			%last_run_h = ("yr" => $ryr, "mo" => $rmo, "day" => $rday);
		} elsif ($ryr == $last_run_h{"yr"}) {
			if ($rmo > $last_run_h{"mo"}) {
				%last_run_h = ("yr" => $ryr, "mo" => $rmo, "day" => $rday);
			} elsif ($rmo == $last_run_h{"mo"} and $rday > $last_run_h{"day"}) {
				%last_run_h = ("yr" => $ryr, "mo" => $rmo, "day" => $rday);
			}
		}
		
		
#		$sarvar{"by_sample_month"}->{"$yr-$mo"} = 0 unless defined $sarvar{"by_sample_month"}->{"$yr-$mo"};
#		$sarvar{"by_sample_month"}->{"$yr-$mo"} = $sarvar{"by_sample_month"}->{"$yr-$mo"}+1;
		
		$sarvar{"by_county"}->{"$county"} = 0 unless defined $sarvar{"by_county"}->{"$county"};
		$sarvar{"by_county"}->{"$county"} = $sarvar{"by_county"}->{"$county"}+1;

		$sarvar{"by_location"}->{"$location"} = 0 unless defined $sarvar{"by_location"}->{"$location"};
		$sarvar{"by_location"}->{"$location"} = $sarvar{"by_location"}->{"$location"}+1;
	
		$sarvar{"by_run"}->{"$sbatch_id"} = 0 unless defined $sarvar{"by_run"}->{"$sbatch_id"};
		$sarvar{"by_run"}->{"$sbatch_id"} = $sarvar{"by_run"}->{"$sbatch_id"}+1;
	}
	
	if (defined $vtw{"$varname"}) {
		$vtw{"$varname"}->{"by_county"} = {} unless defined $vtw{"$varname"}->{"by_county"};
		$vtw{"$varname"}->{"by_county"}->{"$county"} = $datecomp unless defined $vtw{"$varname"}->{"by_county"}->{"$county"};
		$vtw{"$varname"}->{"by_county"}->{"$county"} = $datecomp if $datecomp > $vtw{"$varname"}->{"by_county"}->{"$county"};
		$vtw{"$varname"}->{"by_sample"} = {} unless defined $vtw{"$varname"}->{"by_sample"};
		$vtw{"$varname"}->{"by_sample"}->{"$sample_id"} = 0 unless defined $vtw{"$varname"}->{"by_sample"}->{"$sample_id"};
		$vtw{"$varname"}->{"by_sample"}->{"$sample_id"} = $vtw{"$varname"}->{"by_sample"}->{"$sample_id"} + $prop;
	} elsif (defined $vtw{"$lineage"}) {
		$vtw{"$lineage"}->{"by_county"} = {} unless defined $vtw{"$lineage"}->{"by_county"};
		$vtw{"$lineage"}->{"by_county"}->{"$county"} = $datecomp unless defined $vtw{"$lineage"}->{"by_county"}->{"$county"};
		$vtw{"$lineage"}->{"by_county"}->{"$county"} = $datecomp if $datecomp > $vtw{"$lineage"}->{"by_county"}->{"$county"};
		$vtw{"$lineage"}->{"by_sample"} = {} unless defined $vtw{"$lineage"}->{"by_sample"};
		$vtw{"$lineage"}->{"by_sample"}->{"$sample_id"} = 0 unless defined $vtw{"$lineage"}->{"by_sample"}->{"$sample_id"};
		$vtw{"$lineage"}->{"by_sample"}->{"$sample_id"} = $vtw{"$lineage"}->{"by_sample"}->{"$sample_id"} + $prop;
	}
	
	$sarvar{"by_variant"}->{"$varname"} = 0 unless defined $sarvar{"by_variant"}->{"$varname"};
	$sarvar{"by_variant"}->{"$varname"} = $sarvar{"by_variant"}->{"$varname"}+1;

	$sarvar{"by_lineage"}->{"$lineage"} = 0 unless defined $sarvar{"by_lineage"}->{"$lineage"};
	$sarvar{"by_lineage"}->{"$lineage"} = $sarvar{"by_lineage"}->{"$lineage"}+1;

}
close $fh;

#print Dumper(\%vtw);
#die;

my %sdate = ("yr" => 0, "mo" => 0, "day" => 0);
if ($last_sample =~ m/(\d{4})(\d{2})(\d{2})/) {
	$sdate{"yr"}  = $1;
	$sdate{"mo"}  = $2;
	$sdate{"day"} = $3;
}

# Print summary info to STDOUT
print "#\n# Summary of WaTCH-WV seqr database.\n";
print "# Generated on $NOW.\n#\n";
print "# NOTE: WaTCH level runs from 1 (most concerning) to 5 (least concerning).\n";
print "# NOTE: It indicates the importance of a variant specifically for WV wastewater testing.\n";
print "# NOTE: It is informed by CDC variant classification levels (VOI, VOC, etc.).\n#\n";


print "\n\n########################################\n";
print "Overview:";
print "\n########################################\n\n";
print "Total samples sequenced     : " . scalar(keys %sample_ids) . ".\n";
print "Date last sample collected  : " . $sdate{"mo"} . "/" . $sdate{"day"} . "/" . $sdate{"yr"} . ".\n";
print "Date of last sequencing run : " . $last_run_h{"mo"} . "/" . $last_run_h{"day"} . "/" . $last_run_h{"yr"} . ".\n";
print "Total variants identified   : " . scalar(keys %{$sarvar{"by_variant"}}) . ".\n";
print "Total lineages identified   : " . scalar(keys %{$sarvar{"by_lineage"}}) . ".\n";


# sort by variant name
my @vtw_sorted = sort {$vtw{"$a"} cmp $vtw{"$b"}} keys %vtw;


print "\n\n########################################\n";
print "Distribution of highest concern variants (WaTCH level = 1) across all samples:";
print "\n########################################\n\n";
foreach my $lineage (@vtw_sorted) {
	next if $vtw{"$lineage"}->{"watch_status"} > 1;
	if (defined $vtw{"$lineage"}->{"by_sample"}) {
		print "$vtw{$lineage}->{label} variant $lineage found in " . scalar(keys %{$vtw{"$lineage"}->{"by_sample"}}) . " sample";
		if (scalar(keys %{$vtw{"$lineage"}->{"by_sample"}}) == 1) {
			print ":\n";
		} else {
			print "s:\n";
		}
		foreach my $sample_id (keys %{$vtw{"$lineage"}->{"by_sample"}}) {
			my $sample_ref = $sample_ids{"$sample_id"};
			print "   $sample_ref->{print_date} (run $sample_ref->{sbatch_id}) -> " . sprintf("%.1f", (100*$vtw{"$lineage"}->{"by_sample"}->{"$sample_id"})) . "% of sample $sample_id from $sample_ref->{location} in $sample_ref->{county} county.\n";
		}
		print "\n";
	}
}


print "\n\n########################################\n";
print "Geographic distribution of all variants to watch (WaTCH level <= $warn_level, sorted by most recent detection):";
print "\n########################################\n\n";
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

print "\n\n########################################\n";
print "Top 20 most abundant variants across all samples (any  WaTCH level):";
print "\n########################################\n\n";
my @variants_sorted = sort {$sarvar{"by_variant"}->{"$b"} <=> $sarvar{"by_variant"}->{"$a"}} keys %{$sarvar{"by_variant"}};
my $i = 0;
foreach my $variant (@variants_sorted) {
	$i++;
	print "$sarvar{by_variant}->{$variant} samples contain $variant\n";
	last if $i == 20;
}

print "\n";
print "\n# EOF.\n";


exit;

