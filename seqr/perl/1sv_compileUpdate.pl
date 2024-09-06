#! /usr/bin/env perl

use strict;
use warnings;

use DateTime qw( );
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=  "Using sequence demix proportion files located in RUNDIR, compile an update file ";
$usage   .=   "and write it to STDOUT.\n";
$usage   .=   "\n";

my $rundir;
my $status = 0;

# Uses file from this original source:
# https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json

my $SEQR_ALIASES = "resources/alias_key.json";


while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } else {
		$rundir = "$arg";
	}
}

die "FATAL: A run directory containing demix files must be provided!\n$usage\n" unless defined $rundir;
die "FATAL: Run directory $rundir is not a readable directory.\n$usage\n" unless -d $rundir;


my %aliases  = ();	# Stores info about variant aliases.
my %props    = ();	# Stores the actual variant proportions from the run in $rundir.

# Read in the variant alias key (json)
#
open my $fhA, "<", "$SEQR_ALIASES" or die "Unable to open $SEQR_ALIASES: $!\n";
while (my $line = <$fhA>) {
	chomp $line;
	next if "$line" =~ /^{/ or "$line" =~ /^}/;
	$line =~ s/ //gi;
	my ($primary, $astring) = split /:/, "$line", -1;

	$primary =~ s/"//g;
		
	$astring =~ s/,$//;
	$astring =~ s/"|\[|\]//g;
	
	unless ("$astring" eq "") {
		$aliases{"$primary"} = {} unless defined $aliases{"$primary"};
		foreach my $alias (split /,/, "$astring", -1) {
			$aliases{"$alias"} = {} unless defined $aliases{"$alias"};
			$aliases{"$alias"}->{"$primary"} = 1;
			$aliases{"$primary"}->{"$alias"} = 1;
		}
	}

}
close $fhA;


=cut
# Add the VTW info to the lineages hash, making sure they are de-aliased as needed!
#
foreach my $vtw_k (keys %vtw) {
	my $vtw_set = 0;
	if (defined $lineages{"$vtw_k"}) {
		$lineages{"$vtw_k"}->{"status"} = $vtw{"$vtw_k"}->{"status"};
		$lineages{"$vtw_k"}->{"WHO_label"} = $vtw{"$vtw_k"}->{"WHO_label"};
		$vtw_set = 1;
	} else {
		foreach my $lineage (keys %lineages) {
			if (defined $lineages{"$lineage"}->{"aliases"}->{"$vtw_k"}) {
				$lineages{"$lineage"}->{"status"} = $vtw{"$vtw_k"}->{"status"};
				$lineages{"$lineage"}->{"WHO_label"} = $vtw{"$vtw_k"}->{"WHO_label"};
				$vtw_set = 1;
				last;
			}
		}
	}
	if ($vtw_set == 0) {
		$lineages{"$vtw_k"} = {"aliases" => {}, 
													 "status" => $vtw{"$vtw_k"}->{"status"},
													 "WHO_label" => $vtw{"$vtw_k"}->{"WHO_label"}
													 };
	}
}
=cut


#print Dumper(\%aliases);
#die;



#
# Read in the variant proportion files from $rundir.
#
my @files = glob(qq("$rundir/*.txt"));
die "FATAL: $progname requires at least one txt file with freyja proportions in the run directory.\n$usage\n" unless scalar(@files) > 0;
#print Dumper(\@files);
#die;

foreach my $f (@files) {
	next if "$f" =~ /_NTC.txt/i;
	open my $fh, "<", "$f" or die "Unable to open $f: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		my ($runid, $sample_id, $varname, $prop) = split /\t/, "$line", -1;
		next unless defined $runid and $runid ne "";
#		my $lineage = makeLineage($varname);
#		my $parental = makeParental($varname);
		$props{"$sample_id"} = {} unless defined $props{"$sample_id"};
		$props{"$sample_id"}->{"$runid"} = {} unless defined $props{"$sample_id"}->{"$runid"};
		$props{"$sample_id"}->{"$runid"}->{"$varname"} = $prop;
	}
	close $fh;
}

#print Dumper(\%props);
#die;

print "sample_id\tsbatch_id\tvariant\tproportion\taliases\n";
foreach my $sample_id (keys %props) {
	foreach my $runid (keys %{$props{"$sample_id"}}) {
		foreach my $varname (keys %{$props{"$sample_id"}->{"$runid"}}) {
			my $astring = "$varname";
			$astring .= "," . join(",", keys(%{$aliases{"$varname"}})) if defined $aliases{"$varname"};
			
			print "$sample_id\t";
			print "$runid\t";
			print "$varname\t";
			print "$props{$sample_id}->{$runid}->{$varname}\t";
			print "$astring\n";
		}
	}
}

exit $status;



