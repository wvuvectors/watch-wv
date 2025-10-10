#! /usr/bin/env perl

use strict;
use warnings;

use DateTime qw( );
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options]\n";
$usage   .=  "Compiles a seqr update file from sequence demix proportion files.\n";
$usage   .=   "  [-i --indir]  Path to the folder that contains the demix files.\n";
$usage   .=   "  [-o --outdir] Path to the folder that will contain the update file(s).\n";
$usage   .=   "\n";

my ($indir, $outdir);
my $status = 0;

# Uses file from this original source:
# https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json
my $SEQR_ALIASES = "resources/alias_key.json";


while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } elsif ("$arg" eq "-i" or "$arg" eq "--indir") {
		defined ($indir = shift) or die "FATAL: Malformed argument to -i in $progname.\n$usage\n";
  } elsif ("$arg" eq "-o" or "$arg" eq "--outdir") {
		defined ($outdir = shift) or die "FATAL: Malformed argument to -o in $progname.\n$usage\n";
	}
}

die "FATAL: An input directory containing demix files must be provided!\n$usage\n" unless defined $indir;
die "FATAL: Input dir $indir is not a readable dir.\n$usage\n" unless -d "$indir";

die "FATAL: An output directory to hold the update file(s) must be provided!\n$usage\n" unless defined $outdir;
die "FATAL: Output dir $outdir is not a readable dir.\n$usage\n" unless -d "$outdir";



my %aliases  = ();	# Stores info about variant aliases.
my %props    = ();	# Stores the actual variant proportions from the run in $indir.

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
# Read in the variant proportion files from $indir.
#
my @files = glob(qq("$indir/*.txt"));
die "FATAL: $progname requires at least one txt file with freyja proportions in the input directory.\n$usage\n" unless scalar(@files) > 0;
#print Dumper(\@files);
#die;

foreach my $f (@files) {
	next if "$f" =~ /_NTC.txt/i;
	open my $infh, "<", "$f" or die "Unable to open $f for reading: $!\n";
	while (my $line = <$infh>) {
		chomp $line;
		my ($runid, $sample_id, $varname, $prop) = split /\t/, "$line", -1;
		next unless defined $runid and $runid ne "";
#		my $lineage = makeLineage($varname);
#		my $parental = makeParental($varname);
		$props{"$sample_id"} = {} unless defined $props{"$sample_id"};
		$props{"$sample_id"}->{"$runid"} = {} unless defined $props{"$sample_id"}->{"$runid"};
		$props{"$sample_id"}->{"$runid"}->{"$varname"} = $prop;
	}
	close $infh;
}

#print Dumper(\%props);
#die;

open my $outfh, ">", "$outdir/update.seqr.txt" or die "Unable to open $outdir/update.seqr.txt for writing: $!\n";

print $outfh "sample_id\tsbatch_id\tvariant\tproportion\taliases\n";
foreach my $sample_id (keys %props) {
	foreach my $runid (keys %{$props{"$sample_id"}}) {
		foreach my $varname (keys %{$props{"$sample_id"}->{"$runid"}}) {
			my $astring = "$varname";
			$astring .= "," . join(",", keys(%{$aliases{"$varname"}})) if defined $aliases{"$varname"};
			
			print $outfh "$sample_id\t";
			print $outfh "$runid\t";
			print $outfh "$varname\t";
			print $outfh "$props{$sample_id}->{$runid}->{$varname}\t";
			print $outfh "$astring\n";
		}
	}
}
close $outfh;

exit $status;



