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

my $SEQR_ALIASES = "resources/alias_key.json";
my $SEQR_VOIS    = "resources/vtws.txt";


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


my %vtw      = ();	# Stores info about variant priority (to watch; of concern; of interest).
my %aliases  = ();	# Stores info about variant aliases.
my %props    = ();	# Stores the actual variant proportions from the run in $rundir.

# Read in the file containing SARS variants to watch (VOCs, VBMs, VOIs, etc.)
#
open my $fhI, "<", "$SEQR_VOIS" or die "Unable to open $SEQR_VOIS: $!\n";
while (my $line = <$fhI>) {
	chomp $line;
	if ("$line" =~ /^## (.+?)$/) {
		warn "INFO : List of variants to watch was last updated on $1.\n";
	} elsif ("$line" =~ /^### (.+?)$/) {
		warn "INFO : List of variants to watch is derived from: $1.\n";
	} else {
		next if "$line" =~ /^WHO/i or "$line" =~ /^#/;
		my ($who_label, $pango_lineage, $cdc_status, $watch_status) = split /\t/, "$line", -1;
		$vtw{"$pango_lineage"} = {} unless defined $vtw{"$pango_lineage"};
		$vtw{"$pango_lineage"}->{"cdc_status"}   = "$cdc_status";
		$vtw{"$pango_lineage"}->{"watch_status"} = "$watch_status";
		$vtw{"$pango_lineage"}->{"who_label"}    = "$who_label";
	}
}
close $fhI;

#print Dumper(\%vtw);
#die;


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

print "sbatch_id\tsample_id\tvariant\tproportion\taliases\tcdc_status\twatch_status\n";
foreach my $sample_id (keys %props) {
	foreach my $runid (keys %{$props{"$sample_id"}}) {
		foreach my $varname (keys %{$props{"$sample_id"}->{"$runid"}}) {
			my $astring = "";
			$astring = join(",", keys(%{$aliases{"$varname"}})) if defined $aliases{"$varname"};
			
			my ($status_cdc, $status_watch, $who) = ("NA", "5", "NA");
			if (defined $vtw{"$varname"}) {
				$status_cdc   = "$vtw{$varname}->{cdc_status}";
				$status_watch = "$vtw{$varname}->{watch_status}";
				$who          = "$vtw{$varname}->{who_label}";
			} else {
				foreach my $alias (keys %{$aliases{"$varname"}}) {
					if (defined $vtw{"$alias"}) {
						$status_cdc   = "$vtw{$alias}->{cdc_status}";
						$status_watch = "$vtw{$alias}->{watch_status}";
						$who          = "$vtw{$alias}->{who_label}";
						last;
					}
				}
			}
			print "$runid\t";
			print "$sample_id\t";
			print "$varname\t";
			print "$props{$sample_id}->{$runid}->{$varname}\t";
			print "$astring\t$status_cdc\t$status_watch\n";
		}
	}
}

exit $status;



=cut
sub makeParental {
	my $variant = shift;

	if ("$variant" !~ /\./) {
		return "$variant";
	}

#	if (defined $vtw{"$variant"}) {
		#print "Found $variant (a Variant to Watch) in the new data!\n";
#		return "$variant";
#	}
	
	my @subs = split /\./, "$variant", -1;
	my $variant_base = "$variant";
	$variant_base = "$subs[0]" if scalar(@subs) > 0;
	
	if ("$variant_base" eq "XBB") {
		my $parental = "$variant_base";
		$parental = "$parental.$subs[1]" if scalar(@subs) > 1;
		$parental = "$parental.$subs[2]" if scalar(@subs) > 2;	
		return "$parental";
	}
	
	if (defined $lineages{"$variant_base"}) {

		my $parental = $lineages{"$variant_base"};

		if ("$parental" =~ /^B\.1\.1\.529/i) {
			return "B.1.1.529";
		}

		if ("$parental" =~ /^XBB/i) {
			my @asubs = split /\./, "$parental", -1;
			$parental = "XBB.$asubs[1]" if scalar(@asubs) > 1;
			$parental = "$parental.$asubs[2]" if scalar(@asubs) > 2;	
			return "$parental";
		}

		return "$parental";
	}
	
	my $parental = "$variant_base";
	$parental = "$parental.$subs[1]" if scalar(@subs) > 1;
	return "$parental";
}


sub makeLineage {
	my $variant = shift;
	
	if ("$variant" !~ /\./) {
		return "$variant";
	}
	
	foreach my $vartw (keys %vtw) {
		if ("$vartw" eq "$variant") {
			return "$variant";
		} elsif ("$vartw" =~ /\*$/) {
			my $base = "$vartw";
			#print "Trying to chop $base\n";
			chop($base);
			return "$vartw" if (defined $base and index("$variant", "$base") != -1);
		}
	}
	
#	if (defined $vtw{"$variant"}) {
		#print "Found $variant (a Variant to Watch) in the new data!\n";
#		return "$variant";
#	}
	
	my @subs = split /\./, "$variant", -1;
	my $lineage_base = "$subs[0]";

	my $lineage = "$lineage_base.$subs[1]" if scalar(@subs) > 1;
	$lineage = "$lineage.$subs[2]" if "$subs[0]" eq "XBB" and scalar(@subs) > 2;	
	return "$lineage";
	
}
=cut







