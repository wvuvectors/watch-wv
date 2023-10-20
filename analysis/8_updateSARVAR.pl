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
die "FATAL: $progname requires at least one txt file with freyja proportions in the run directory.\n$usage\n" unless scalar(@files) > 0;
#print Dumper(\@files);
#die;

my %metadata = ();
my %props    = ();
my %lineages = ();


my $WATCHDB_MAIN  = "updates/watchdb.LATEST.txt";
my $SARVARDB_MAIN = "updates/sarvardb.LATEST.txt";
my $SARVARDB_INCR = "updates/sarvardb/sarvardb.$NOW.txt";

my $SARVAR_ALIASES = "resources_seq/alias_key.json";
my $SARVAR_VOIS    = "resources_seq/VTW.txt";


# Read in the file containing SARS variants to watch (VOCs, VBMs, VOIs, etc.)
#
my %vtw = ();
print "Reading variants to watch from $SARVAR_VOIS.\n";
open my $fhI, "<", "$SARVAR_VOIS" or die "Unable to open $SARVAR_VOIS: $!\n";
while (my $line = <$fhI>) {
	chomp $line;
	if ("$line" =~ /^## (.+?)$/) {
		print "Note: list of variants to watch was last updated on $1.\n";
	} elsif ("$line" =~ /^### (.+?)$/) {
		print "Taken from: $1.\n";
	} else {
		next if "$line" =~ /^WHO/i or "$line" =~ /^#/;
		my ($who_label, $lineage, $status) = split /\t/, "$line", -1;
		$vtw{"$lineage"} = {} unless defined $vtw{"$lineage"};
		$vtw{"$lineage"}->{"status"} = "$status";
		$vtw{"$lineage"}->{"WHO_label"} = "$who_label";
	}
}
close $fhI;
print "Done.\n";

#print Dumper(\%vtw);
#die;

# Read in the SARS variant alias key (json)
#
print "Reading variant aliases from $SARVAR_ALIASES.\n";
open my $fhA, "<", "$SARVAR_ALIASES" or die "Unable to open $SARVAR_ALIASES: $!\n";
while (my $line = <$fhA>) {
	chomp $line;
	next if "$line" =~ /^{/ or "$line" =~ /^}/;
	$line =~ s/ //gi;
	my ($lineage, $aliases) = split /:/, "$line", -1;

	$lineage =~ s/"//g;
	
	$aliases =~ s/,$//;
	$aliases =~ s/"|\[|\]//g;

	unless ("$aliases" eq "") {
		my @alist = split /,/, "$aliases", -1;
		my %alist_unique = ();
		foreach my $a (@alist) {
			$alist_unique{"$a"} = 1;
		}
		my $astring = join ",", keys %alist_unique;
		$lineages{"$lineage"} = "$astring";
	}

}
close $fhA;
print "Done.\n";


=cut
# Add the VTW info to the lineages hash, making sure they are de-aliased as needed!
#
print "Adding info about variants to watch to the lineages hash.\n";
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

print "Done.\n";
=cut


#print Dumper(\%lineages);
#die;


# Read in the existing SARS variant database file
#
print "Reading the existing SARS variant database $SARVARDB_MAIN.\n";
open my $fhV, "<", "$SARVARDB_MAIN" or die "Unable to open $SARVARDB_MAIN: $!\n";
while (my $line = <$fhV>) {
	next if "$line" =~ /^sample_id/;
	chomp $line;
	my ($at, $runid, $facility, $county, $start, $end, $parental, $lineage, $varname, $prop) = split /\t/, "$line", -1;
	next unless defined $runid and $runid ne "";
	#my $parental = makeParental($varname);
	$props{"$at"} = {} unless defined $props{"$at"};
	$props{"$at"}->{"$runid"} = {} unless defined $props{"$at"}->{"$runid"};
	$props{"$at"}->{"$runid"}->{"$varname"} = {"proportion" => $prop, "parental" => "$parental", "lineage" => "$lineage"};
}
close $fhV;
print "Done.\n";



# Read in the variant proportion files from the rundir dashboard directory
#
print "Reading new variant proportion files from $rundir.\n";
foreach my $f (@files) {
	open my $fh, "<", "$f" or die "Unable to open $f: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		my ($runid, $at, $varname, $prop) = split /\t/, "$line", -1;
		next unless defined $runid and $runid ne "";
		my $lineage = makeLineage($varname);
		my $parental = makeParental($varname);
		$props{"$at"} = {} unless defined $props{"$at"};
		$props{"$at"}->{"$runid"} = {} unless defined $props{"$at"}->{"$runid"};
		$props{"$at"}->{"$runid"}->{"$varname"} = {"proportion" => $prop, "lineage" => "$lineage", "parental" => "$parental"};
		$metadata{"$at"} = {} unless defined $metadata{"$at"};
	}
	close $fh;
}
print "Done.\n";

#print Dumper(\%props);
#die;



# Read in the existing WaTCH database file
# Restricted to those samples in either the existing sarvarDB or the new run
#
print "Reading the existing WaTCH database from $WATCHDB_MAIN.\n";
print "This allows us to link sample IDs from the variant data to metadata such as county, dates, etc.\n";
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
print "Done.\n";


#print Dumper(\%metadata);
#die;


print "Everything seems to have been read ok, so preparing to write updated sarsvarDB now.\n";
print "Backing up $SARVARDB_MAIN to $SARVARDB_MAIN.OLD.\n";
`cp $SARVARDB_MAIN $SARVARDB_MAIN.OLD`;

print "Printing updated database to $SARVARDB_INCR.\n";
open my $fho, ">", "$SARVARDB_INCR" or die "Unable to open $SARVARDB_INCR for writing: $!\n";
print $fho "sample_id\tseqrun_id\tfacility\tcounty\tstart_datetime\tend_datetime\tparental\tlineage\tvariant\tproportion\n";
foreach my $at (keys %props) {
	foreach my $rid (keys %{$props{"$at"}}) {
		foreach my $varname (keys %{$props{"$at"}->{"$rid"}}) {
			print $fho "$at\t";
			print $fho "$rid\t";
			print $fho "$metadata{$at}->{facility}\t";
			print $fho "$metadata{$at}->{county}\t";
			print $fho "$metadata{$at}->{start_datetime}\t";
			print $fho "$metadata{$at}->{end_datetime}\t";
			print $fho "$props{$at}->{$rid}->{$varname}->{parental}\t";
			print $fho "$props{$at}->{$rid}->{$varname}->{lineage}\t";
			print $fho "$varname\t";
			print $fho "$props{$at}->{$rid}->{$varname}->{proportion}\n";
		}
	}
}
close $fho;

print "Copying updated database to $SARVARDB_MAIN.\n";
`cp $SARVARDB_INCR $SARVARDB_MAIN`;

print "Done!\n";
exit;




sub makeParental {
	my $variant = shift;
#	if ("$variant" !~ /\./) {
#		return "$variant";
#	}

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

	if (defined $vtw{"$variant"}) {
		#print "Found $variant (a Variant to Watch) in the new data!\n";
		return "$variant";
	}
	
	my @subs = split /\./, "$variant", -1;
	my $lineage_base = "$subs[0]";

	my $lineage = "$lineage_base.$subs[1]" if scalar(@subs) > 1;
	$lineage = "$lineage.$subs[2]" if "$subs[0]" eq "XBB" and scalar(@subs) > 2;	
	return "$lineage";
	
}







