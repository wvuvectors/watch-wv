#! /usr/bin/env perl

use strict;
use warnings;

use DateTime qw( );
use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options]\n";
$usage   .=  "Using the sequence demix proportion file (-i), update the seqr database in the given outdir (-o).\n";
$usage   .=   "  [-d --dbdir]    Path to the folder that contains the full seqr database file(s).\n";
$usage   .=   "  [-u --update]    Path to the folder that contains the update file(s).\n";
$usage   .=   "  [-o --outdir]    Path to the output directory.\n";
$usage   .=   "  [-w --watchdb]   Path to the watch database directory to inform the update (default: ../patchr/data/latest).\n";
$usage   .=   "  [-r --relineage] Recalculate the parental and lineage assignments for all entries.\n";
$usage   .=   "\n";

my ($dbdir, $updir, $outdir);
my $watchdir = "../patchr/data/latest";
my $relin = 0;

my %aliases  = ();	# Stores info about variant aliases.
my %vtw      = ();	# Stores info about variants of concern.
my %metadata = ();	# Stores sample metadata from watchdb.
my %props    = ();	# Stores variant proportion data.

my $status = 0;

my $SEQR_ALIASES = "resources/alias_key.json";
my $SEQR_VOIS    = "resources/vtws.txt";


while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } elsif ("$arg" eq "-d" or "$arg" eq "--dbdir") {
		defined ($dbdir = shift) or die "FATAL: Malformed argument to -d in $progname.\n$usage\n";
  } elsif ("$arg" eq "-u" or "$arg" eq "--update") {
		defined ($updir = shift) or die "FATAL: Malformed argument to -u in $progname.\n$usage\n";
  } elsif ("$arg" eq "-o" or "$arg" eq "--outdir") {
		defined ($outdir = shift) or die "FATAL: Malformed argument to -o in $progname.\n$usage\n";
  } elsif ("$arg" eq "-w" or "$arg" eq "--watchdb") {
		defined ($watchdir = shift) or die "FATAL: Malformed argument to -w in $progname.\n$usage\n";
  } elsif ("$arg" eq "-r" or "$arg" eq "--relineage") {
		$relin = 1;
	}
}

die "FATAL: A directory containing the FULL seqr datbase must be provided (-d).\n$usage\n" unless defined $dbdir;
die "FATAL: Seqr database dir $dbdir is not a readable dir.\n$usage\n" unless -d "$dbdir";
die "FATAL: Seqr database file $dbdir/watchdb.seqr.txt is not a readable file.\n$usage\n" unless -f "$dbdir/watchdb.seqr.txt";

die "FATAL: A directory containing the seqr datbase UPDATE must be provided (-u).\n$usage\n" unless defined $updir;
die "FATAL: Update dir $updir is not a readable dir.\n$usage\n" unless -d "$updir";
die "FATAL: Update file $updir/update.seqr.txt is not a readable file.\n$usage\n" unless -f "$updir/update.seqr.txt";

die "FATAL: An output directory for the new db version must be provided (-o).\n$usage\n" unless defined $outdir;
die "FATAL: Output directory $outdir is not a readable directory.\n$usage\n" unless -d $outdir;


die "FATAL: watchDB directory $watchdir is not a readable directory.\n$usage\n" unless -d "$watchdir";




# Read in the variant alias key (json) in case we need it.
# Only necessary if the $relin flag is set, since the update file comes in with aliases already.
#
# This is a 1-to-1 map between variant alias pairs.
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



# Read in the latest seqr file to populate the props hash.
# sample_id
# sbatch_id
# location
# county
# sample_collection_start_datetime
# sample_collection_end_datetime
# lineage_group
# variant
# variant_aliases
# variant_proportion
#
open my $fhV, "<", "$dbdir/watchdb.seqr.txt" or die "Unable to open $dbdir/watchdb.seqr.txt for reading: $!\n";
while (my $line = <$fhV>) {
	next if "$line" =~ /^sample_id/;
	chomp $line;

	my @cols = split /\t/, "$line", -1;

	# Recalcuate the lineage group for previous variants, if that flag is set by the user.
	if ($relin == 1) {
		my $alist = makeAliases("$cols[7]");
		$cols[8] = "$alist";
		$cols[6] = makeLineageGroup("$alist");
	}
	
	# Add the data to the props hash, keyed on sample id, run id, and variant.
	$props{"$cols[0]"} = {} unless defined $props{"$cols[0]"};
	$props{"$cols[0]"}->{"$cols[1]"} = {} unless defined $props{"$cols[0]"}->{"$cols[1]"};
	$props{"$cols[0]"}->{"$cols[1]"}->{"$cols[7]"} = {
		"aliases"       => "$cols[8]", 
		"proportion"    => $cols[9], 
		"lineage_group" => "$cols[6]"
	};
	
	# Add the sample metadata to the metadata hash, to minimize redundancy.
	$metadata{"$cols[0]"} = {
		"sample_id" => "$cols[0]",
		"location"  => "$cols[2]",
		"county"    => "$cols[3]",
		"sample_collection_start_datetime" => "$cols[4]",
		"sample_collection_end_datetime"   => "$cols[5]"
	} unless defined $metadata{"$cols[0]"};
}
close $fhV;



# Read in the seqr update file.
# sample_id
# sbatch_id
# variant
# proportion
# aliases
#
my $linenum = 0;
open my $fhU, "<", "$updir/update.seqr.txt" or die "Unable to open $updir/update.seqr.txt for reading: $!\n";
while (my $line = <$fhU>) {
	chomp $line;
	
	my @cols = split /\t/, "$line", -1;
	
	if ($linenum == 0) {
		$linenum++;
	} else {
		$props{"$cols[0]"} = {} unless defined $props{"$cols[0]"};
		$props{"$cols[0]"}->{"$cols[1]"} = {} unless defined $props{"$cols[0]"}->{"$cols[1]"};
		$props{"$cols[0]"}->{"$cols[1]"}->{"$cols[2]"} = {
			"proportion"    => $cols[3],
			"aliases"       => "$cols[4]",
			"lineage_group" => ""
		};
		
		$props{"$cols[0]"}->{"$cols[1]"}->{"$cols[2]"}->{"lineage_group"} = makeLineageGroup("$cols[4]");
		
		# Add the sample id to the metadata hash for lookup, if it doesn't already exist.
		$metadata{"$cols[0]"} = {} unless defined $metadata{"$cols[0]"};
	}
}
close $fhU;

#print Dumper(\%props);
#die;

#print Dumper(\%metadata);
#die;



# Read in the existing WaTCH database to populate (or update) sample metadata.
# We only need the sample and location tables for now.
#
my %tables = (
	"sample" => {}
);

my %resources = (
	"location" => {},
	"wwtp"		 => {},
	"county"	 => {},
	"lab"			 => {}
);

foreach my $table (keys %tables) {
	$linenum = 0;
	my @headers = ();
	open (my $tFH, "<", "$watchdir/watchdb.$table.txt") or die "Unable to open $watchdir/watchdb.$table.txt for reading: $!\n";
	while (my $line = <$tFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @cols = split "\t", "$line", -1;

		if ($linenum == 0) {
			$linenum++;
			for (@cols) {
				push @headers, trim("$_");
			}
		} else {
			# first field is always the uid regardless of table
			my $uid = trim("$cols[0]");
			
			# skip samples that aren't in the metadata hash, to save time.
			next unless defined $metadata{"$uid"};
			
			$tables{"$table"}->{"$uid"} = {};
			for (my $i=0; $i < scalar(@cols); $i++) {
				$tables{"$table"}->{"$uid"}->{"$headers[$i]"} = trim("$cols[$i]");
			}
		}
	}
	close $tFH;
}

#print Dumper(\%tables);
#die;


#
# Read WaTCHdb resources from resources/ dir.
#
# Read resource tables into hash
my $resource_wkbk = ReadData("../patchr/resources/watchdb.all_tables.xlsx", dtfmt => "mm/dd/yy");
#print Dumper($resource_wkbk);
#die;

foreach my $sheet_name (keys %{$resource_wkbk->[0]->{"sheet"}}) {
#	print "$sheet_name\n";
	my @field_names = ();

	$resources{"$sheet_name"} = {} unless defined $resources{"$sheet_name"};
	my $pos = $resource_wkbk->[0]->{"sheet"}->{"$sheet_name"};
	my $sheetRef = $resource_wkbk->[$pos];
	
	# The field labels are in row A (or, index 1 of each array in 'cell' ). Store these in an array.
	for (my $i=1; $i < scalar(@{$sheetRef->{"cell"}}); $i++) {
		#print "$sheetRef->{cell}->[$i]->[1]\n";
		next unless defined $sheetRef->{"cell"}->[$i]->[1] and "$sheetRef->{cell}->[$i]->[1]" ne "";
		my $col_name = trim("$sheetRef->{cell}->[$i]->[1]");
		push @field_names, "$col_name";
	}
	#print Dumper(\@field_names);
	#last;

	# The keys for each resource sheet are in column A (or, array 1 of 'cell'). 
	# Loop over these to add entries to the resources hash.
	for (my $i=2; $i < scalar(@{$sheetRef->{"cell"}->[1]}); $i++) {
		#print "$sheetRef->{cell}->[1]->[$i]\n";
		next unless defined $sheetRef->{"cell"}->[1]->[$i] and "$sheetRef->{cell}->[1]->[$i]" ne "";
		my $key = trim("$sheetRef->{cell}->[1]->[$i]");
		$resources{"$sheet_name"}->{"$key"} = {};
	}
	#last;
	
	# To get the actual values, loop over the 'cell' array.
	# Get the key from array 1.
	# Loop over the field_names array to get the values.
	for (my $i=0; $i < scalar(@{$sheetRef->{"cell"}})-1; $i++) {
		next unless defined $sheetRef->{"cell"}->[$i+1] and defined $field_names[$i] and "$field_names[$i]" ne "";
		my $colRef = $sheetRef->{"cell"}->[$i+1];
		my $field = trim("$field_names[$i]");
		my $resourceRef = $resources{"$sheet_name"};
		# Loop over the values in this column. Add each to the resources hash.
		for (my $j=2; $j < scalar(@{$colRef}); $j++) {
			next unless defined $sheetRef->{"cell"}->[1]->[$j];
			my $key = trim("$sheetRef->{cell}->[1]->[$j]");
			my $value = "";
			$value = trim("$colRef->[$j]") if defined $colRef->[$j];
			$resources{"$sheet_name"}->{"$key"}->{"$field"} = "$value";
		}
	}
}

#print Dumper(\%resources);
#die;



# Update the metadata hash with new or changed sample metadata from the watchdb.
#
foreach my $sample_id (keys %metadata) {
	next unless defined $tables{"sample"}->{"$sample_id"};
	
	my $sampleRef = $tables{"sample"}->{"$sample_id"};
	
	my $locid  = "";
	$locid = $sampleRef->{"location_id"} if defined $sampleRef->{"location_id"};
	my $start  = "";
	$start = $sampleRef->{"sample_collection_start_datetime"} if defined $sampleRef->{"sample_collection_start_datetime"};
	my $end    = "";
	$end = $sampleRef->{"sample_collection_end_datetime"} if defined $sampleRef->{"sample_collection_end_datetime"};
	my $county = "";
	$county = $resources{"location"}->{"$locid"}->{"location_counties_served"} if defined $resources{"location"}->{"$locid"}->{"location_counties_served"};
	
	$metadata{"$sample_id"}->{"location"} = "$locid";
	$metadata{"$sample_id"}->{"county"}   = "$county";
	$metadata{"$sample_id"}->{"sample_collection_start_datetime"} = "$start";
	$metadata{"$sample_id"}->{"sample_collection_end_datetime"}   = "$end";
}



# Write the updated seqr file to $outdir.
#
open my $fho, ">", "$outdir/watchdb.seqr.txt" or die "Unable to open $outdir/watchdb.seqr.txt for writing: $!\n";
print $fho "sample_id\tsbatch_id\tlocation\tcounty\tsample_collection_start_datetime\tsample_collection_end_datetime\tlineage_group\tvariant\tvariant_aliases\tvariant_proportion\n";
foreach my $sample_id (keys %props) {
	foreach my $sbid (keys %{$props{"$sample_id"}}) {
		foreach my $varname (keys %{$props{"$sample_id"}->{"$sbid"}}) {
			print $fho "$sample_id\t";
			print $fho "$sbid\t";
			print $fho "$metadata{$sample_id}->{location}\t";
			print $fho "$metadata{$sample_id}->{county}\t";
			print $fho "$metadata{$sample_id}->{sample_collection_start_datetime}\t";
			print $fho "$metadata{$sample_id}->{sample_collection_end_datetime}\t";
			print $fho "$props{$sample_id}->{$sbid}->{$varname}->{lineage_group}\t";
			print $fho "$varname\t";
			print $fho "$props{$sample_id}->{$sbid}->{$varname}->{aliases}\t";
			print $fho "$props{$sample_id}->{$sbid}->{$varname}->{proportion}\n";
		}
	}
}
close $fho;


exit $status;



sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}



sub makeLineageGroup {
	my $aliases_csv = shift;

	# input contains a comma-sep list of all aliases for this variant (including the "primary" name).
	# Use @varlist to go through each alias.
	#
	my @varlist = split /,/, "$aliases_csv", -1;

	my @voclist = keys %vtw;

	# Check each alias against our list of variants of concern (keys of the %vtw hash).
	# If found, use the pango_lineage as the LG.
	# Variants in %vtw may include regex terms so it requires a N^n approach.
	# Fortunately both arrays are short!
	#
	foreach my $v (@varlist) {
		foreach my $voc (@voclist) {
			if ("$voc" =~ /\*$/) {
				$voc =~ s/\*$//;
				return "$voc" if "$v" =~ /^$voc/;
			} else {
				return "$voc" if "$v" eq "$voc";
			}
		}
	}
	
	# If we get here, it isn't a known VOC.
	#
	foreach my $v (@varlist) {

	# If it's an XB+ lineage, use XB+.[0-9].[0-9]
#	if ("$v" =~ /(^XB.+?\.\d+\.\d+).*/) {
	# If it's an XB+ lineage, use XB+.[0-9]
	if ("$v" =~ /(^XB.+?\.\d+)\..*/) {
		return "$1";
	}
	
	# If it's a CR lineage, use CR.[0-9]
	if ("$v" =~ /(^CR\.\d+).*/) {
		return "$1";
	}
	
	# If it's a BA lineage, use BA.[0-9]
	if ("$v" =~ /(^BA\.\d+).*/) {
		return "$1";
	}
	
	# If it's a BQ lineage, use BQ.[0-9].[0-9]
	if ("$v" =~ /(^BQ\.\d+\.\d+).*/) {
		return "$1";
	}
	
	# If it's a B lineage, use B.[0-9].[0-9]
	if ("$v" =~ /(^B\.\d+\.\d+).*/) {
		return "$1";
	}
	
}

	# If we get here, it's an oddball.
	# Let's use the shortest aliases.
	#
	my $lg = "";
	foreach my $v (@varlist) {
		$lg = "$v" if length "$lg" == 0 or length "$v" < length "$lg";
	}
	return "$lg";
	
}



sub makeAliases {
	my $varname = shift;

	my $astring = "$varname";
	$astring .= "," . join(",", keys(%{$aliases{"$varname"}})) if defined $aliases{"$varname"};
	
	return "$astring";
}





