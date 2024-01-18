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
$usage   .=   "  [-i --infile]    Path to the update table.\n";
$usage   .=   "  [-o --outdir]    Path to the output directory.\n";
$usage   .=   "  [-w --watchdb]   Path to the watch database directory to inform the update (default: ../patchr/data/latest).\n";
$usage   .=   "  [-r --relineage] Recalculate the parental and lineage assignments for all entries.\n";
$usage   .=   "\n";

my ($infile, $outdir);
my $dbdir = "../patchr/data/latest";
my $relin = 0;

my $status = 0;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } elsif ("$arg" eq "-i" or "$arg" eq "--infile") {
		defined ($infile = shift) or die "FATAL: Malformed argument to -i in $progname.\n$usage\n";
  } elsif ("$arg" eq "-o" or "$arg" eq "--outdir") {
		defined ($outdir = shift) or die "FATAL: Malformed argument to -o in $progname.\n$usage\n";
  } elsif ("$arg" eq "-w" or "$arg" eq "--watchdb") {
		defined ($dbdir = shift) or die "FATAL: Malformed argument to -w in $progname.\n$usage\n";
  } elsif ("$arg" eq "-r" or "$arg" eq "--relineage") {
		$relin = 1;
	}
}

die "FATAL: An update file containing demix data must be provided (-i).\n$usage\n" unless defined $infile;
die "FATAL: Update file $infile is not a readable file.\n$usage\n" unless -f $infile;

die "FATAL: An output directory for the seqrdb tables must be provided (-o).\n$usage\n" unless defined $outdir;
die "FATAL: Output directory $outdir is not a readable directory.\n$usage\n" unless -d $outdir;

die "FATAL: watchDB directory $dbdir is not a readable directory.\n$usage\n" unless -d $dbdir;

my %metadata = ();
my %props    = ();


# Read in the latest seqr database file to populate the props hash.
# sample_id
# sbatch_id
# location
# county
# sample_collection_start_datetime
# sample_collection_end_datetime
# lineage_group
# variant
# variant_proportion
#
if (-f "data/latest/seqrdb.txt") {
	open my $fhV, "<", "data/latest/seqrdb.txt" or die "Unable to open data/latest/seqrdb.txt for reading: $!\n";
	while (my $line = <$fhV>) {
		next if "$line" =~ /^sample_id/;
		chomp $line;

		my @cols = split /\t/, "$line", -1;

		# Recalcuate the lineage group for previous variants, if that flag is set by the user.
		$cols[6] = makeLineageGroup("$cols[7]") if $relin == 1;
		
		# Add the data to the props hash, keyed on sample id, run id, and variant.
		$props{"$cols[0]"} = {} unless defined $props{"$cols[0]"};
		$props{"$cols[0]"}->{"$cols[1]"} = {} unless defined $props{"$cols[0]"}->{"$cols[1]"};
		$props{"$cols[0]"}->{"$cols[1]"}->{"$cols[7]"} = {
			"proportion" => $cols[8], 
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
}



# Read in the seqrdb update file.
# sbatch_id
# sample_id
# variant
# proportion
# aliases
# cdc_status
# watch_status
#
my $linenum = 0;
open my $fhI, "<", "$infile" or die "Unable to open $infile: $!\n";
while (my $line = <$fhI>) {
	chomp $line;
	
	my @cols = split /\t/, "$line", -1;
	
	if ($linenum == 0) {
		$linenum++;
	} else {
		$props{"$cols[1]"} = {} unless defined $props{"$cols[1]"};
		$props{"$cols[1]"}->{"$cols[0]"} = {} unless defined $props{"$cols[1]"}->{"$cols[0]"};
		$props{"$cols[1]"}->{"$cols[0]"}->{"$cols[2]"} = {
			"proportion"    => $cols[3],
			"aliases"       => "$cols[4]",
			"cdc_status"    => "$cols[5]",
			"watch_status"  => $cols[6],
			"lineage_group" => ""
		};

		if ($cols[6] < 3) {
			$props{"$cols[1]"}->{"$cols[0]"}->{"$cols[2]"}->{"lineage_group"} = "$cols[2]";
		} else {
			$props{"$cols[1]"}->{"$cols[0]"}->{"$cols[2]"}->{"lineage_group"} = makeLineageGroup("$cols[2]");
		}
		
		# Add the sample id to the metadata hash for lookup, if it doesn't already exist.
		$metadata{"$cols[1]"} = {} unless defined $metadata{"$cols[1]"};
	}
}
close $fhI;

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
	open (my $tFH, "<", "$dbdir/watchdb.$table.txt") or die "Unable to open $dbdir/watchdb.$table.txt for reading: $!\n";
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


# Write the updated seqrdb file to $outdir.
#
open my $fho, ">", "$outdir/seqrdb.txt" or die "Unable to open $outdir/seqrdb.txt for writing: $!\n";
print $fho "sample_id\tsbatch_id\tlocation\tcounty\tsample_collection_start_datetime\tsample_collection_end_datetime\tlineage_group\tvariant\tvariant_proportion\n";
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
	my $variant = shift;

	if ("$variant" !~ /\./) {
		return "$variant";
	}
	
	my $group = "";
	
	my @subs = split /\./, "$variant", -1;
	my $base = "$variant";
	$base = "$subs[0]" if scalar(@subs) > 0;
	
	if ("$base" eq "XBB") {
		$group = "$base";
		$group = "$group.$subs[1]" if scalar(@subs) > 1;
		$group = "$group.$subs[2]" if scalar(@subs) > 2;	
	} elsif ("$base" =~ /^B\.1\.1\.529/i) {
		$group = "B.1.1.529";
	} else {
		my $group = "$base";
		$group = "$group.$subs[1]" if scalar(@subs) > 1;
	}
	
	return "$group";
}




