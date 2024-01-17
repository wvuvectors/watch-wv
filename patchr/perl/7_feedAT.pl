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

my $feed         = "AT";
my @feed_headers = ();
my %feed_cvmap   = ();

# Get the required feed fields (column headers), in order.
# These are located in the resources/ folder.
open (my $fFH, "<", "resources/${feed}_fields.txt") or die "FATAL: Unable to open resources/${feed}_fields.txt for reading: $!.";
while (my $line = <$fFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	push @feed_headers, "$line";
}
close $fFH;


# Get the controlled vocabulary map for this feed, if one exists.
# This is located in the resources/ folder.
my $linenum = 0;
if (-f "resources/${feed}_cvm.txt") {
	open (my $vFH, "<", "resources/${feed}_cvm.txt") or die "FATAL: Unable to open resources/${feed}_cvm.txt for reading: $!.";
	while (my $line = <$vFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		#my @cols = split /\t/, "$line", -1;
		my ($extref_field, $extref_value, $watch_table, $watch_field, $watch_value, $reference) = split /\t/, "$line", -1;
		if ($linenum == 0) {
			$linenum++;
		} else {
			$feed_cvmap{"$extref_field"} = {} unless defined $feed_cvmap{"$extref_field"};
			$feed_cvmap{"$extref_field"}->{"$watch_field"} = {} unless defined $feed_cvmap{"$extref_field"}->{"$watch_field"};
			$feed_cvmap{"$extref_field"}->{"$watch_field"}->{"$watch_value"} = {"value" => "$extref_value", "reference" => "$reference"};
		}
	}
	close $vFH;
}


my %tables = (
	"sample"   => {}
);

my %resources = (
	"location" => {}
);


#
# Read WaTCH database tables from dbdir.
#

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
				push @headers, "$_";
			}
		} else {
			# first field is always the uid regardless of table
			my $uid = $cols[0];

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
my $resource_wkbk = ReadData("resources/watchdb.all_tables.xlsx", dtfmt => "mm/dd/yy");
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
		my $col_name = "$sheetRef->{cell}->[$i]->[1]";
		push @field_names, "$col_name";
	}
	#print Dumper(\@field_names);
	#last;

	# The keys for each resource sheet are in column A (or, array 1 of 'cell'). 
	# Loop over these to add entries to the resources hash.
	for (my $i=2; $i < scalar(@{$sheetRef->{"cell"}->[1]}); $i++) {
		#print "$sheetRef->{cell}->[1]->[$i]\n";
		next unless defined $sheetRef->{"cell"}->[1]->[$i] and "$sheetRef->{cell}->[1]->[$i]" ne "";
		my $key = "$sheetRef->{cell}->[1]->[$i]";
		$resources{"$sheet_name"}->{"$key"} = {};
	}
	#last;
	
	# To get the actual values, loop over the 'cell' array.
	# Get the key from array 1.
	# Loop over the field_names array to get the values.
	for (my $i=0; $i < scalar(@{$sheetRef->{"cell"}})-1; $i++) {
		next unless defined $sheetRef->{"cell"}->[$i+1] and defined $field_names[$i] and "$field_names[$i]" ne "";
		my $colRef = $sheetRef->{"cell"}->[$i+1];
		my $field = "$field_names[$i]";
		my $resourceRef = $resources{"$sheet_name"};
		# Loop over the values in this column. Add each to the resources hash.
		for (my $j=2; $j < scalar(@{$colRef}); $j++) {
			next unless defined $sheetRef->{"cell"}->[1]->[$j];
			my $key = "$sheetRef->{cell}->[1]->[$j]";
			my $value = "";
			$value = "$colRef->[$j]" if defined $colRef->[$j];
			$resources{"$sheet_name"}->{"$key"}->{"$field"} = "$value";
		}
	}
}

#print Dumper(\%resources);
#die;


#
# Read data from WaTCH-AT mapping table. Keyed by AT column header.
#
my %map = ();

my $map_wkbk = ReadData("resources/${feed}_map.xlsx", dtfmt => "mm/dd/yy");
my @map_rows = Spreadsheet::Read::rows($map_wkbk->[1]);

my %map_headers = ();
for (my $i=0; $i < scalar(@{$map_rows[0]}); $i++) {
	$map_headers{$i} = $map_rows[0][$i];
}

for (my $i=1; $i < scalar(@map_rows); $i++) {
	my $extref_field = $map_rows[$i][1];
	$map{$extref_field} = {};
	for (my $j=0; $j < scalar(@{$map_rows[$i]}); $j++) {
		$map{$extref_field}->{"$map_headers{$j}"} = $map_rows[$i][$j];
	}
}
#print Dumper(\%map);
#die;

#print Dumper(\@feed_headers);
#die;


# To assemble the output file, we loop over entries in the sample table.
# For each sample, we loop over @feed_headers in order, assemble the line, and print it.

#open(my $oFH, ">", "$OUTFILE_BK") or die "Unable to open $OUTFILE_BK for writing: $!";
#print $oFH "\"" . join("\",\"", @feed_headers) . "\"\n";
print "\"" . join("\",\"", @feed_headers) . "\"\n";

RSLTLOOP:
foreach my $sample_id (keys %{$tables{"sample"}}) {
	my $resultRef = $tables{"sample"}->{"$sample_id"};
	
	# Skip unless sample_collection_datetime >= TODAY - days
	my $dt_string = $resultRef->{"sample_collection_start_datetime"};	# MM-DD-YYYY HH:MM:SS
	next if !defined $dt_string or "$dt_string" eq "";

	my ($d_str, $t_str) = split /\s/, "$resultRef->{sample_collection_start_datetime}", -1;
	my ($mo, $dd, $yyyy) = split /\//, "$d_str", -1;
	#my ($hh, $mi, $ss) = split /:/, "$t_str", -1;

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
				
	next if DateTime->compare($dt, $THRESHOLD_DATE) == -1;
	
	# Standardize samples with empty flow.
	$resultRef->{"sample_flow"} = "" if "$resultRef->{sample_flow}" eq "NA" or "$resultRef->{sample_flow}" eq "0";

	my @result_vals = ();
	my $val_pass    = 1;
	HDRLOOP:
	foreach my $header (@feed_headers) {
		# These are also the keys of our %map	hash.
		my %mapEntry = %{$map{"$header"}};
		my $req_ne = "yes";
		$req_ne = "no" unless lc($mapEntry{"required_nonempty"}) eq "yes";
		
		#print Dumper(\%mapEntry);
		#die;
		
		my $val = "";
		#print "$mapEntry{WaTCH_field}\n";
		
		if ("$mapEntry{WaTCH_field}" eq "[empty]") {
			# Print the empty string.
			$val = "";
		} elsif ("$mapEntry{WaTCH_field}" eq "[constant]") {
			# Print the value of the constant field.
			$val = "$mapEntry{constant_value}";
		} elsif ("$mapEntry{WaTCH_field}" eq "[indirect]") {
			# Look up value in a linked table.
			my $table = "$mapEntry{WaTCH_table}";
			$val = lookupLinked("$sample_id", "$table", "$mapEntry{linked_field}");
		} elsif ("$mapEntry{WaTCH_field}" eq "[calculated]") {
			# Call the calcExtref sub to figure this out.
			my $lfield = "";
			$lfield = "$mapEntry{linked_field}" if defined $mapEntry{"linked_field"};
			$val = calcExtref("$sample_id", "$header", "$lfield");
		} elsif ("$mapEntry{WaTCH_table}" eq "sample") {
			# Look up this field in the results table directly.
			my $table = $mapEntry{"WaTCH_table"};
			my $field = $mapEntry{"WaTCH_field"};
			$val = $resultRef->{"$field"} if defined $resultRef->{"$field"};
		}
		
		#last HDRLOOP if "$val" eq "" and "$req_ne" eq "yes";
		
		# convert date or time to correct format
		if ("$mapEntry{is_date_or_time}" eq "yes") {
			my $dtRef = formatDT("$val");
			$val = $dtRef->{"date"};
			$val = $dtRef->{"time"} if "$header" =~ /_time$/;
		}
		
		push @result_vals, "$val";
#		$val_pass = 1 unless "$val" eq "";
	}
#	print $oFH "\"" . join("\",\"", @result_vals) . "\"\n";
	print "\"" . join("\",\"", @result_vals) . "\"\n" unless $val_pass == 0;
}


exit $status;





sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


sub lookupLinked {
	my $sid    = shift;
	my $tname  = shift;
	my $tfield = shift;
	
	my $table_cat = "table";
	my $linkedId;

	if ("$tname" eq "location") {
		$table_cat = "resource";
		$linkedId = $tables{"sample"}->{"$sid"}->{"location_id"};
	}
	
	my $linked_val = "";
	if ("$table_cat" eq "table") {
		$linked_val = $tables{"$tname"}->{"$linkedId"}->{"$tfield"} if defined $tables{"$tname"}->{"$linkedId"}->{"$tfield"};
	} elsif ("$table_cat" eq "resource") {
		if ("$linkedId" =~ /,/) {
			my @a = split /,\s*/, "$linkedId", -1;
			my @b = ();
			foreach (@a) {
				push(@b, $resources{"$tname"}->{"$_"}->{"$tfield"}) if defined $resources{"$tname"}->{"$_"}->{"$tfield"};
			}
			$linked_val = join(",", @b);
		} else {
			$linked_val = $resources{"$tname"}->{"$linkedId"}->{"$tfield"} if defined $resources{"$tname"}->{"$linkedId"}->{"$tfield"};
		}
	}
	
	#print "$aid\t$tname\t$tfield\t$linked_val\n";
	return "$linked_val";
}


sub calcExtref {
	my $sid          = shift;
	my $extref_field = shift;
	my $linked_field = shift;
	
	my $val = "";
	
	if ("$extref_field" eq "Sample Collection Method") {
		# location.location_category but 'not institution specific' if wwtp
		my $locid = $tables{"sample"}->{"$sid"}->{"location_id"};
		$val = $resources{"location"}->{"$locid"}->{"location_collection_basis"} if defined $resources{"location"}->{"$locid"}->{"location_collection_basis"};
#		if (!defined $val) {
#			print "$sid\t$locid\n";
#			die;
#		}
		if (lc($val) eq "grab") {
			$val = "Grab";
		} else {
			$val = "Composite";
		}
	}
	
	return "$val";
}


sub formatDT {
	# Format the date and time appropriately. I hate date hacking.
	# DateTime::Format::Strptime doesn't work although I have not given up hope.
	#
	my $dtstring = shift;
	
	my %dtHash = ("date" => "", "time" => "");
	my ($mon, $day, $yr, $hrs, $min);
	
	$dtstring =~ s/ AM//gi;
	if (($mon, $day, $yr, $hrs, $min) = $dtstring =~ /(.+?)\/(.+?)\/(.+?) (.+?):(.+?)/) {
		$yr = "20$yr" if length $yr == 2;
		$mon = "0$mon" if length $mon == 1;
		$day = "0$day" if length $day == 1;
		$dtHash{"date"} = "$yr-$mon-$day";

		$hrs = "0$hrs" if length $hrs == 1;
		$min = "0$min" if length $min == 1;
		$dtHash{"time"} = "$hrs:$min";

	} elsif (($mon, $day, $yr) = $dtstring =~ /(.+?)\/(.+?)\/(.+)/) {
		$yr = "20$yr" if length $yr == 2;
		$mon = "0$mon" if length $mon == 1;
		$day = "0$day" if length $day == 1;
		$dtHash{"date"} = "$yr-$mon-$day";
	}
	
	return \%dtHash;
}



