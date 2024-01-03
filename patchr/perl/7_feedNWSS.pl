#! /usr/bin/env perl

use strict;
use warnings;

use Date::Format;
use DateTime qw( );
#use DateTime::Format::Strptime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [DBDIR]\n";
$usage   .=  "Build feed table for NWSS using WaTCHdb tables in DBDIR (data/latest/ by default) ";
$usage   .=  "and writes to STDOUT.\n";
$usage   .=   "\n";


my $status = 0;
my $dbdir  = "data/latest";


while (@ARGV) {
  my $arg = shift;
  if ("$arg" eq "-h") {
		die "$usage";
  } else {
		$dbdir = "$arg";
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



my $feed         = "NWSS";
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


# Get the controlled vocabulary map for this feed.
# This is located in the resources/ folder.
my $linenum = 0;
open (my $vFH, "<", "resources/${feed}_cvm.txt") or die "FATAL: Unable to open resources/${feed}_cvm.txt for reading: $!.";
while (my $line = <$vFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	#my @cols = split /\t/, "$line", -1;
	my ($nwss_field, $nwss_value, $watch_table, $watch_field, $watch_value, $reference) = split /\t/, "$line", -1;
	if ($linenum == 0) {
		$linenum++;
	} else {
		$feed_cvmap{"$nwss_field"} = {} unless defined $feed_cvmap{"$nwss_field"};
		$feed_cvmap{"$nwss_field"}->{"$watch_field"} = {} unless defined $feed_cvmap{"$nwss_field"}->{"$watch_field"};
		$feed_cvmap{"$nwss_field"}->{"$watch_field"}->{"$watch_value"} = {"value" => "$nwss_value", "reference" => "$reference"};
	}
}
close $vFH;


my %tables = (
	"abatch"				=> {},
	"assay" 				=> {},
	"cbatch" 				=> {},
	"concentration" => {},
	"ebatch" 				=> {},
	"extraction" 		=> {},
	"result" 				=> {},
	"sample"				=> {}
);

my %resources = (
	"location" => {},
	"wwtp"		 => {},
	"county"	 => {},
	"lab"			 => {}
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
# Read data from WaTCH-NWSS mapping table. Keyed by NWSS column header.
#
my %map = ();

my $map_wkbk = ReadData("resources/${feed}_map.xlsx", dtfmt => "mm/dd/yy");
my @map_rows = Spreadsheet::Read::rows($map_wkbk->[1]);

my %map_headers = ();
for (my $i=0; $i < scalar(@{$map_rows[0]}); $i++) {
	$map_headers{$i} = $map_rows[0][$i];
}

for (my $i=1; $i < scalar(@map_rows); $i++) {
	my $nwss_field = $map_rows[$i][1];
	$map{$nwss_field} = {};
	for (my $j=0; $j < scalar(@{$map_rows[$i]}); $j++) {
		$map{$nwss_field}->{"$map_headers{$j}"} = $map_rows[$i][$j];
	}
}
#print Dumper(\%map);
#die;

#print Dumper(\@feed_headers);
#die;

# To assemble the output file, we loop over entries in the result table.
# For each result, we loop over @feed_headers in order, assemble the line, and print it.

#open(my $oFH, ">", "$OUTFILE_BK") or die "Unable to open $OUTFILE_BK for writing: $!";
#print $oFH "\"" . join("\",\"", @feed_headers) . "\"\n";
print "\"" . join("\",\"", @feed_headers) . "\"\n";

RSLTLOOP:
foreach my $assay_id (keys %{$tables{"result"}}) {
	my $resultRef = $tables{"result"}->{"$assay_id"};
	
	# A few short-circuits:
	#
	# Skip results with NTC contamination.
	next if lc("$resultRef->{target_result_validated}") eq "ntc above threshold";
	#
	# Skip if pcr_target is not in NWSS controlled vocabulary.
	next unless defined $feed_cvmap{"pcr_target"}->{"assay_target"}->{"$resultRef->{target}"};
	#
	# Skip if pcr_gene_target is not in NWSS controlled vocabulary.
	next unless defined $feed_cvmap{"pcr_target"}->{"assay_target"}->{"$resultRef->{target_genetic_locus}"};
	#	
	# Standardize results with no flow.
	$resultRef->{"sample_flow"} = "" if "$resultRef->{sample_flow}" eq "NA" or "$resultRef->{sample_flow}" eq "0";

	my @result_vals = ();
	my $val_pass    = 0;
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
			$val = lookupLinked("$assay_id", "$table", "$mapEntry{linked_field}");
		} elsif ("$mapEntry{WaTCH_field}" eq "[calculated]") {
			# Call the calcNWSS sub to figure this out.
			my $lfield = "";
			$lfield = "$mapEntry{linked_field}" if defined $mapEntry{"linked_field"};
			$val = calcNWSS("$assay_id", "$header", "$lfield");
		} elsif ("$mapEntry{WaTCH_table}" eq "result") {
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
		$val_pass = 1 unless "$val" eq "";
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
	my $aid    = shift;
	my $tname  = shift;
	my $tfield = shift;
	
	my $table_cat = "table";
	my $linkedId;

	if ("$tname" eq "abatch") {
		$linkedId = $tables{"assay"}->{"$aid"}->{"assay_batch_id"};
	} elsif ("$tname" eq "cbatch") {
		my $eid = $tables{"assay"}->{"$aid"}->{"extraction_id"};
		my $cid = $tables{"extraction"}->{"$eid"}->{"concentration_id"};
		$linkedId = $tables{"concentration"}->{"$cid"}->{"concentration_batch_id"};
	} elsif ("$tname" eq "ebatch") {
		my $eid = $tables{"assay"}->{"$aid"}->{"extraction_id"};
		$linkedId = $tables{"extraction"}->{"$eid"}->{"extraction_batch_id"};
	} elsif ("$tname" eq "location") {
		$table_cat = "resource";
		$linkedId = $tables{"result"}->{"$aid"}->{"location_id"};
	} elsif ("$tname" eq "wwtp") {
		$table_cat = "resource";
		my $locid = $tables{"result"}->{"$aid"}->{"location_id"};
		$linkedId = $resources{"location"}->{"$locid"}->{"location_primary_wwtp_id"};
	} elsif ("$tname" eq "county") {
		$table_cat = "resource";
		my $locid = $tables{"result"}->{"$aid"}->{"location_id"};
		$linkedId = $resources{"location"}->{"$locid"}->{"location_counties_served"};
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


sub calcNWSS {
	my $aid          = shift;
	my $nwss_field   = shift;
	my $linked_field = shift;
	
	my $val = "";
	
	if ("$nwss_field" eq "institution_type") {
		# location.location_category but 'not institution specific' if wwtp
		my $locid = $tables{"result"}->{"$aid"}->{"location_id"};
		$val = $resources{"location"}->{"$locid"}->{"location_category"};
		$val = "not institution specific" if lc($val) eq "wwtp";
	} elsif ("$nwss_field" eq "major_lab_method") {
		# Test abatch.assay_quantification_method
		my $abid = $tables{"assay"}->{"$aid"}->{"assay_batch_id"};
		if (lc "$tables{abatch}->{$abid}->{assay_quantification_method}" eq "ddpcr") {
			$val = 1;
		} else {
			$val = 2;
		}
	} elsif ("$nwss_field" eq "major_lab_method_desc") {
		# Test abatch.assay_quantification_method
		my $abid = $tables{"assay"}->{"$aid"}->{"assay_batch_id"};
		if (lc "$tables{abatch}->{$abid}->{assay_quantification_method}" eq "ddpcr") {
			$val = "";
		} else {
			$val = "This sample was quantified by real-time qPCR";
		}
	} elsif ("$nwss_field" eq "ntc_amplify") {
		# if result.target_result_validated eq 'NTC above threshold' then yes; otherwise no
		if (lc("$tables{result}->{$aid}->{target_result_validated}") eq "ntc above threshold") {
			$val = "yes";
		} else {
			$val = "no";
		}
	} elsif ("$nwss_field" =~ /^pcr.+target(.*)$/) {
		my $watch_val = $tables{"assay"}->{"$aid"}->{"$linked_field"};
		if ("$1" eq "_ref") {
			$val = "$feed_cvmap{$nwss_field}->{$linked_field}->{$watch_val}->{reference}" if defined $feed_cvmap{"$nwss_field"}->{"$linked_field"}->{"$watch_val"}->{"reference"};
		} else {
			$val = "$feed_cvmap{$nwss_field}->{$linked_field}->{$watch_val}->{value}" if defined $feed_cvmap{"$nwss_field"}->{"$linked_field"}->{"$watch_val"}->{"value"};
		}
	} elsif ("$nwss_field" eq "sample_location") {
		my $locid = $tables{"result"}->{"$aid"}->{"location_id"};
		$val = $resources{"location"}->{"$locid"}->{"location_category"};
		$val = "upstream" unless lc($val) eq "wwtp";
	} elsif ("$nwss_field" eq "sample_location_specify") {
		my $locid = $tables{"result"}->{"$aid"}->{"location_id"};
		my $cat  = $resources{"location"}->{"$locid"}->{"location_category"};
		my $name = $resources{"location"}->{"$locid"}->{"location_common_name"};
		$val = "$name" if lc($cat) ne "wwtp";
	} elsif ("$nwss_field" eq "sample_type") {
		my $locid = $tables{"result"}->{"$aid"}->{"location_id"};
		my $hrs = $resources{"location"}->{"$locid"}->{"location_collection_window_hrs"};
		my $bas = $resources{"location"}->{"$locid"}->{"location_collection_basis"};
		$val = "$hrs" . "-hr " . "$bas" . "-weighted composite";
	} elsif ("$nwss_field" eq "extraction_method") {
		my $eid  = $tables{"assay"}->{"$aid"}->{"extraction_id"};
		my $ebid = $tables{"extraction"}->{"$eid"}->{"extraction_batch_id"};
		my $wval = $tables{"ebatch"}->{"$ebid"}->{"$linked_field"};
		$val = "$feed_cvmap{$nwss_field}->{$linked_field}->{$wval}->{value}" if defined $feed_cvmap{"$nwss_field"}->{"$linked_field"}->{"$wval"}->{"value"};
	} elsif ("$nwss_field" eq "concentration_method") {
		my $eid  = $tables{"assay"}->{"$aid"}->{"extraction_id"};
		my $cid  = $tables{"extraction"}->{"$eid"}->{"concentration_id"};
		my $cbid = $tables{"concentration"}->{"$cid"}->{"concentration_batch_id"};
		my $wval = $tables{"cbatch"}->{"$cbid"}->{"$linked_field"};
		$val = "$feed_cvmap{$nwss_field}->{$linked_field}->{$wval}->{value}" if defined $feed_cvmap{"$nwss_field"}->{"$linked_field"}->{"$wval"}->{"value"};
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



