#! /usr/bin/env perl

use strict;
use warnings;
use DateTime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use Statistics::LineFit;

use Data::Dumper;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] DBDIR\n";
$usage   .=   "Creates/updates the result table in DBDIR, using database files found in DBDIR.\n";
$usage   .=   "\n";


my $dbdir;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } else {
		$dbdir = $arg;
	}
}

die "FATAL: $progname requires a valid data directory.\n$usage\n" unless defined $dbdir and -d "$dbdir";

my $status = 0;

my $LOD   = 2855;



# key columns for WaTCH data
my %table2key = ("abatch"        => "assay_batch_id", 
								 "assay"         => "assay_id", 
								 "cbatch"        => "concentration_batch_id", 
								 "concentration" => "concentration_id", 
								 "ebatch"        => "extraction_batch_id", 
								 "extraction"    => "extraction_id", 
								 "sample"		   	 => "sample_id", 
								 "result"        => "assay_id");

# Existing WaTCH data keyed by WaTCH table
# Sub-hashes for each table are keyed by uid
my %table2watch = ("abatch"        => {}, 
									 "assay"         => {}, 
									 "cbatch"        => {}, 
									 "concentration" => {}, 
									 "ebatch"        => {}, 
									 "extraction"    => {},
									 "sample"		     => {},
									 "result"        => {});



# Read data tables into hash
foreach my $table (keys %table2key) {
	#print "DEBUG::: $progname Reading $table.\n";
	my $keyname = $table2key{"$table"};
	my @colnames = ();
	my $keycol  = -1;
	my $linenum = 0;
	if (-f"$dbdir/watchdb.${table}.txt") {
		open (my $dbFH, "<", "$dbdir/watchdb.${table}.txt") or die "Unable to open $dbdir/watchdb.${table}.txt for reading: $!\n";
		while (my $line = <$dbFH>) {
			chomp $line;
			next if $line =~ /^\s*$/;
			my @cols = split "\t", "$line", -1;
			if ($linenum == 0) {
				# First line of the file contains the column names
				# Store them and determine which column matches the key
				foreach (my $i=0; $i<scalar(@cols); $i++) {
					push @colnames, "$cols[$i]";
					$keycol = $i if "$keyname" eq "$cols[$i]";
				}
			} else {
				# Extract the data for this row, keyed by the ID
				my $thisId = "$cols[$keycol]";
				$table2watch{"$table"}->{"$thisId"} = {};
				for (my $i=0; $i<scalar(@cols); $i++) {
					if (!defined $colnames[$i]) {
						print "DEBUG:: $table col name ($i) does not exist!\n";
						print "$line\n";
						print Dumper(@colnames);
						die;
					}
					$table2watch{"$table"}->{"$thisId"}->{"$colnames[$i]"} = "$cols[$i]";
				}
			}
			$linenum++;
		}
		close $dbFH;
	} else {
		print "WARN : File $dbdir/watchdb.${table}.txt does not exist or can not be read.\n";
	}
}


# WaTCH resource tables read into their own hash
my %resources = ();

# Read resource tables into hash
my $resource_wkbk = ReadData("resources/watchdb.all_tables.xlsx", dtfmt => "mm/dd/yy");
#print Dumper($resource_wkbk);
#die;

foreach my $sheet_name (keys %{$resource_wkbk->[0]->{"sheet"}}) {
	#print "$sheet_name\n";
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

	# The keys for each resource sheet are in column A (or, array 1 of 'cell'). 
	# Loop over these to add entries to the resources hash.
	for (my $i=2; $i < scalar(@{$sheetRef->{"cell"}->[1]}); $i++) {
		#print "$sheetRef->{cell}->[1]->[$i]\n";
		next unless defined $sheetRef->{"cell"}->[1]->[$i] and "$sheetRef->{cell}->[1]->[$i]" ne "";
		my $key = "$sheetRef->{cell}->[1]->[$i]";
		$resources{"$sheet_name"}->{"$key"} = {};
	}
	
	# To get the actual values, loop over the 'cell' array.
	# Get the key from array 1.
	# Loop over the field_names array to get the values.
	for (my $i=1; $i < scalar(@{$sheetRef->{"cell"}}); $i++) {
		next unless defined $sheetRef->{"cell"}->[$i];
		my $colRef = $sheetRef->{"cell"}->[$i];
		my $field = "$field_names[$i-1]";
		my $resourceRef = $resources{"$sheet_name"};
		# Loop over the values in this column. Add each to the resources hash.
		for (my $j=2; $j < scalar(@{$colRef}); $j++) {
			my $key = "$sheetRef->{cell}->[1]->[$j]";
			my $value = "";
			$value = "$colRef->[$j]" if defined $colRef->[$j];
			$resources{"$sheet_name"}->{"$key"}->{"$field"} = "$value";
		}
	}

}

#print Dumper(\%resources);
#die;




# Every assay has one and only one corresponding result
# Samples may have more than one assay (and result)
#
# Required input data:
# assay.assay_target (eg SARS-CoV-2)
# assay.assay_target_category (eg Pathogen)
# assay.assay_target_genetic_locus (eg N1)
# assay.assay_target_copies_per_ul_reaction
# assay.assay_input_ul
# abatch.assay_reaction_ul
# cbatch.concentration_input_ml
# cbatch.concentration_output_ml
# ebatch.extraction_input_ul
# ebatch.extraction_output_ul
# sample.sample_qc
# sample.location_id
# sample.sample_flow
# location.location_population_served

# LOOKUPS
# assay.extraction_id -> extraction.extraction_id
# extraction.concentration_id -> concentration.concentration_id
# concentration.sample_id -> sample.sample_id
# sample.location_id -> location.location_id

# R dashboard will import 5 tables:
# result, location, county, wwtp, and sample



#
# Make sure all required keys have some value for each asset
#
=cut
foreach my $asset_id (keys %asset2data) {
	next unless scalar keys %{$asset2data{"$asset_id"}} > 0;	# The delete function leaves the key so ignore empty assets.
	my $asset = $asset2data{"$asset_id"};

	foreach my $key (keys %required_sample_keys) {
		$asset->{"$key"} = "NA" unless defined $asset->{"$key"};

		# Re-format all dates
		if ($asset->{"$key"} =~ /(\d{1,2})\/(\d{1,2})\/(\d{2,4}).*/) {
			$asset->{"$key"} = "$1/$2/$3";
		} elsif ($asset->{"$key"} =~ /(\d{1,2})\-(\d{1,2})\-(\d{2,4}).*/) {
			$asset->{"$key"} = "$1/$2/$3";
		}
	}
	
	foreach my $key (keys %required_sampleQC_keys) {
		$asset->{"$key"} = "NA" unless defined $asset->{"$key"};
	}
}
=cut

#print Dumper(\%table2watch);
#die;


# Calculate the results (copies per L of wastewater and normalized derivatives) for each assay
foreach my $assay_id (keys %{$table2watch{"assay"}}) {
	
	# If this assay id is already present in the result table, and has been validated, skip it.
	next if defined $table2watch{"result"}->{"$assay_id"} and $table2watch{"result"}->{"$assay_id"}->{"target_result_validated"} == 1;
	
	# Generate various data hashes related to this assay
	my %assay = %{$table2watch{"assay"}->{"$assay_id"}};
	my $exid = $assay{"extraction_id"};
	my $abid = $assay{"assay_batch_id"};

	my %abatch = %{$table2watch{"abatch"}->{"$abid"}};

	my %extraction = %{$table2watch{"extraction"}->{"$exid"}};
	my $cnid = $extraction{"concentration_id"};
	my $ebid = $extraction{"extraction_batch_id"};

	my %ebatch = %{$table2watch{"ebatch"}->{"$ebid"}};

	my %concentration = %{$table2watch{"concentration"}->{"$cnid"}};
	my $smid = $concentration{"sample_id"};
	my $cbid = $concentration{"concentration_batch_id"};
	
	if (!defined $table2watch{"sample"}->{"$smid"}) {
		print "##########################################################################\n";
		print "#\n";
		print "WARN: $smid was not found in the SAMPLE table! This assay will be skipped.\n";
		print "      assay id        : $assay_id\n";
		print "      extraction id   : $exid\n";
		print "      concentration id: $cnid\n";
		print "#\n";
		print "##########################################################################\n";
		next;
	}
	
	my %sample = %{$table2watch{"sample"}->{"$smid"}};
	my $lcid = $sample{"location_id"};
	
	if (!defined $resources{"location"}->{"$lcid"}) {
		print "##########################################################################\n";
		print "#\n";
		print "WARN: $lcid was not found in the LOCATION resource table! This assay will be skipped.\n";
		print "      assay id        : $assay_id\n";
		print "      extraction id   : $exid\n";
		print "      concentration id: $cnid\n";
		print "      sample id       : $smid\n";
		print "#\n";
		print "##########################################################################\n";
		next;
	}

	my %cbatch = %{$table2watch{"cbatch"}->{"$cbid"}};
	
	# Initialize the result entry, keyed by assay id
	$table2watch{"result"}->{"$assay_id"} = {"assay_id"                  => "$assay_id", 
																					 "sample_id"                 => $smid,
																					 "target"                    => $assay{"assay_target"},
																					 "target_category"           => $assay{"assay_target_category"},
																					 "target_genetic_locus"      => $assay{"assay_target_genetic_locus"}, 
																					 "lab_id"   								 => "ZooWVU",
																					 "location_id"               => "$lcid",
																					 "collection_start_datetime" => $table2watch{"sample"}->{"$smid"}->{"sample_collection_start_datetime"}, 
																					 "collection_end_datetime"   => $table2watch{"sample"}->{"$smid"}->{"sample_collection_end_datetime"}, 
																					 "sample_flow" => $table2watch{"sample"}->{"$smid"}->{"sample_flow"}, 
																					 "assay_target_copies_per_ul_reaction" => "NA",
																					 "assay_reaction_ul"									 => "NA",
																					 "extraction_output_ul"								 => "NA",
																					 "concentration_output_ml"						 => "NA",
																					 "assay_input_ul"											 => "NA",
																					 "extraction_input_ul"								 => "NA",
																					 "concentration_input_ml"							 => "NA",
																					 "target_copies_per_l"								 => "NA",
																					 "target_copies_per_lday"							 => "NA",
																					 "target_copies_flownorm"							 => "NA",
																					 "target_copies_flownorm_per_person"	 => "NA",
																					 "target_copies_per_ldayperson"				 => "NA",
																					 "target_result_validated"					   => 0
																					};
	
	## NEED to check all of these values exist. Otherwise set status and skip the calculation.
	if (defined $assay{"assay_target_copies_per_ul_reaction"} and isEmpty($assay{"assay_target_copies_per_ul_reaction"}) == 0) {
		my $copies_per_lww = 1000 * 1000 * 
			($assay{"assay_target_copies_per_ul_reaction"} * $abatch{"assay_reaction_ul"} * $ebatch{"extraction_output_ul"} * $cbatch{"concentration_output_ml"}) / 
			($assay{"assay_input_ul"} * $ebatch{"extraction_input_ul"} * $cbatch{"concentration_input_ml"});

#		$copies_per_lww = 1 + int rand($LOD) if $copies_per_lww < 1;
		
		$table2watch{"result"}->{"$assay_id"}->{"assay_target_copies_per_ul_reaction"} = $assay{"assay_target_copies_per_ul_reaction"};
		$table2watch{"result"}->{"$assay_id"}->{"target_copies_per_l"} = $copies_per_lww;

		$resources{"location"}->{"$lcid"}->{"location_population_served"} =~ s/,//i;
		$table2watch{"sample"}->{"$smid"}->{"sample_flow"} =~ s/,//i;
		
		$table2watch{"result"}->{"$assay_id"}->{"target_copies_per_lday"} = 
			$copies_per_lww * 24 * (1/$resources{"location"}->{"$lcid"}->{"location_collection_window_hrs"});
	
		if (defined $table2watch{"sample"}->{"$smid"}->{"sample_flow"} and $table2watch{"sample"}->{"$smid"}->{"sample_flow"} ne "" and $table2watch{"sample"}->{"$smid"}->{"sample_flow"} ne "NA") {
			# copies/liter x 3.78541 liters/gallon x 10^6 gallon/million_gallons x million_gallons/day x 1 day/24 hrs x collection_time_hrs

			$table2watch{"result"}->{"$assay_id"}->{"target_copies_flownorm"} = 
				$copies_per_lww * 3.78541 * 1000000 * 
				$table2watch{"sample"}->{"$smid"}->{"sample_flow"} * 
				($resources{"location"}->{"$lcid"}->{"location_collection_window_hrs"} / 24);

			if (defined $resources{"location"}->{"$lcid"}->{"location_population_served"} and $resources{"location"}->{"$lcid"}->{"location_population_served"} ne "" and $resources{"location"}->{"$lcid"}->{"location_population_served"} ne "NA") {
				$table2watch{"result"}->{"$assay_id"}->{"target_copies_flownorm_per_person"} = 
					$table2watch{"result"}->{"$assay_id"}->{"target_copies_flownorm"} / $resources{"location"}->{"$lcid"}->{"location_population_served"};
			}
		}

		# copies/liter x 24 hrs/day x 1/collection_time_hrs x 1/population_served
		if (defined $resources{"location"}->{"$lcid"}->{"location_population_served"} and $resources{"location"}->{"$lcid"}->{"location_population_served"} ne "" and $resources{"location"}->{"$lcid"}->{"location_population_served"} ne "NA") {
			$table2watch{"result"}->{"$assay_id"}->{"target_copies_per_ldayperson"} = 
					$copies_per_lww * 24 * 
					(1/$resources{"location"}->{"$lcid"}->{"location_collection_window_hrs"}) * 
					(1/$resources{"location"}->{"$lcid"}->{"location_population_served"});
		}
		## NEED some better way to validate this result. Compare to control?
		$table2watch{"result"}->{"$assay_id"}->{"target_result_validated"}   = 1;
	}
}

#my $cn_per_l = 1000 * 
#(CN/RXN * $calc_vals{"Extraction Output Volume (mL)"} * $calc_vals{"Concentration Output Volume (mL)"}) / 
#($calc_vals{"Assay Volume (mL)"} * $calc_vals{"Extraction Input Volume (mL)"} * $calc_vals{"Concentration Input Volume (mL)"});


# Add data from Marshall University to result file
=cut
my @keysWmu = ();
$count      = 0;

if (-f "$WATCHFILE_MU") {
	open (my $watchMuInFH, "<", "$WATCHFILE_MU") or die "Unable to open WATCHFILE_MU for reading: $!\n";
	while (my $line = <$watchMuInFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @values = split "\t", "$line", -1;

		if ($count == 0) {
			# populate the keys array with values on the first line
			for (my $j=0; $j < scalar(@values); $j++) {
				$keysWmu[$j] = convertKey($values[$j]);
			}
		} else {
			# first field is the asset ID
			my $asset_id = "MU-$values[0]";
			# populate the hash entry for this asset with the remaining values
			$asset2data{$asset_id} = {} unless defined $asset2data{$asset_id};
			for (my $j=0; $j < scalar(@values); $j++) {
				$values[$j] = "NA" if isEmpty("$values[$j]") == 1;
				$asset2data{"$asset_id"}->{"$keysWmu[$j]"} = "$values[$j]" if defined $required_sample_keys{"$keysWmu[$j]"} or defined $required_sampleQC_keys{"$keysW[$j]"};
			}
			$asset2data{$asset_id}->{"responsible_lab"} = "Marshall University Infectious Disease Surveillance Lab";
			$asset2data{$asset_id}->{"responsible_lab_abbrev"} = "MUIDSL";
		}
		$count++;
	}
	close $watchMuInFH;
}
=cut



#
# Write the result file
#

my @result_colnames = ("assay_id",
											 "sample_id",
											 "collection_start_datetime",
											 "collection_end_datetime",
											 "sample_flow",
											 "sample_qc", 
											 "location_id",
											 "target",
											 "target_category",
											 "target_genetic_locus", 
											 "lab_id",
											 "target_copies_per_l",
											 "target_copies_per_lday",
											 "target_copies_flownorm",
											 "target_copies_flownorm_per_person",
											 "target_copies_per_ldayperson",
											 "target_result_validated");
											 
open (my $rFH, ">", "$dbdir/watchdb.result.txt") or die "Unable to open $dbdir/watchdb.result.txt for writing: $!";
print $rFH join("\t", @result_colnames) . "\n";
foreach my $assay_id (keys %{$table2watch{"result"}}) {
	print $rFH "$assay_id";
	for (my $i=1; $i < scalar(@result_colnames); $i++) {
		my $colname = $result_colnames[$i];
		my $colval  = "NA";
		if (defined $table2watch{"result"}->{"$assay_id"}->{"$colname"} and isEmpty($table2watch{"result"}->{"$assay_id"}->{"$colname"}) == 0) {
			$colval = trim($table2watch{"result"}->{"$assay_id"}->{"$colname"});
		}
		print $rFH "\t$colval";
	}
	print $rFH "\n";
}
close $rFH;


exit $status;


=cut
sub extract {
	my $target = shift;
	my $scope  = shift;
	my $endpts = shift;
	my $asset_id = shift;
	my $asset2sum_ref = shift;
	my $val_by_position_ref = shift;
	
	my @vals = ();
	for (my $j=$endpts->[1]; $j>=$endpts->[0]; $j--) {
		unless ("$val_by_position_ref->[$j]->{$target}" eq "NA") {
			push @vals, $val_by_position_ref->[$j]->{"$target"};
		}
	}

	my $pos = $endpts->[0]-1;
	while (scalar(@vals) <= ($endpts->[1] - $endpts->[0]) and $pos >= 0) {
		unless ("$val_by_position_ref->[$pos]->{$target}" eq "NA") {
			push @vals, $val_by_position_ref->[$pos]->{"$target"};
		}
		$pos--;
	}
	
	if (scalar @vals > 1) {
		my ($m, $ci, $slope) = (0, 0, 0);
		print "$asset_id\t$target\t" if $debug == 1;
		if (scalar(@vals) > 0) {
			$m = mean(\@vals);
			$ci = ci(\@vals, 0.98);
			$slope = regress(\@vals);
		}
		print "\n" if $debug == 1;
		$asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.mean"}  = $m;
		$asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.ci"}    = $ci;
		#$asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.slope"} = $slope;
	} else {
		print "$asset_id: $target: $scope: $endpts->[0], $endpts->[1]. Less than 2 values for df calculation.\n" unless "$target" =~ /load/;
		#print Dumper($val_by_position_ref);
		#die;
	}
	
#	if ($debug == 1 and $asset2data{"$asset_id"}->{"Location"} eq "StarCityWWTP-01") {
#		print $asset2data{"$asset_id"}->{"Sample Composite End"} . "\t$scope.$target\t";
#		print "mean: " . $asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.mean"} . "\tCI: " . $asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.ci"};
#		print "\t" . scalar(@vals) . " vals: " . join(",", @vals) . "\n";
#	}
}
=cut

sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


sub isEmpty {
	my $val = shift;
	
	my $is_empty = 0;
	$is_empty = 1 if "$val" eq "NaN" or "$val" eq "-" or "$val" eq "none" or "$val" eq "" or "$val" eq "TBD" or "$val" eq "NA";
	
	return $is_empty;
}

=cut
sub convertKey {
	my $keyIn = shift;
	my %converter = ("N1" => "Assay Target 1 Result (CN/L)",
									 "N2" => "Assay Target 2 Result (CN/L)",
									 "sample composite start" => "Sample Composite Start",
									 "sample composite end" => "Sample Composite End",
									 "sample received date" => "Sample Received Date",
									 "Sample flow (MGD)" => "Sample Flow (MGD)",
									 "location" => "Location");
	
	if (defined $converter{"$keyIn"}) {
		return $converter{"$keyIn"};
	} else {
		return "$keyIn";
	}
	
}

sub mean {
	my $vals = shift;
	my $n = scalar(@$vals);
	my $sum = 0;
	foreach my $val (@$vals) {
		$sum += $val;
	}
	my $m = $sum/$n;
	return $m;
}

sub sd {
	my $vals = shift;

	my $n = scalar(@$vals);
	my $m = mean($vals);
	print "m=$m\t" if $debug == 1;
	my $sum = 0;
	
	foreach my $val (@$vals) {
		$sum += ($val - $m) * ($val - $m);
	}
	my $sdev = sqrt($sum/$n);
	return $sdev;
}

sub ci {
	my $vals = shift;
	my $level = shift;
	
	my $n = scalar(@$vals);
	my $df = $n-1;
	my $alpha = 0.5 * (1 - $level);	# 90% = 0.05, 95% = 0.025, 98% = 0.002
	my $t = $TDIST_2TAIL{$df}->{$alpha};
	if (!defined $t) {
		print "df: $df\t, alpha: $alpha\n";
		die;
	}
	my $sd = sd($vals);
	my $q = $sd/sqrt($n);
	my $conf = $t * $q;
	
	if ($debug == 1) {
		print "n=$n\t";
		print "a=$alpha\t";
		print "s=$sd\t";
		print "q=$q\t";
		print "t=$t\t";
		print "conf=$conf\t";
	}
		
	return $conf;
}

sub regress {
	my $vals = shift;
	
	my $lineFit = Statistics::LineFit->new();
	my @x = ();
	my @y = ();
	for (my $i=0; $i<scalar(@$vals); $i++) {
		my $val = $vals->[$i];
		push @x, $i;
		push @y, $val;
	}
	@y = reverse(@y);
	
	$lineFit->setData(\@x, \@y);
  my ($intercept, $slope) = $lineFit->coefficients();

	return $slope;
}
=cut



