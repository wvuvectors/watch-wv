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
$usage   .=   "Generates a result table using database files found in DBDIR.\n";
$usage   .=   "\n";


my $dbdir;

while (@ARGV) {
  my $arg = shift;
  if ("$arg" eq "-h") {
		die $usage;
  } else {
		$dbdir = "$arg";
	}
}

die "FATAL: $progname requires a valid data directory.\n$usage\n" unless defined $dbdir and -d "$dbdir";

my $status   = 0;
my $errlimit = 20;

my $LOD           = 2855;	# copies of target per L
my $NTC_THRESHOLD = 10;		# copies of target per PCR reaction



# key columns for WaTCH data
my %table2key = ("abatch"        => "assay_batch_id", 
								 "assay"         => "assay_id", 
								 "cbatch"        => "concentration_batch_id", 
								 "concentration" => "concentration_id", 
								 "ebatch"        => "extraction_batch_id", 
								 "extraction"    => "extraction_id", 
								 "result"        => "assay_id", 
								 "sample"		   	 => "sample_id");

# Existing WaTCH data keyed by WaTCH table
# Sub-hashes for each table are keyed by uid
my %table2watch = ("abatch"        => {}, 
									 "assay"         => {}, 
									 "cbatch"        => {}, 
									 "concentration" => {}, 
									 "ebatch"        => {}, 
									 "extraction"    => {},
									 "result"        => {},
									 "sample"		     => {});


# Assay controls are stored in their own hash.
# Each control is an array of copies per ul reaction, keyed by AB id, dye, and control type.
my %controls = ();

# Read data tables into hash
foreach my $table (keys %table2key) {
	#print "DEBUG::: $progname Reading $table.\n";
	my $keyname = $table2key{"$table"};
	my @colnames = ();
	my $keycol  = -1;
	my $linenum = 0;
	if (-f "$dbdir/watchdb.${table}.txt") {
		open (my $dbFH, "<", "$dbdir/watchdb.${table}.txt") or die "Unable to open $dbdir/watchdb.${table}.txt for reading: $!\n";
		while (my $line = <$dbFH>) {
			chomp $line;
			next if "$line" =~ /^\s*$/;
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
						exit -1;
					}
					$table2watch{"$table"}->{"$thisId"}->{"$colnames[$i]"} = "$cols[$i]";
				}
			}
			$linenum++;
		}
		close $dbFH;
	} else {
		if ("$table" eq "result") {
			print "INFO : There is no previous result table in $dbdir.\n";
			print "INFO : That's not fatal, because I generate a complete result file anyway.\n";
			print "INFO : However, it means I may be changing any previous validated results.\n";
			print "INFO : Just so you know.\n";
		} else {
			print "WARN : File $dbdir/watchdb.${table}.txt does not exist or can not be read.\n";
		}
	}
}


#print Dumper(\%table2watch);
#die;


# WaTCH resource tables read into their own hash
my %resources = ();

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
			next unless defined $sheetRef->{"cell"}->[1]->[$j] and "$sheetRef->{cell}->[1]->[$j]" ne "";
			my $key = "$sheetRef->{cell}->[1]->[$j]";
			my $value = "";
			$value = "$colRef->[$j]" if defined $colRef->[$j];
			$resources{"$sheet_name"}->{"$key"}->{"$field"} = "$value";
		}
	}
}

#print Dumper(\%resources);
#die;




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
foreach my $assay_id (keys %{$table2watch{"assay"}}) {
	my $ahashRef = $table2watch{"assay"}->{"$assay_id"};

	# By default, this assay will be processed into a result.
	$ahashRef->{"reject"} = "no";

	# KLUDGE! Skip if this assay is a pool.
	$ahashRef->{"reject"} = "pool" if "$assay_id" =~ /^.+\-P.+$/;
	
	# If this assay is a control, add it to the controls hash and skip.
	if ($ahashRef->{"assay_class"} eq "control") {
		$ahashRef->{"reject"} = "control";

		my $abid  = $ahashRef->{"assay_batch_id"};
		my $atype = $ahashRef->{"assay_type"};
		my $dye   = $ahashRef->{"assay_target_fluorophore"};
		my $cpr   = $ahashRef->{"assay_target_copies_per_ul_reaction"};
		unless ($cpr =~ /^(\d|\.)+$/gi) {
			print Dumper($ahashRef);
			die;
		}
		
		$controls{"$abid"} = {} unless defined $controls{"$abid"};
		$controls{"$abid"}->{"$dye"} = {} unless defined $controls{"$abid"}->{"$dye"};
		$controls{"$abid"}->{"$dye"}->{"$atype"} = [] unless defined $controls{"$abid"}->{"$dye"}->{"$atype"};
		push @{$controls{"$abid"}->{"$dye"}->{"$atype"}}, $cpr;
		
	}
	
	# Assay id is already present in the result table and has been validated.
	$ahashRef->{"reject"} = "pre-existing" if defined $table2watch{"result"}->{"$assay_id"} and lc "$table2watch{result}->{$assay_id}->{target_result_validated}" eq "yes";
	
	next unless "$ahashRef->{reject}" eq "no";

	# Sample for this assay has been rejected or failed initial QC.
	$ahashRef->{"reject"} = "sample rejected" if "$table2watch{sample}->{$ahashRef->{sample_id}}->{sample_status}" =~ /Rejected/i;
	$ahashRef->{"reject"} = "sample rejected" if "$table2watch{sample}->{$ahashRef->{sample_id}}->{sample_qc}" =~ /Fail/i;
	
	next unless $ahashRef->{"reject"} eq "no";
	
	$ahashRef->{"reject"} = "abatch id problem" if 
			!defined $ahashRef->{"assay_batch_id"} or 
			"$ahashRef->{assay_batch_id}" eq "" or 
			!defined $table2watch{"abatch"}->{$ahashRef->{"assay_batch_id"}} or 
			"$table2watch{abatch}->{$ahashRef->{assay_batch_id}}" eq "";

	$ahashRef->{"reject"} = "extraction id problem" if 
			!defined $ahashRef->{"extraction_id"} or 
			"$ahashRef->{extraction_id}" eq "" or 
			!defined $table2watch{"extraction"}->{$ahashRef->{"extraction_id"}} or 
			"$table2watch{extraction}->{$ahashRef->{extraction_id}}" eq "";

	next unless $ahashRef->{"reject"} eq "no";
	

	my %ehash = %{$table2watch{"extraction"}->{$ahashRef->{"extraction_id"}}};
	
	$ahashRef->{"reject"} = "ebatch id problem" if 
			!defined $ehash{"extraction_batch_id"} or 
			"$ehash{extraction_batch_id}" eq "" or 
			!defined $table2watch{"ebatch"}->{$ehash{"extraction_batch_id"}} or 
			"$table2watch{ebatch}->{$ehash{extraction_batch_id}}" eq "";

	$ahashRef->{"reject"} = "concentration id problem" if 
			!defined $ehash{"concentration_id"} or 
			"$ehash{concentration_id}" eq "" or 
			!defined $table2watch{"concentration"}->{$ehash{"concentration_id"}} or 
			"$table2watch{concentration}->{$ehash{concentration_id}}" eq "";

	next unless $ahashRef->{"reject"} eq "no";
	

	my %chash = %{$table2watch{"concentration"}->{$ehash{"concentration_id"}}};
	
	$ahashRef->{"reject"} = "cbatch id problem" if 
			!defined $chash{"concentration_batch_id"} or 
			"$chash{concentration_batch_id}" eq "" or 
			!defined $table2watch{"cbatch"}->{$chash{"concentration_batch_id"}} or 
			"$table2watch{cbatch}->{$chash{concentration_batch_id}}" eq "";

	$ahashRef->{"reject"} = "sample id problem" if 
			!defined $chash{"sample_id"} or 
			"$chash{sample_id}" eq "" or 
			!defined $table2watch{"sample"}->{$chash{"sample_id"}} or 
			"$table2watch{sample}->{$chash{sample_id}}" eq "";

	next unless $ahashRef->{"reject"} eq "no";
	

	my %shash = %{$table2watch{"sample"}->{$chash{"sample_id"}}};
	
	$ahashRef->{"reject"} = "location id problem" if 
			!defined $shash{"location_id"} or 
			"$shash{location_id}" eq "" or 
			!defined $resources{"location"}->{$shash{"location_id"}} or 
			"$resources{location}->{$shash{location_id}}" eq "";

}

#print Dumper(\%table2watch);
#print Dumper(\%resources);
#print Dumper(\%controls);
#die;


# Establish a few error logs for this process.
open (my $resFH, ">", "$dbdir/_result.rejected_assays.txt");
print $resFH "assay_id\treason for rejection\n";

open (my $ctlFH, ">", "$dbdir/_result.missing_controls.txt");
print $ctlFH "assay_batch_id\tdye\tassay_date\n";


# Calculate the results (copies per L of wastewater and normalized derivatives) for each assay
foreach my $assay_id (keys %{$table2watch{"assay"}}) {
	my %assay = %{$table2watch{"assay"}->{"$assay_id"}};
	
	if ("$assay{reject}" ne "no") {
		if ("$assay{reject}" =~ "problem") {
			print $resFH "$assay_id\t$assay{reject}\n";
			$status++;
		}
		next;
	}
		
	# Generate various data hashes related to this assay.

	my $exid = $assay{"extraction_id"};
	my %extraction = %{$table2watch{"extraction"}->{"$exid"}};

	my $abid = $assay{"assay_batch_id"};
	my %abatch = %{$table2watch{"abatch"}->{"$abid"}};

	my $ebid = $extraction{"extraction_batch_id"};
	my %ebatch = %{$table2watch{"ebatch"}->{"$ebid"}};

	my $cnid = $extraction{"concentration_id"};
	my %concentration = %{$table2watch{"concentration"}->{"$cnid"}};

	my $cbid = $concentration{"concentration_batch_id"};
	my %cbatch = %{$table2watch{"cbatch"}->{"$cbid"}};
	
	my $smid = $concentration{"sample_id"};
	my %sample = %{$table2watch{"sample"}->{"$smid"}};

	my $lcid = $sample{"location_id"};
	my %location = %{$resources{"location"}->{"$lcid"}};
	
	
	# Initialize the result entry, keyed by assay id
	$table2watch{"result"}->{"$assay_id"} = {"assay_id"                  => "$assay_id", 
																					 "sample_id"                 => $smid,
																					 "target"                    => $assay{"assay_target"},
																					 "target_genetic_locus"      => $assay{"assay_target_genetic_locus"}, 
																					 "lab_id"   								 => "ZooWVU",
																					 "location_id"               => "$lcid",
																					 "event_type"                => $table2watch{"sample"}->{"$smid"}->{"sample_event"}, 
																					 "collection_start_datetime" => $table2watch{"sample"}->{"$smid"}->{"sample_collection_start_datetime"}, 
																					 "collection_end_datetime"   => $table2watch{"sample"}->{"$smid"}->{"sample_collection_end_datetime"}, 
																					 "sample_flow" 							 => $table2watch{"sample"}->{"$smid"}->{"sample_flow"}, 
																					 "sample_qc" 							   => $table2watch{"sample"}->{"$smid"}->{"sample_qc"}, 
																					 "target_copies_per_l"			 => "NA",
																					 "target_copies_per_ld" 		 => "NA",
																					 "target_copies_per_ldcap"   => "NA",
																					 "target_copies_flownorm"		 => "NA",
																					 "target_copies_fn_per_cap"  => "NA",
																					 "target_per_capita_basis"   => 100,
																					 "nc_copies_per_rxn"	       => "NA",
																					 "pc_copies_per_rxn"	       => "NA",
																					 "target_result_validated"	 => ""
																					};
	
	## NEED to check that all of the required values exist; otherwise skip to the next assay.
	next unless doCalcBase(\%assay, \%abatch, \%ebatch, \%cbatch, \%location) == 1;
	
	my $copies_per_lww = calcCopiesPerL(\%assay, \%abatch, \%ebatch, \%cbatch);
	$copies_per_lww = 1 + int rand($LOD) if $copies_per_lww < $LOD;
	$table2watch{"result"}->{"$assay_id"}->{"target_copies_per_l"} = $copies_per_lww;
	
	my $copies_per_ld = $copies_per_lww / ($location{"location_collection_window_hrs"} / 24);
	$table2watch{"result"}->{"$assay_id"}->{"target_copies_per_ld"} = $copies_per_ld;

	if (doCalcPop(\%location) == 1) {
		my $copies_per_ldcap = calcCopiesPop($copies_per_ld, 100, \%location);
		$table2watch{"result"}->{"$assay_id"}->{"target_copies_per_ldcap"} = $copies_per_ldcap;
	}
	
	if (doCalcFlowNorm(\%sample, \%location) == 1) {
		my $copies_flownorm = calcCopiesFlowNorm($copies_per_lww, \%sample, \%location);
		$table2watch{"result"}->{"$assay_id"}->{"target_copies_flownorm"} = $copies_flownorm;
		
		if (doCalcPop(\%location) == 1) {
			my $copiesfn_per_lcap = calcCopiesPop($copies_flownorm, 100, \%location);
			$table2watch{"result"}->{"$assay_id"}->{"target_copies_fn_per_cap"} = $copiesfn_per_lcap;
		}
	}
	
	# Finally, add any controls. Here we calculate the mean of each control type (PC or NC).
	# Controls are reported as copies per PCR reaction.
	my $assdye  = $assay{"assay_target_fluorophore"};
	my $assdate = $abatch{"assay_date"};
	my $assvol  = $abatch{"assay_reaction_ul"};
	
	if (defined $controls{"$abid"}) {
		if (defined $controls{"$abid"}->{"$assdye"}) {
			foreach my $ctype (keys %{$controls{"$abid"}->{"$assdye"}}) {
				my $cpr = calcScaledMean($controls{"$abid"}->{"$assdye"}->{"$ctype"}, $assvol);
				if ("$ctype" =~ /positive/i) {
					$table2watch{"result"}->{"$assay_id"}->{"pc_copies_per_rxn"} = $cpr;
				} elsif ("$ctype" =~ /negative/i) {
					$table2watch{"result"}->{"$assay_id"}->{"nc_copies_per_rxn"} = $cpr;
					if ($cpr <= $NTC_THRESHOLD) {
						$table2watch{"result"}->{"$assay_id"}->{"target_result_validated"} = "yes";
					} else {
						$table2watch{"result"}->{"$assay_id"}->{"target_result_validated"} = "NTC above threshold";
					}
				} else {
					print "WARN : Unknown type of control '$ctype' attempted for assay batch $abid.\n";
					print "WARN : This is odd, but not fatal so I'm just ignoring this control and moving on.\n";
					print "WARN : You may want to check the controls for $abid in the assay table.\n";
				}
			}
		} else {
			$table2watch{"result"}->{"$assay_id"}->{"target_result_validated"} = "Missing controls: $assdye";
			print $ctlFH "$abid\t$assdye\t$assdate\n";
		}
	} else {
		$table2watch{"result"}->{"$assay_id"}->{"target_result_validated"} = "Missing controls: all";
		print $ctlFH "$abid\tALL\t$assdate\n";
	}
	
}

close $resFH;

#print Dumper($table2watch{"result"});
#die;


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
											 "event_type", 
											 "sample_flow",
											 "sample_qc", 
											 "location_id",
											 "target",
											 "target_genetic_locus", 
											 "lab_id",
											 "target_copies_per_l",
											 "target_copies_per_ld",
											 "target_copies_per_ldcap",
											 "target_copies_flownorm",
											 "target_copies_fn_per_cap",
											 "target_per_capita_basis", 
											 "nc_copies_per_rxn",
											 "pc_copies_per_rxn",
											 "target_result_validated");
											 
open (my $rFH, ">", "$dbdir/watchdb.result.txt") or die "Unable to open $dbdir/watchdb.result.txt for writing: $!";
print $rFH join("\t", @result_colnames) . "\n";
foreach my $assay_id (keys %{$table2watch{"result"}}) {
	print $rFH "$assay_id";
	for (my $i=1; $i < scalar(@result_colnames); $i++) {
		my $colname = $result_colnames[$i];
		my $colval  = "NA";
		if (defined $table2watch{"result"}->{"$assay_id"}->{"$colname"}) {
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

sub doCalcBase {
	my $a   = shift;
	my $ab  = shift;
	my $eb  = shift;
	my $cb  = shift;
	my $loc = shift;
	
	my $do = 1;
	$do = 0 unless defined $a->{"assay_target_copies_per_ul_reaction"} and 
												 "$a->{assay_target_copies_per_ul_reaction}" ne "" and 
												 "$a->{assay_target_copies_per_ul_reaction}" ne "NA"; 

	$do = 0 unless defined $a->{"assay_input_ul"} and 
												 "$a->{assay_input_ul}" ne "" and 
												 "$a->{assay_input_ul}" ne "NA";
												 
	$do = 0 unless defined $ab->{"assay_reaction_ul"} and 
												 "$ab->{assay_reaction_ul}" ne "" and 
												 "$ab->{assay_reaction_ul}" ne "NA";

	$do = 0 unless defined $eb->{"extraction_output_ul"} and 
												 "$eb->{extraction_output_ul}" ne "" and 
												 "$eb->{extraction_output_ul}" ne "NA";
												  
	$do = 0 unless defined $eb->{"extraction_input_ul"} and 
												 "$eb->{extraction_input_ul}" ne "" and 
												 "$eb->{extraction_input_ul}" ne "NA";

	$do = 0 unless defined $cb->{"concentration_output_ml"} and 
												 "$cb->{concentration_output_ml}" ne "" and 
												 "$cb->{concentration_output_ml}" ne "NA";
												  
	$do = 0 unless defined $cb->{"concentration_input_ml"} and 
												 "$cb->{concentration_input_ml}" ne "" and 
												 "$cb->{concentration_input_ml}" ne "NA";

	$do = 0 unless defined $loc->{"location_collection_window_hrs"} and 
												 "$loc->{location_collection_window_hrs}" ne "" and 
												 "$loc->{location_collection_window_hrs}" ne "NA";
												  
	$do = 0 unless defined $loc->{"location_collection_basis"} and 
												 "$loc->{location_collection_basis}" ne "" and 
												 "$loc->{location_collection_basis}" ne "NA";
	return $do;
}

sub doCalcFlowNorm {
	my $s   = shift;
	my $loc = shift;
	
	my $do = 1;
	$do = 0 if !defined $s->{"sample_flow"} or 
												 "$s->{sample_flow}" eq "0" or 
												 "$s->{sample_flow}" eq "" or 
												 "$s->{sample_flow}" eq "NA";

	$do = 0 if !defined $loc->{"location_collection_window_hrs"} or 
												 "$loc->{location_collection_window_hrs}" eq "" or 
												 "$loc->{location_collection_window_hrs}" eq "NA"; 

	$do = 0 if !defined $loc->{"location_collection_basis"} or 
												 "$loc->{location_collection_basis}" eq "" or 
												 "$loc->{location_collection_basis}" eq "NA";

	return $do;
}

sub doCalcPop {
	my $loc = shift;
	
	my $do = 1;
	$do = 0 if !defined $loc->{"location_population_served"} or 
												 "$loc->{location_population_served}" eq "" or 
												 "$loc->{location_population_served}" eq "NA";

	return $do;
}

sub calcCopiesPerL {

	my $a  = shift;
	my $ab = shift;
	my $eb = shift;
	my $cb = shift;
	
	my $val = "NA";
	$val = 
		1000 * 
		1000 * 
		($a->{"assay_target_copies_per_ul_reaction"} * 
		 $ab->{"assay_reaction_ul"} * 
		 $eb->{"extraction_output_ul"} * 
		 $cb->{"concentration_output_ml"}
		) / 
		($a->{"assay_input_ul"} * 
		 $eb->{"extraction_input_ul"} * 
		 $cb->{"concentration_input_ml"});
	
	return $val;
}

sub calcCopiesPop {
	my $cn    = shift;
	my $basis = shift;
	my $loc   = shift;

	my $val = "-";

	my $pop = $loc->{"location_population_served"};
	$pop =~ s/,//i;
	$val = $cn / ($pop/$basis);
	
	return $val;
}

sub calcCopiesFlowNorm {
	# copies/liter x 3.78541 liters/gallon x 10^6 gallon/million_gallons x flow million_gallons/day x 1 day/24 hrs x collection_time_hrs
	my $cpl = shift;
	my $s   = shift;
	my $loc = shift;

	my $val = "NA";
	if ("$loc->{location_collection_basis}" =~ /flow/) {
		$val = $cpl;
	} else {
		my $flow = $s->{"sample_flow"};
		$flow =~ s/,//i;
		$val = $cpl * 3.78541 * 1000000 * $flow * ($loc->{"location_collection_window_hrs"} / 24);
	}
	
	return $val;
}


sub calcScaledMean {
	my $arrRef = shift;
	my $mult   = shift;
	
	return 0 unless scalar @{$arrRef} > 0;
	$mult = 1 unless defined $mult;
	
	my $val = 0;
	foreach (@{$arrRef}) {
		$val += ($_ * $mult);
	}
	$val /= scalar(@{$arrRef});
	
	return $val;
}



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



