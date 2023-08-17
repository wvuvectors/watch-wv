#! /usr/bin/env perl


use strict;
use warnings;
use Text::CSV;
use DateTime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use Data::Dumper;


# read in plate files from directory and convert to single table

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=  "Build update table for WaTCH-WV database.\n";
$usage   .=   "\n";

my $rundir;

my ($assetsF, $runMetaF, $assayDataF) = ("assets_all.csv", "run_metadata.csv", "assay_plate.csv");

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } else {
		$rundir = $arg;
	}
}

die "FATAL: $progname requires a valid run directory.\n$usage\n" unless defined $rundir and -d "$rundir";

die "FATAL: $progname requires a valid $runMetaF file. This is produced by the workflow.sh script.\n$usage\n" unless -f "$rundir/$runMetaF";
die "FATAL: $progname requires a valid $assayDataF file. This is produced by the workflow.sh script.\n$usage\n" unless -f "$rundir/$assayDataF";


# keyed master hash for all compiled input data
my %asset2data = ();
my %ctl2data   = ();

# keyed hash for controls in the current run
my %current_ctl = ();

# list of targets and their database field keys
my %TARGETS = ();
my $REACTION_VOL_UL = 20;

#
# Get time stamp for now
#
my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');



# Resource and DB filepaths
my $TARGETSFILE    = "resources/molecular_targets.txt";
my $FIELDSFILE     = "resources/fields_bsajam2023.txt";

my $TECHREPFILE    = "resources/fields_watchdb_tech.txt";
my $WATCHFILE_MAIN = "updates/bsajam2023.LATEST.txt";
my $WATCHFILE_INCR = "updates/bsajam2023/bsajam2023.$NOW.txt";

#`cp $WATCHFILE_MAIN $WATCHFILE_MAIN.OLD`;



# Get the master list of targets
# target id corresponds to the column name in the assay data file
# target name corresponds to the value for that field
#
open (my $tFH, "<", "$TARGETSFILE") or die "Unable to open $TARGETSFILE for reading: $!\n";
while (my $line = <$tFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my ($target_id, $target_name) = split /\t/, "$line", -1;
	$TARGETS{$target_id} = $target_name;
}
close $tFH;

print scalar(keys %TARGETS) . " molecular targets added.\n";


#
# Get the master list of WaTCHdb field names for BSa Jam 2023. These are ALREADY SORTED.
#
my @watchdb_fields = ();
open (my $watchFieldsFH, "<", "$FIELDSFILE") or die "Unable to open $FIELDSFILE for reading: $!";
while (<$watchFieldsFH>) {
	chomp;
	push @watchdb_fields, "$_";
}
close $watchFieldsFH;


#
# Read the run metadata file, which contains all of the merged plate file information
# This goes into the master asset2data hash
# Also parses out a map of well id to sample id for later use
#
my @r_field_names = ();
my %sample2well   = ();

my $csvR = Text::CSV->new({auto_diag => 4, binary => 1});
open (my $runMetaFH, "<", "$rundir/$runMetaF") or die "Unable to open $rundir/$runMetaF for reading: $!\n";
my $count = 0;
while (my $line = $csvR->getline($runMetaFH)) {
	my @values = @$line;

	if ($count == 0) {
		# populate the keys array with values on the first line
		for (my $i=0; $i < scalar(@values); $i++) {
			my $key = $values[$i];
#			$key = "Control, Internal (Type)" if "$key" eq "Assay Internal Control";
#			$key = "Control, Process (Type)" if "$key" eq "Assay Process Control";
			$r_field_names[$i] = "$key";
		}
	} else {
		# first field is the sample ID == asset ID
		my $asset_id = $values[0];

		# PCR control data get stored in a separate hash; they are processed differently from samples
		if ("$asset_id" =~ /^PCR/) {
			$current_ctl{$asset_id} = 1;
			$ctl2data{$asset_id} = {} unless defined $ctl2data{$asset_id};
			for (my $i=1; $i < scalar(@values); $i++) {
				$ctl2data{"$asset_id"}->{"$r_field_names[$i]"} = "$values[$i]" unless "$values[$i]" eq "";
			}
		} else {
			# populate/update the master hash entry for this asset with the sample metadata
				# Add assay data to the entry.
				$asset2data{"$asset_id"} = {};
				
				# set run metadata
				for (my $i=1; $i < scalar(@values); $i++) {
					next if "$values[$i]" eq "";
					$asset2data{"$asset_id"}->{"$r_field_names[$i]"} = "$values[$i]";
					if ("$r_field_names[$i]" eq "Assay Plate Location") {
						$sample2well{"$asset_id"} = {} unless defined $sample2well{"$asset_id"};
						$sample2well{"$asset_id"} = "$values[$i]";
					}
				}
		}

	}
	$count++;
}
close $runMetaFH;

#print Dumper(\%asset2data);
#die;


#
# Read and parse the ddPCR data file for the current run
# These data are keyed by plate location, not sample id, so they are stored in a separate hash, %well2assay
#
my %well2assay = ();
my %p_field_pos = ("Target" => -1, "Status" => -1, "Experiment" => -1, "DyeName(s)" => -1, "Conc(copies/uL)" => -1, "Accepted Droplets" => -1);

my $csvP = Text::CSV->new({auto_diag => 4, binary => 1});
open (my $assayDataFH, "<", "$rundir/$assayDataF") or die "Unable to open $rundir/$assayDataF for reading: $!\n";
$count = 0;
while (my $line = $csvP->getline($assayDataFH)) {
	my @values = @$line;
	next if scalar @values < 1;
	
	if ($count == 0) {
		# populate the p_field_pos hash with values on the first line (as keys) and their column positions (as values)
		for (my $i=0; $i < scalar(@values); $i++) {
			$p_field_pos{"$values[$i]"} = $i if defined $p_field_pos{"$values[$i]"};
		}
	} else {
		# first column is the well ID
		# there will be multiple rows for each well in a multiplex file
		my $well_id = $values[0];
		
		# populate this well2assay hash entry for this well using the remaining column values
		$well2assay{$well_id} = {} unless defined $well2assay{$well_id};
		$well2assay{$well_id}->{"Status"} = $values[$p_field_pos{"Status"}];
		$well2assay{$well_id}->{"Experiment"} = $values[$p_field_pos{"Experiment"}];
		$well2assay{$well_id}->{"Accepted Droplets"} = $values[$p_field_pos{"Accepted Droplets"}];
		
		# sort out each target in the well
		my $target_name = $values[$p_field_pos{"Target"}];
		$well2assay{$well_id}->{"$target_name"} = {} unless defined $well2assay{$well_id}->{"$target_name"};
		$well2assay{$well_id}->{"$target_name"}->{"DyeName(s)"} = $values[$p_field_pos{"DyeName(s)"}];
		
		# A negative value here indicates an assay failure
		my $cn_per_ul = $values[$p_field_pos{"Conc(copies/uL)"}];
		#print "$well_id\t$cn_per_ul\n";
		$cn_per_ul = -1 if $cn_per_ul eq "" or $cn_per_ul eq "No Call";
		$well2assay{$well_id}->{"$target_name"}->{"Concentration (CN/Rxn)"} = $cn_per_ul * $REACTION_VOL_UL;
		
		#foreach my $key (keys %p_field_pos) {
		#	$well2assay{"$well_id"}->{"$values[$p_field_pos{Target}]"}->{"$key"} = $values[$p_field_pos{"$key"}] unless "$key" eq "Target";
		#}
		
	}
	$count++;
}
close $assayDataFH;

#print Dumper(\%well2assay);
#die;


my $dowrite = 1;

#
# merge data from current run - stored in %well2assay and keyed by well id - into the main asset2data hash - keyed by asset id
# also calculate the CN/L for each target in the current run
#
foreach my $asset_id (keys %sample2well) {

	# get the well id for this sample/asset
	my $well_id = $sample2well{"$asset_id"};
	next unless defined $well_id and $well_id ne "";
	#print "$asset_id\t$well_id\n";

	# All of the assay results for this sample/asset should be stored under the corresponding well id in the well2assay hash
	# Transfer it over into the master asset2data hash, which is keyed by the sample/asset id
	#
	unless (defined $well2assay{$well_id}) {
		$dowrite = 0;
		my $END =
			 DateTime
					->now( time_zone => 'local' )
					->set_time_zone('floating')
					->strftime('%Y-%m-%d.%H-%M-%S');

		print "HIGH: During run data merge, well ID $well_id is not in the well2assay hash. This is fatal.\n";
		print "    Here are the well2assay hash keys: " . join(", ", keys(%well2assay)) . "\n";
		print "# WaTCHdb update ABORTED.\n";
		print "# $END.\n";
		die;
	}

	my %sample_data = %{$well2assay{$well_id}};
	$asset2data{"$asset_id"}->{"Assay Droplet Quantification Method"} = $sample_data{"Status"};
	$asset2data{"$asset_id"}->{"Assay Experiment Type"} = $sample_data{"Experiment"};
	$asset2data{"$asset_id"}->{"Assay Accepted Droplets"} = $sample_data{"Accepted Droplets"};
	
	# Calculate the CN/L for each target in the current run found in the TARGETS hash
	# If the concentration is negative it means the assay returned a "No Call"
	# This should be replaced by a ND or something equivalent, but at the moment it is just given as 0
	#
	foreach my $target_id (keys %TARGETS) {
		my $target_name = $TARGETS{$target_id};
		if (defined $sample_data{"$target_name"}) {
			#print "$asset_id\t$target_name\t";
			my $assay_result = 0;
			$assay_result = calc_result("$asset_id", $sample_data{"$target_name"}->{"Concentration (CN/Rxn)"}) unless $sample_data{"$target_name"}->{"Concentration (CN/Rxn)"} < 0;
			$asset2data{"$asset_id"}->{"$target_id"} = "$target_name";
			$asset2data{"$asset_id"}->{"$target_id Result (CN/L)"} = $assay_result;
			$asset2data{"$asset_id"}->{"$target_id Concentration (CN/Rxn)"} = $sample_data{"$target_name"}->{"Concentration (CN/Rxn)"};
		
			#
			# Record results of the PCR controls for the current run (keys in %current_ctl)
			# Each sample in the current run will receive the PCR control data
			#
			if (defined $ctl2data{"PCR NC"} and defined $ctl2data{"PCR NC"}->{"Assay Plate Location"} and defined $well2assay{$ctl2data{"PCR NC"}->{"Assay Plate Location"}}) {
				my $well_id_nc = $ctl2data{"PCR NC"}->{"Assay Plate Location"};
				my %nc_data = %{$well2assay{$well_id_nc}};
				$asset2data{"$asset_id"}->{"$target_id PCR NC Result (CN/Rxn)"} = $nc_data{"$target_name"}->{"Concentration (CN/Rxn)"};
			} else {
				$asset2data{"$asset_id"}->{"$target_id PCR NC Result (CN/Rxn)"} = "NA";
			}

			if (defined $ctl2data{"PCR PC"} and defined $ctl2data{"PCR PC"}->{"Assay Plate Location"}  and defined $well2assay{$ctl2data{"PCR PC"}->{"Assay Plate Location"}}) {
				my $well_id_pc = $ctl2data{"PCR PC"}->{"Assay Plate Location"};
				my %pc_data = %{$well2assay{$well_id_pc}};
				$asset2data{"$asset_id"}->{"$target_id PCR PC Result (CN/Rxn)"} = $pc_data{"$target_name"}->{"Concentration (CN/Rxn)"};
			} else {
				$asset2data{"$asset_id"}->{"$target_id PCR PC Result (CN/Rxn)"} = "NA";
			}
			
		}
	}
}

=cut
	# store the PCR control result for each asset in the current run
	foreach my $target_id (keys %TARGETS) {
		my $target_name = $TARGETS{$target_id};
		next if $target_id =~ /Control/;
		if (defined $ctl2data{"$target"}) {
			my $assay_result = $assay_data{"$target"}->{"Concentration (CN/Rxn)"};
			
			# validate PCR controls using threshold values
			if ("$target" =~ /PCR PC/ and $current_validation{"PCR PC"}->{"Status"} eq "Pass" and $assay_result < $ctl_thresholds{"PCR PC"}) {
				$current_validation{"$ctl_id"}->{"Status"} = "Fail";
				$current_validation{"$ctl_id"}->{"Reason"} = "PCR Positive Below Threshold";
			} elsif ("$target" =~ /PCR NC/ and $current_validation{"PCR NC"}->{"Status"} eq "Pass" and $assay_result > $ctl_thresholds{"PCR NC"}) {
				$current_validation{"$ctl_id"}->{"Status"} = "Fail";
				$current_validation{"$ctl_id"}->{"Reason"} = "PCR Negative Above Threshold";
			}
			
			foreach my $asset_id (keys %current_run) {
				#$asset2data{$ctl_id}->{"$TARGETS{$target}"} = "$target";
				#$asset2data{$ctl_id}->{"$TARGETS{$target} Result (Copies/Reaction)"} = $assay_result;
				$asset2data{$asset_id}->{"$TARGETS{$target} $ctl_id Result (CN/Rxn)"} = $assay_result;
				$asset2data{$asset_id}->{"QC Check, Assay Positive Control"} = $current_validation{"PCR PC"}->{"Status"};
				$asset2data{$asset_id}->{"QC Check, Assay Positive Control"} = $asset2data{$asset_id}->{"QC Check, Assay Positive Control"} . ": " . $current_validation{"PCR PC"}->{"Reason"} unless $current_validation{"PCR PC"}->{"Status"} eq "Pass";
				$asset2data{$asset_id}->{"QC Check, Assay Negative Control"} = $current_validation{"PCR NC"}->{"Status"};
				$asset2data{$asset_id}->{"QC Check, Assay Negative Control"} = $asset2data{$asset_id}->{"QC Check, Assay Negative Control"} . ": " . $current_validation{"PCR NC"}->{"Reason"} unless $current_validation{"PCR NC"}->{"Status"} eq "Pass";
			}
		}
	}
	
	# remove the control from the %asset2data hash
	delete $asset2data{$ctl_id};
}

#print Dumper(\%asset2data);
#die;

=cut


#
# %asset2data now contains all the asset data
# Write to WATCHFILE_INCR

if ($dowrite == 0) {
	print "\n###\n##\n#\nThere was at least one fatal or high-level warning during the update. As a result, $WATCHFILE_MAIN has NOT been changed!\n#\n##\n###\n\n"; 
	#`rm $WATCHFILE_MAIN.OLD`;
} else {
	open (my $watchFH, ">", "$WATCHFILE_INCR") or die "Unable to open $WATCHFILE_INCR for writing: $!";
	print $watchFH join("\t", @watchdb_fields) . "\n";
	foreach my $asset_id (keys %asset2data) {
		my $asset = $asset2data{"$asset_id"};
		print $watchFH "$asset_id";
		foreach my $field (@watchdb_fields) {
			next if "$field" eq "Sample ID" or "$field" eq "Asset ID" or "$field" eq "Replicate (Bio.Tech)";
			my $value = "";
			$value = "$asset->{$field}" if defined $asset->{"$field"};
			print $watchFH "\t$value";
		}
		print $watchFH "\n";
	}
	close $watchFH;

	`cp $WATCHFILE_INCR $WATCHFILE_MAIN`;
}


exit 0;

sub calc_result {
	my $aid = shift;
	my $CN = shift;
	
	my %calc_vals = ();
	$calc_vals{"Extraction Output Volume (mL)"} = 0.075;
	$calc_vals{"Concentration Output Volume (mL)"} = 0.475;
	$calc_vals{"Assay Volume (mL)"} = 0.005;
	$calc_vals{"Extraction Input Volume (mL)"} = .400;
	$calc_vals{"Concentration Input Volume (mL)"} = 10.0;
	
	my $valid = 1;
	foreach my $k (keys %calc_vals) {
		if (!defined $calc_vals{$k} or $calc_vals{$k} eq "") {
			print "HIGH: I have assay data for $aid in this run, but no value for $k. This is concerning.\n";
			$dowrite = 0;
			$valid = 0;
		} elsif ($k eq "Assay Volume (mL)" or $k eq "Extraction Input Volume (mL)" or $k eq "Concentration Input Volume (mL)") {
			if ($calc_vals{$k} == 0) {
				print "HIGH: Sample $aid in this run has a $k value of 0. This is concerning.\n";
				$dowrite = 0;
				$valid = 0;
			}
		}
	}
	unless ($valid == 1) {
		print "    I am unable to calculate a result for $aid. This should be confirmed.\n";
		$dowrite = 0;
		return -1;
	}
	
	my $cn_per_l = 1000 * 
								 ($CN * $calc_vals{"Extraction Output Volume (mL)"} * $calc_vals{"Concentration Output Volume (mL)"}) / 
								 ($calc_vals{"Assay Volume (mL)"} * $calc_vals{"Extraction Input Volume (mL)"} * $calc_vals{"Concentration Input Volume (mL)"});
	#print "INFO: Sample $aid has a calculated value of $cn_per_l.\n";
	return $cn_per_l;
	
	#return 0;
}


sub trim {
	my $input = shift;
	
	if (!defined ref($input)) {
		$input =~ s/^\s*//gi;
		$input =~ s/\s*$//gi;
	} elsif (ref($input) eq "ARRAY") {
		for (my $i=0; $i < scalar(@$input); $i++) {
			$input->[$i] =~ s/^\s*//gi;
			$input->[$i] =~ s/\s*$//gi;
		}
	} elsif (ref($input) eq "HASH") {
		foreach my $key (keys %$input) {
			$input->{$key} =~ s/^\s*//gi;
			$input->{$key} =~ s/\s*$//gi;
		}
	}

	return $input;
}







