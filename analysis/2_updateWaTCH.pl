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

die "FATAL: $progname requires a valid $assetsF file. This is produced by the workflow.sh script.\n$usage\n" unless -f "$rundir/$assetsF";
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
my $FIELDSFILE     = "resources/fields_watchdb.txt";
my $ATDOMFILE      = "resources/fields_watchdb_sample.txt";
my $BIOREPFILE     = "resources/fields_watchdb_biol.txt";
my $TECHREPFILE    = "resources/fields_watchdb_tech.txt";
my $WATCHFILE_MAIN = "updates/watchdb.LATEST.txt";
my $WATCHFILE_INCR = "updates/watchdb/watchdb.$NOW.txt";
my $LOGFILE        = "logs/watchdb/watchdb.$NOW.log";


`cp $WATCHFILE_MAIN $WATCHFILE_MAIN.OLD`;


# Open the log file for writing.
open (my $logFH, ">", "$LOGFILE") or die "Unable to open $LOGFILE for reading: $!\n";
print $logFH "# WaTCHdb update log.\n";
print $logFH "# Command: $progname.\n";
print $logFH "# Data output written to $WATCHFILE_INCR and $WATCHFILE_MAIN.\n";
print $logFH "# Fatal errors and system warnings are NOT logged here; see STDERR for those.\n";
print $logFH "#\n";
print $logFH "# Run started $NOW.\n";
print $logFH "#\n";

print "# Run started $NOW.\n";


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

print $logFH scalar(keys %TARGETS) . " molecular targets added.\n";
print scalar(keys %TARGETS) . " molecular targets added.\n";


#
# Get the master list of WaTCHdb field names. These are ALREADY SORTED.
#
my @watchdb_fields = ();
open (my $watchFieldsFH, "<", "$FIELDSFILE") or die "Unable to open $FIELDSFILE for reading: $!";
while (<$watchFieldsFH>) {
	chomp;
	push @watchdb_fields, "$_";
}
close $watchFieldsFH;


#
# Get the list of WaTCHdb field names that AT is allowed to overwrite (sample metadata).
#
my %atdom_fields = ();
open (my $atdomFH, "<", "$ATDOMFILE") or die "Unable to open $ATDOMFILE for reading: $!";
while (<$atdomFH>) {
	chomp;
	$atdom_fields{"$_"} = 1;
}
close $atdomFH;



#
# Get the list of WaTCHdb field names associated with biological reps.
#
my %biorep_fields = ();
open (my $biorepFH, "<", "$BIOREPFILE") or die "Unable to open $BIOREPFILE for reading: $!";
while (<$biorepFH>) {
	chomp;
	$biorep_fields{"$_"} = 1;
}
close $biorepFH;



#
# Get the list of WaTCHdb field names associated with technical reps (assay).
#
my %techrep_fields = ();
open (my $techrepFH, "<", "$TECHREPFILE") or die "Unable to open $TECHREPFILE for reading: $!";
while (<$techrepFH>) {
	chomp;
	$techrep_fields{"$_"} = 1;
}
close $techrepFH;



#
# Fetch input from latest WaTCH database.
# Use this to populate the keys of the master hash %asset2data
# Skip if there is no database file
#
my @w_field_names = ();
my $count = 0;

if (-f "$WATCHFILE_MAIN") {
	open (my $watchInFH, "<", "$WATCHFILE_MAIN") or die "Unable to open $WATCHFILE_MAIN for reading: $!\n";
	while (my $line = <$watchInFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @values = split "\t", "$line", -1;

	if ($count == 0) {
		# populate the keys array with values on the first line
		for (my $i=0; $i < scalar(@values); $i++) {
			$w_field_names[$i] = $values[$i];
		}
	} else {
			# first field is the asset ID. second is the replicate id.
			my $asset_id = $values[0];
			my $rep_id = $values[1];
			
			# populate the hash entry for this asset with the remaining values
			$asset2data{"$asset_id"} = {} unless defined $asset2data{"$asset_id"};
			$asset2data{"$asset_id"}->{"$rep_id"} = {} unless defined $asset2data{"$asset_id"}->{"$rep_id"};
			for (my $i=2; $i < scalar(@values); $i++) {
				$asset2data{"$asset_id"}->{"$rep_id"}->{"$w_field_names[$i]"} = $values[$i];
			}
		}
		$count++;
	}
	close $watchInFH;
}
print $logFH "INFO: Finished WaTCHdb import. " . scalar(keys(%asset2data)) . " total assets present.\n";
print "INFO: Finished WaTCHdb import. " . scalar(keys(%asset2data)) . " total assets present.\n";

#print Dumper(\%asset2data);
#die;


#
# Read in assets from AssetTiger file.
# This will populate IDs and metadata for samples in the current run and any other samples not in the master db yet.
# This data goes into the asset2data hash.
#
# AssetTiger field metadata overwrites WaTCH data by default, to allow field team to update info. These names are keys of %atdom_fields.
#
# WaTCH may contain multiple replicates on a single asset id, while AT contains only the latest replicate. So if an incoming AT 
# id is already present in WaTCH, overwrite with the metadata (%atdom_fields) but otherwise, ignore the AT entry.
#
#
my @a_field_names = ();
my $new_from_AT = 0;

my $csvA = Text::CSV->new({auto_diag => 4, binary => 1});
open (my $assetsInFH, "<", "$rundir/$assetsF") or die "Unable to open $rundir/$assetsF for reading: $!\n";
$count = 0;
while (my $line = $csvA->getline($assetsInFH)) {
	my @values = @$line;
	@values = @{trim(\@values)};
	
	if ($count == 0) {
		# populate the fields array with values on the first line
		for (my $i=0; $i < scalar(@values); $i++) {
			$a_field_names[$i] = $values[$i];
		}
	} else {
		# value in the first field is the asset ID / sample ID
		my $asset_id = $values[0];

		if (defined $asset2data{"$asset_id"}) {
			# incoming asset from AT found in WaTCHdb
			# overwrite WaTCHdb metadata from AT (%atdom_fields membership required) but otherwise ignore
			foreach my $rep_id (keys %{$asset2data{"$asset_id"}}) {
				for (my $i=1; $i < scalar(@values); $i++) {
					next if !defined $values[$i] or "$values[$i]" eq "" or "$values[$i]" eq "NA" or "$values[$i]" eq "-";
					$asset2data{"$asset_id"}->{"$rep_id"}->{"$a_field_names[$i]"} = "$values[$i]" if defined $atdom_fields{"$a_field_names[$i]"};
				}
			}
		} else {
			# incoming asset from AT is absent from WaTCHdb
			# add as a new id to WatchDB and write AT data
			$new_from_AT++;
			$asset2data{"$asset_id"} = {};
			# Since this is a new db entry, replicate id is 1.1 by default
			$asset2data{"$asset_id"}->{"1.1"} = {};
			
			# init new entry with WaTCHdb fields
			for (my $i=2; $i < scalar(@w_field_names); $i++) {
				$asset2data{"$asset_id"}->{"1.1"}->{"$w_field_names[$i]"} = "";
			}
			
			# add AT values for this new entry. Should only be sample info at this point.
			for (my $i=2; $i < scalar(@values); $i++) {
				# Set the Status: NWSS and Status: WaTCH of new entries to "Processing"
				if ("Status: NWSS" eq "$a_field_names[$i]" or "Status: WaTCH" eq "$a_field_names[$i]") {
					$asset2data{"$asset_id"}->{"1.1"}->{"$a_field_names[$i]"} = "Processing" if !defined $asset2data{"$asset_id"}->{"1.1"}->{"$a_field_names[$i]"} or $asset2data{"$asset_id"}->{"1.1"}->{"$a_field_names[$i]"} eq "" or "values[$i]" eq "";
				} else {
					next if !defined $values[$i] or "$values[$i]" eq "" or "$values[$i]" eq "NA" or "$values[$i]" eq "-";
					$asset2data{"$asset_id"}->{"1.1"}->{"$a_field_names[$i]"} = "$values[$i]" if defined $atdom_fields{"$a_field_names[$i]"};
				}
			}
		}
	}
	$count++;
}
close $assetsInFH;

print $logFH "INFO: Finished AssetTiger import. $new_from_AT new assets added. " . scalar(keys(%asset2data)) . " total assets present.\n";
print "INFO: Finished AssetTiger import. $new_from_AT new assets added. " . scalar(keys(%asset2data)) . " total assets present.\n";

#print Dumper(\%asset2data);
#die;


#
# Read the run metadata file, which contains all of the merged plate file information
# This goes into the master asset2data hash as well
# Also parses out a map of well id to sample id for later use
#
my @r_field_names = ();
my %sample2well   = ();

my $csvR = Text::CSV->new({auto_diag => 4, binary => 1});
open (my $runMetaFH, "<", "$rundir/$runMetaF") or die "Unable to open $rundir/$runMetaF for reading: $!\n";
$count = 0;
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
			if (defined $asset2data{"$asset_id"}) {
				# Add assay data to the entry.
				# If a complete rep already exists, add as a new rep. Completed assay requires an Assay Date.
				my ($brep, $trep) = (1,1);
				
				while (defined $asset2data{"$asset_id"}->{"${brep}.${trep}"} and defined $asset2data{"$asset_id"}->{"${brep}.${trep}"}->{"Assay Date"} and $asset2data{"$asset_id"}->{"${brep}.${trep}"}->{"Assay Date"} ne "") {
					$trep++;
				}
				
				# initialize the rep if it is being added here
				unless (defined $asset2data{"$asset_id"}->{"${brep}.${trep}"}) {
					$asset2data{"$asset_id"}->{"${brep}.${trep}"} = {};
					for (my $i=2; $i < scalar(@w_field_names); $i++) {
						$asset2data{"$asset_id"}->{"${brep}.${trep}"}->{"$w_field_names[$i]"} = "";
					}
					# copy over all the data from rep 1.1 except Assay data (%techrep_fields).
					for (my $i=2; $i < scalar(@w_field_names); $i++) {
						my $key = "$w_field_names[$i]";
						next if defined $techrep_fields{"$key"};
						$asset2data{"$asset_id"}->{"${brep}.${trep}"}->{"$key"} = $asset2data{"$asset_id"}->{"1.1"}->{"$key"};
					}
					
				}
				
				# set run metadata
				for (my $i=1; $i < scalar(@values); $i++) {
					next if "$values[$i]" eq "";
					$asset2data{"$asset_id"}->{"${brep}.${trep}"}->{"$r_field_names[$i]"} = "$values[$i]";
					if ("$r_field_names[$i]" eq "Assay Plate Location") {
						$sample2well{"$asset_id"} = {} unless defined $sample2well{"$asset_id"};
						$sample2well{"$asset_id"}->{"${brep}.${trep}"} = "$values[$i]";
					}
				}
				
			} else {
				print $logFH "HIGH: Mysterious sample ID $asset_id found during Run Metadata import, does not exist in AT or WaTCHdb. This is concerning.\nHIGH: $asset_id will be ignored!\n";
				print "HIGH: Mysterious sample ID $asset_id found during Run Metadata import, does not exist in AT or WaTCHdb. This is concerning.\nHIGH: $asset_id will be ignored!\n";
			}
			
		}

	}
	$count++;
}
close $runMetaFH;


#
# At this point, the asset2data hash should contain all data for all samples except samples assayed in the current run
# (and possibly future runs if they've been logged into AssetTiger already)
#
#print $logFH "INFO: Finished run metadata import. " . scalar(keys(%asset2data)) . " total assets present.\n";
#print "INFO: Finished run metadata import. " . scalar(keys(%asset2data)) . " total assets present.\n";

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

	foreach my $rep_id (keys %{$sample2well{"$asset_id"}}) {
		# get the well id for this sample/asset
		my $well_id = $sample2well{"$asset_id"}->{"$rep_id"};
		next unless defined $well_id and $well_id ne "";
		#print "$asset_id\t$well_id\n";
	
		# All of the assay results for this sample/asset should be stored under the corresponding well id in the well2assay hash
		# Transfer it over into the master asset2data hash, which is keyed by the sample/asset id
		#
		unless (defined $well2assay{$well_id}) {
			$dowrite = 0;
			print $logFH "HIGH: During run data merge, well ID $well_id is not in the well2assay hash. This is fatal.\n";
			print $logFH "    Here are the well2assay hash keys: " . join(", ", keys(%well2assay)) . "\n";
			my $END =
				 DateTime
						->now( time_zone => 'local' )
						->set_time_zone('floating')
						->strftime('%Y-%m-%d.%H-%M-%S');

			print $logFH "#\n";
			print $logFH "# WaTCHdb update ABORTED.\n";
			print $logFH "# $END.\n";
			close $logFH;
			print "HIGH: During run data merge, well ID $well_id is not in the well2assay hash. This is fatal.\n";
			print "    Here are the well2assay hash keys: " . join(", ", keys(%well2assay)) . "\n";
			print "# WaTCHdb update ABORTED.\n";
			print "# $END.\n";
			die;
		}
	
		my %sample_data = %{$well2assay{$well_id}};
		$asset2data{"$asset_id"}->{"$rep_id"}->{"Assay Droplet Quantification Method"} = $sample_data{"Status"};
		$asset2data{"$asset_id"}->{"$rep_id"}->{"Assay Experiment Type"} = $sample_data{"Experiment"};
		$asset2data{"$asset_id"}->{"$rep_id"}->{"Assay Accepted Droplets"} = $sample_data{"Accepted Droplets"};
		
		# Calculate the CN/L for each target in the current run found in the TARGETS hash
		# If the concentration is negative it means the assay returned a "No Call"
		# This should be replaced by a ND or something equivalent, but at the moment it is just given as 0
		#
		foreach my $target_id (keys %TARGETS) {
			my $target_name = $TARGETS{$target_id};
			if (defined $sample_data{"$target_name"}) {
				#print "$asset_id\t$target_name\t";
				my $assay_result = 0;
				$assay_result = calc_result("$asset_id", "$rep_id", $sample_data{"$target_name"}->{"Concentration (CN/Rxn)"}) unless $sample_data{"$target_name"}->{"Concentration (CN/Rxn)"} < 0;
				$asset2data{"$asset_id"}->{"$rep_id"}->{"$target_id"} = "$target_name";
				$asset2data{"$asset_id"}->{"$rep_id"}->{"$target_id Result (CN/L)"} = $assay_result;
				$asset2data{"$asset_id"}->{"$rep_id"}->{"$target_id Concentration (CN/Rxn)"} = $sample_data{"$target_name"}->{"Concentration (CN/Rxn)"};
			
				#
				# Record results of the PCR controls for the current run (keys in %current_ctl)
				# Each sample in the current run will receive the PCR control data
				#
				my $well_id_nc = $ctl2data{"PCR NC"}->{"Assay Plate Location"};
				my %nc_data = %{$well2assay{$well_id_nc}};
				$asset2data{"$asset_id"}->{"$rep_id"}->{"$target_id PCR NC Result (CN/Rxn)"} = $nc_data{"$target_name"}->{"Concentration (CN/Rxn)"};
				my $well_id_pc = $ctl2data{"PCR PC"}->{"Assay Plate Location"};
				my %pc_data = %{$well2assay{$well_id_pc}};
				$asset2data{"$asset_id"}->{"$rep_id"}->{"$target_id PCR PC Result (CN/Rxn)"} = $pc_data{"$target_name"}->{"Concentration (CN/Rxn)"};
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
=cut

#print Dumper(\%asset2data);
#die;



#
# %asset2data now contains all the asset data
# Write to the WaTCHdb file
#
if ($dowrite == 0) {
	print "\n###\n##\n#\nThere was at least one fatal or high-level warning during the update. As a result, $WATCHFILE_MAIN has NOT been changed!\n#\n##\n###\n\n"; 
	`rm $WATCHFILE_MAIN.OLD`;
} else {
	open (my $watchFH, ">", "$WATCHFILE_INCR") or die "Unable to open $WATCHFILE_INCR for writing: $!";
	print $watchFH join("\t", @watchdb_fields) . "\n";
	foreach my $asset_id (keys %asset2data) {
		foreach my $rep_id (keys %{$asset2data{"$asset_id"}}) {
			my $asset = $asset2data{"$asset_id"}->{"$rep_id"};
			print $logFH "LOW : WaTCH status of $asset_id $rep_id has been changed to " . $asset->{"Status: WaTCH"} . ".\n" if change_status($asset) == 1;
			print $watchFH "$asset_id\t$rep_id";
			foreach my $field (@watchdb_fields) {
				next if "$field" eq "Sample ID" or "$field" eq "Asset ID" or "$field" eq "Replicate (Bio.Tech)";
				my $value = "";
				$value = "$asset->{$field}" if defined $asset->{"$field"};
				print $watchFH "\t$value";
			}
			print $watchFH "\n";
		}
	}
	close $watchFH;

	`cp $WATCHFILE_INCR $WATCHFILE_MAIN`;
}


my $END =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');


print $logFH "#\n";
print $logFH "# WaTCHdb update complete.\n";
print $logFH "# $END.\n";
close $logFH;

print "INFO: WaTCHdb update complete.\n";
print "INFO: Run finished $END.\n";


exit 0;

=cut
sub do_replace {
	my $asset   = shift;
	my $a_field = shift;
	my $newval  = shift;
		
	# Test incoming new value from AT for existence. Do not replace if it is not a valid value.
	return 0 unless defined $newval and "$newval" ne "" and "$newval" ne "NaN" and "$newval" ne "-";

	# Test existing asset value for existence. Can't replace if it doesn't exist - indicates a spurious field was exported from AT and can be ignored.
	return 0 unless defined $asset->{"$a_field"};
	
	# If the new value is identical to the existing value, there is no need to replace.
	return 0 if "$newval" eq "$asset->{$a_field}";
	
	# Ignore dates; AT, Excel, and perl format date strings differently.
#	return 0 if "$a_field" =~ /date/i or "$a_field" =~ /start/i or "$a_field" =~ /end/i;
	
	# If the value is a Result or a Status, NEVER replace WaTCHdb data with AT data.
	return 0 if "$a_field" =~ /result/i or "$a_field" =~ /status:/i;
	
	# If we reach this line, the new and existing values are valid, not dates, and different from each other.
	# In this case, we want to replace the WaTCHdb value with the new value from AT.
	return 1;
}
=cut


sub calc_result {
	my $aid = shift;
	my $rid = shift;
	my $CN = shift;
	
	my %calc_vals = ();
	$calc_vals{"Extraction Output Volume (mL)"} = $asset2data{"$aid"}->{"$rid"}->{"Extraction Output Volume (mL)"};
	$calc_vals{"Concentration Output Volume (mL)"} = $asset2data{"$aid"}->{"$rid"}->{"Concentration Output Volume (mL)"};
	$calc_vals{"Assay Volume (mL)"} = $asset2data{"$aid"}->{"$rid"}->{"Assay Volume (mL)"};
	$calc_vals{"Extraction Input Volume (mL)"} = $asset2data{"$aid"}->{"$rid"}->{"Extraction Input Volume (mL)"};
	$calc_vals{"Concentration Input Volume (mL)"} = $asset2data{"$aid"}->{"$rid"}->{"Concentration Input Volume (mL)"};
	
	my $valid = 1;
	foreach my $k (keys %calc_vals) {
		if (!defined $calc_vals{$k} or $calc_vals{$k} eq "") {
			print $logFH "HIGH: I have assay data for $aid $rid in this run, but no value for $k. This is concerning.\n";
			print "HIGH: I have assay data for $aid $rid in this run, but no value for $k. This is concerning.\n";
			$dowrite = 0;
			$valid = 0;
		} elsif ($k eq "Assay Volume (mL)" or $k eq "Extraction Input Volume (mL)" or $k eq "Concentration Input Volume (mL)") {
			if ($calc_vals{$k} == 0) {
				print $logFH "HIGH: Sample $aid $rid in this run has a $k value of 0. This is concerning.\n";
				print "HIGH: Sample $aid $rid in this run has a $k value of 0. This is concerning.\n";
				$dowrite = 0;
				$valid = 0;
			}
		}
	}
	unless ($valid == 1) {
		print $logFH "    I am unable to calculate a result for $aid $rid. This should be confirmed.\n";
		print "    I am unable to calculate a result for $aid $rid. This should be confirmed.\n";
		$dowrite = 0;
		return -1;
	}
	
	my $cn_per_l = 1000 * 
								 ($CN * $calc_vals{"Extraction Output Volume (mL)"} * $calc_vals{"Concentration Output Volume (mL)"}) / 
								 ($calc_vals{"Assay Volume (mL)"} * $calc_vals{"Extraction Input Volume (mL)"} * $calc_vals{"Concentration Input Volume (mL)"});
	print $logFH "INFO: Sample $aid $rid has a new calculated value of $cn_per_l.\n";
	#print "INFO: Sample $aid has a calculated value of $cn_per_l.\n";
	return $cn_per_l;
	
	#return 0;
}


sub change_status {
	my $asset_hash = shift;
	
	# If NWSS status is empty, set to the default of Processing. This will be checked during NWSS output (separate script).
	$asset_hash->{"Status: NWSS"} = "Processing" unless defined $asset_hash->{"Status: NWSS"} and $asset_hash->{"Status: NWSS"} ne "";
	# If WaTCH status is empty, set to the default of Processing.
	$asset_hash->{"Status: WaTCH"} = "Processing" unless defined $asset_hash->{"Status: WaTCH"} and $asset_hash->{"Status: WaTCH"} ne "";
	
	# Do nothing if the WaTCH status is already Complete or Rejected. These are EOL status settings.
	return 0 if $asset_hash->{"Status: WaTCH"} eq "Complete" or $asset_hash->{"Status: WaTCH"} =~ /Rejected/;
	
	# Only if there is an existing droplet count, the status can change to either Confirming or Complete.
	if (defined $asset_hash->{"Assay Accepted Droplets"} and $asset_hash->{"Assay Accepted Droplets"} ne "" and $asset_hash->{"Assay Accepted Droplets"} >= 8000) {
		if (defined $asset_hash->{"Assay Target 1 Result (CN/L)"} and $asset_hash->{"Assay Target 1 Result (CN/L)"} ne "" and $asset_hash->{"Assay Target 1 Result (CN/L)"} ne "-1") {
			$asset_hash->{"Status: WaTCH"} = "Complete";
		} else {
			$asset_hash->{"Status: WaTCH"} = "Confirming";
		}
		return 1;
	}
	
	return 0;
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







