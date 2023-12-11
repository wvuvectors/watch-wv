#! /usr/bin/env perl


use strict;
use warnings;
use Text::CSV qw(csv);
#use DateTime qw( );
use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use Data::Dumper;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;


my $usage = "\n";
$usage   .= "Usage: $progname [options] UPDIR\n";
$usage   .=   "Compile WaTCH lab run files into database update tables and write them to UPDIR.\n";
$usage   .=   "The file update.batch_files.txt in UPDIR is queried for the list of batch files to include in this update.\n";
$usage   .=   "If update.batch_files.txt does not exist or is empty, $progname will exit with an error status.\n";
$usage   .=   "\n";

my $updir;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } else {
		$updir = $arg;
	}
}
die "FATAL: $progname requires a valid update directory to write its files.\n$usage\n" unless defined $updir and -d "$updir";

# Store paths to batches in this update in an array @batch_files.
my @batch_files = ();
open(my $upFH, "<", "$updir/update.batch_files.txt") or die "Unable to open $updir/update.batch_files.txt to read batch files for processing: $!";
while (my $line = <$upFH>) {
	chomp $line;
	push @batch_files, "$line";
}
close $upFH;

#print Dumper(\@batch_files);
#die;


# Store the column names that we want to print to the set of update files.
my %updatecols    = (
	"assay"         => ["assay_id", "extraction_id", "assay_batch_id", "assay_location_in_batch", "assay_input_ul", "assay_target", "assay_target_genetic_locus", "assay_target_macromolecule", "assay_target_fluorophore", "assay_accepted_droplets", "assay_target_predicted_copies_per_ul_reaction", "assay_target_copies_per_ul_reaction", "assay_comment"], 
	"concentration" => ["concentration_id", "sample_id", "concentration_batch_id", "concentration_location_in_batch", "concentration_comment"], 
	"control"       => ["control_id", "assay_batch_id", "control_type", "control_location_in_batch", "control_input_ul", "control_template", "control_macromolecule", "control_fluorophore", "control_accepted_droplets", "control_predicted_copies_per_ul_reaction", "control_copies_per_ul_reaction", "control_comment"], 
	"control_k"     => ["control_id", "assay_batch_id", "control_type", "assay_location_in_batch", "assay_input_ul", "control_template", "control_macromolecule", "assay_target_fluorophore", "assay_accepted_droplets", "control_predicted_copies_per_ul_reaction", "assay_target_copies_per_ul_reaction", "assay_comment"], 
	"extraction"    => ["extraction_id", "concentration_id", "extraction_batch_id", "extraction_location_in_batch", "extraction_location_in_storage", "extraction_comment"],
	"archive"       => ["archive_id", "sample_id", "archive_batch_id", "archive_location_in_batch", "archive_location_in_storage", "archive_comment"]
);



##############################################
# Now read in the list of batch file paths.
# Key off the file extension to determine if it is data (csv) or metadata (excel) file.

# Assay data are stored in the %well2data hash.
my %well2data = ();

# Metadata are stored in the batches, batchcols, control_wells, and sample2id hashes.
my %batchmeta     = ("assay" => {}, "concentration" => {}, "extraction" => {}, "archive" => {});
my %batchcols     = ("assay" => {}, "concentration" => {}, "extraction" => {}, "archive" => {});
my %control_wells = ();
my %sample2id     = ();


foreach my $filepath (@batch_files) {
	chomp $filepath;
	if ("$filepath" =~ /\.csv$/) {
		# csv files contain the ddPCR data.
		# This populates the well2data hash.
		procAssayData("$filepath", \%well2data);
	} elsif ("$filepath" =~ /\.xlsx$/) {
		# xlsx files contain the run metadata
		my $batch_wkbk = ReadData("$filepath", dtfmt => "mm/dd/yy");
	
		# Parse data based on batch type, which is recorded in cell B1 (page 1, row 1, column 2) of each batch Excel file.
		my @sheet  = Spreadsheet::Read::rows($batch_wkbk->[1]);
		my @typerow = Spreadsheet::Read::row($batch_wkbk->[1], 1);
		my $btype   = lc $typerow[1];
	
		next unless defined $batchmeta{"$btype"};
	
		# Process the batch metadata, including the batch id. This is stored on sheet 1 of each batch Excel file.
		my $batch_id = procBatchMetadata("$btype", \@sheet, $batchmeta{"$btype"}, $batchcols{"$btype"}, \%control_wells);
	
		# Process the batch sample IDs. These are stored on sheet 2 of each batch Excel file.
		procBatchSamples("$btype", $batch_id, $batch_wkbk, $batchmeta{"$btype"}->{"$batch_id"}, \%sample2id);
	}
	else {
		warn "$filepath is an unknown file type and will be ignored!\n";
	}
}
#print Dumper(\%batchmeta);
#print Dumper(\%well2data);
#print Dumper(\%control_wells);
#die;



# Write batch data to tab-delimited files.
my %batch2bfile = ("assay" => "abatch", "concentration" => "cbatch", "extraction" => "ebatch", "archive" => "rbatch");
foreach my $btype (keys %batchmeta) {
	my $bfilename = $batch2bfile{"$btype"};
	open (my $pfh, ">", "$updir/update.$bfilename.txt") or die "Unable to open $updir/update.$bfilename.txt for writing: $!";
	my @cols = sort keys %{$batchcols{"$btype"}};

	# print the column headers
	for (my $i=0; $i<scalar(@cols); $i++) {
		my $col = $cols[$i];
		next if "assay_input_ul" eq "$col";
		print $pfh "\t" unless $i == 0;
		print $pfh "$col";
	}
	print $pfh "\n";

	#print each row of the batch metadata
	foreach my $bid (keys %{$batchmeta{"$btype"}}) {
		for (my $i=0; $i<scalar(@cols); $i++) {
			my $col = $cols[$i];
			next if "assay_input_ul" eq "$col";
			print $pfh "\t" unless $i == 0;
			print $pfh "$batchmeta{$btype}->{$bid}->{$col}";
		}
		print $pfh "\n";
	}
	close $pfh;
}


# Write the easy, one-to-one update files (concentration, extract, and archive sets).
foreach my $btype (keys %batchmeta) {
	next if "$btype" eq "assay";
	#my $lead = substr $btype, 0, 1;
	open (my $ufh, ">", "$updir/update.$btype.txt") or die "Unable to open $updir/update.$btype.txt for writing: $!";
	print $ufh join("\t", @{$updatecols{"$btype"}}) . "\n";
	foreach my $bid (keys %{$batchmeta{"$btype"}}) {
		foreach my $well (keys %{$batchmeta{"$btype"}->{"$bid"}->{"wells"}}) {
			my $sample_id = $batchmeta{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"sample_id"};
			my $uid = $sample2id{"$sample_id"}->{"${btype}_id"};
			print $ufh "$uid";
			for (my $i=1; $i<scalar(@{$updatecols{"$btype"}}); $i++) {
				my $key = $updatecols{"$btype"}->[$i];
				if (defined $batchmeta{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"}) {
					print $ufh "\t" . $batchmeta{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"};
				} elsif (defined $batchmeta{"$btype"}->{"$bid"}->{"$key"}) {
					print $ufh "\t" . $batchmeta{"$btype"}->{"$bid"}->{"$key"};
				} elsif (defined $sample2id{"$sample_id"}->{"$key"}) {
					print $ufh "\t" . $sample2id{"$sample_id"}->{"$key"};
				} else {
					print $ufh "\t";
				}
			}
			print $ufh "\n";
		}
	}
	close $ufh;
}


# Build simplified fluorophore codes.
# These are appended to the control IDs to make them unique.
my %fabbrevs = ("FAM" 				 => "FA",
								"HEX" 				 => "HE",
								"Cy5" 				 => "CY",
								"Cy5.5" 			 => "CZ",
								"ROX" 				 => "RO", 
								"ATTO 590" 		 => "AT",
								"FAM/HEX" 		 => "FH",
								"Cy5/Cy5.5" 	 => "CC",
								"ROX/ATTO 590" => "RA");


# Write the assay update file.
# This is a more complicated, many-to-one set.
open (my $ufh, ">", "$updir/update.assay.txt") or die "Unable to open $updir/update.assay.txt for writing: $!";
print $ufh join("\t", @{$updatecols{"assay"}}) . "\n";

# Get the first letter of the batch type to use in the unique IDs.
my $lead = substr "assay", 0, 1;

foreach my $bid (keys %{$batchmeta{"assay"}}) {	# Each batch of this type in the metadata hash
	foreach my $well (keys %{$batchmeta{"assay"}->{"$bid"}->{"wells"}}) {	# Each well in the metadata hash
		# Skip wells that we've identified as controls.
		# These will be printed in the next block.
		next if defined $control_wells{$well};
		
		# Get the sample id (aka asset id).
		my $sample_id = $batchmeta{"assay"}->{"$bid"}->{"wells"}->{"$well"}->{"sample_id"};
		# Init an incrementor, for the uid.
		my $count = 1;
		
		# Get the data for this well from the well2data hash.
		# Each well contains a hash keyed on the dye.
		my %well_data = %{$well2data{"$well"}};

		# Loop over the dyes for this well.
		# Each dye represents a single assay and prints to a row.
		foreach my $dye (keys %well_data) {
			# Construct a uid based on sample id and the incrementor.
			my $uid = "$sample_id.${lead}$count";
			$count++;
			# Print the uid.
			print $ufh "$uid";
			
			# Loop over the update column names for this batch type, stored in the updatecols hash.
			# These are in array form to control the print order to the output.
			foreach my $colname (@{$updatecols{"assay"}}) {
				# We pull update values from several hashes.
				# This block tests for the existence in these hashes in order to find the value.
				my $val = $batchmeta{"assay"}->{"$bid"}->{"wells"}->{"$well"}->{"$colname"};
				$val = $batchmeta{"assay"}->{"$bid"}->{"$colname"} unless defined $val;
				$val = $sample2id{"$sample_id"}->{"$colname"} unless defined $val;
				$val = $batchmeta{"assay"}->{"$bid"}->{"fluorophores"}->{"$dye"}->{"$colname"} unless defined $val;
				$val = $well2data{"$well"}->{"$dye"}->{"$colname"} unless defined $val;
				
				print $ufh "\t";
				print $ufh "$val" if defined $val;
			}
			
			print $ufh "\n";
		}
	}
}
close $ufh;



# Write the controls update file.
# This is a more complicated, many-to-one set.
# A control is linked to an assay batch (ie many assays), so we record the control-to-abatch relationship here.
open (my $xfh, ">", "$updir/update.control.txt") or die "Unable to open $updir/update.control.txt for writing: $!";
print $xfh join("\t", @{$updatecols{"control"}}) . "\n";

# Get the first letter of the batch type to use in the unique IDs.
$lead = substr "control", 0, 1;

# Each assay batch in the metadata hash has its own set of controls, keyed by dye.
# Control data are stored as an array of hashes.
foreach my $bid (keys %{$batchmeta{"assay"}}) {	# Each assay batch in the metadata hash
	foreach my $dye (keys %{$batchmeta{"assay"}->{"$bid"}->{"fluorophores"}}) {	# Each dye in this assay batch

		# Init an incrementor, for the uid.
		my $count = 1;
		my $dshrt = $fabbrevs{"$dye"};

		foreach my $controlRef (@{$batchmeta{"assay"}->{"$bid"}->{"fluorophores"}->{"$dye"}}) {	# Each control for this dye in this assay batch
			
			my @control_wells = split /,\s*/, "$controlRef->{well}", -1;
			foreach my $well (@control_wells) {
				my $uid = "$bid.${lead}.${dshrt}.$count";
				$count++;
				# Print the uid.
				print $xfh "$uid";
		
				# Loop over the update column names for this batch type, stored in the updatecols hash.
				# These are in array form to control the print order to the output.
				foreach my $colname (@{$updatecols{"assay"}}) {
					# We pull update values from several hashes.
					# This block tests for the existence in these hashes in order to find the value.
					my $val = $controlRef{"control_$colname"};
					$val = $batchmeta{"assay"}->{"$bid"}->{"wells"}->{"$well"}->{"$colname"} unless defined $val;
					$val = $batchmeta{"assay"}->{"$bid"}->{"$colname"} unless defined $val;
					$val = $sample2id{"$sample_id"}->{"$colname"} unless defined $val;
					$val = $batchmeta{"assay"}->{"$bid"}->{"fluorophores"}->{"$dye"}->{"$colname"} unless defined $val;
					$val = $well2data{"$well"}->{"$dye"}->{"$colname"} unless defined $val;
				
					print $xfh "\t";
					print $xfh "$val" if defined $val;
				}	# update columns
				print $ufh "\n";
			}	# control wells
		}	# controls array
	}	# dyes
}	# assay batches

close $xfh;

exit 0;


sub procAssayData {
	my $fp           = shift;
	my $well2dataRef = shift;

	my $csvA  = Text::CSV->new({auto_diag => 4, binary => 1});
	
	open (my $afh, "<", "$fp") or die "procData sub is unable to open $fp for reading: $!\n";
	my %fieldHash = ();
	my $lineCount = 0;

	while (my $line = $csvA->getline($afh)) {
		my @cols = @$line;
		if ($lineCount == 0) {
			if ("$cols[0]" =~ /Well/) {
				$cols[0] = "Well";
				# Record the cell headers so we can access the data by hash.
				my $n = 0;
				foreach (@cols) {
					$fieldHash{"$_"} = $n;
					$n++;
				}
			} else {
				# Not a data file so close the fh and dump out.
				close $afh;
				last;
			}
		} else {
#			my ($cell, $result, $dye, $accepted_droplets) = ($cols[0], $cols[6], $cols[12], $cols[13]);
			my $well   = $cols[$fieldHash{"Well"}] if defined $fieldHash{"Well"};
			my $t_conc = $cols[$fieldHash{"Conc(copies/uL)"}] if defined $fieldHash{"Conc(copies/uL)"};
			my $dye    = $cols[$fieldHash{"DyeName(s)"}] if defined $fieldHash{"DyeName(s)"};
			my $axdrop = $cols[$fieldHash{"Accepted Droplets"}] if defined $fieldHash{"Accepted Droplets"};
			if (defined $well and defined $t_conc and defined $dye and defined $axdrop) {
				$well2dataRef->{"$well"} = {} unless defined $well2dataRef->{"$well"};
				$well2dataRef->{"$well"}->{"$dye"} = {
					"assay_target_copies_per_ul_reaction" => $t_conc, 
					"assay_accepted_droplets"             => $axdrop
				};
			} else {
				print "***** ERROR : Malformed line in data file. *****\n";
				print "File: $fp.\n";
				print "Line: $lineCount.\n";
				print "Data: " . join(", ", @$line) . "\n";
				print "***** ACTION: This line will be ignored. *****\n";
			}
		}
		$lineCount++;
	}
	close $afh;
	
}


sub procBatchMetadata {
	my $btype           = shift;
	my $sheetRef        = shift;
	my $batchmetaRef    = shift;
	my $batchcolsRef    = shift;
	my $controlwellsRef = shift;
	
	$btype = lc $btype;
	
	my $bid;
	my @dyes = ();
	
	my %local_batch  = ();
	my %local_fluoro = ();
	
	for (my $i=1; $i < scalar(@$sheetRef); $i++) {
		next if "" eq "$sheetRef->[$i][0]";
		my $key = lc "$sheetRef->[$i][0]";
		my $val = "$sheetRef->[$i][1]";

		push @dyes, "$val" if "$key" eq "target fluorophore";

		if ("$key" eq "batch id") {
			$bid = "$val";
		} elsif ("$key" eq "date") {
			my $datetime = DateTime::Format::Excel->parse_datetime($val);
			$val = $datetime->ymd('-');
		}
		$key =~ s/ /_/gi;
		$key = "${btype}_$key";
		if (scalar @dyes > 0) {
			$local_fluoro{$dyes[-1]} = {} unless defined $local_fluoro{$dyes[-1]};
			$local_fluoro{$dyes[-1]}->{"$key"} = "$val";
		} else {
			$local_batch{"$key"} = "$val";
		}
	}
	
	#if ("$btype" eq "assay") {
	#	print Dumper(\%local_fluoro);
	#	die;
	#}
	
	die "FATAL: No batch id recognized: $btype\n" unless defined $bid;
	
	$batchmetaRef->{"$bid"} = {} unless defined $batchmetaRef->{"$bid"};
	
	# placeholder to hold sample data (keyed per well)	
	$batchmetaRef->{"$bid"}->{"wells"} = {};
	
	foreach my $key (keys %local_batch) {
		$batchmetaRef->{"$bid"}->{"$key"} = "$local_batch{$key}";
		$batchcolsRef->{"$key"} = 1;
	}
	
	if (scalar keys %local_fluoro > 0) {
		$batchmetaRef->{"$bid"}->{"fluorophores"} = {};
		foreach my $dye (keys %local_fluoro) {
			$batchmetaRef->{"$bid"}->{"fluorophores"}->{$dye} = {"controls" => []};
			my @controls = ();
			foreach my $key (keys %{$local_fluoro{$dye}}) {
				if ("$key" =~ /^assay_control_(\d+)_(.+)$/) {
					my $rep   = $1;
					my $field = lc $2;
					$controls[$rep] = {} unless defined $controls[$rep];
					$controls[$rep]->{"$field"} = "$local_fluoro{$dye}->{$key}";
					# Store the control wells as keys in the control_wells hash.
					# This makes it easier to split assay wells from control wells on output.
					if ("$field" eq "well") {
						my @ctl_wells = split /,\s*/, "$local_fluoro{$dye}->{$key}", -1;
						for (my $i=0; $i<scalar(@ctl_wells); $i++) {
							my $ctl_well = $ctl_wells[$i];
							if ("$ctl_well" =~ /([A-Za-z])(\d+)/) {
								my ($row, $col) = ($1, $2);
								$row = uc $row;
								$col = "0$col" if length $col < 2;
								$ctl_well = "$row$col";
								$controlwellsRef->{"$ctl_well"} = 1;
							}
						}
					}
				} else {
					$batchmetaRef->{"$bid"}->{"fluorophores"}->{$dye}->{"$key"} = "$local_fluoro{$dye}->{$key}";
				}
			}
			
			# Prune controls that are undefined or have type NA before storing in the batch meta hash.
			# First we store the indices of controls that are not defined (the first) or have type NA.
			my @nas = ();
			for (my $i=0; $i<scalar(@controls); $i++) {
				push @nas, $i if !defined $controls[$i] or $controls[$i]->{"type"} eq "NA";
			}
			# Then we splice out those indices from the controls array.
			# This must be done in reverse order to preserve the index order!
			foreach my $index (reverse @nas) {
				splice @controls, $index, 1;
			}
			# Store the pruned controls array in the batch meta hash.
			$batchmetaRef->{"$bid"}->{"fluorophores"}->{$dye}->{"controls"} = \@controls unless scalar(@controls) == 0;
		}
	}
	
	return $bid;
}	



sub procBatchSamples {
	my $btype         = shift;
	my $bid           = shift;
	my $batch_wkbkRef = shift;
	my $batchmetaRef  = shift;
	my $sample2idRef  = shift;
	
	$btype = lc $btype;
	
	# The location (order) of sheets in each batch excel file.
	my %sheet2batch = ("assay"         => {"metadata" => 1, "sample ids" => 2, "backrefs" => 3, "comments" => 4, "volume overrides" => 5}, # backrefs are extraction ids, defaults to 1
										 "concentration" => {"metadata" => 1, "sample ids" => 2, "spikes" => 3, "comments" => 4}, 
										 "extraction"    => {"metadata" => 1, "sample ids" => 2, "backrefs" => 3, "comments" => 4, "storage ids" => 5}, # backrefs are concentration ids, defaults to 1
										 "archive"       => {"metadata" => 1, "sample ids" => 2, "comments" => 3, "storage ids" => 4});


	# Map the source of each set.
	my %ref2source = ("assay" => "extraction", "extraction" => "concentration", "concentration" => "sample");

	my @row2name = ("", "A", "B", "C", "D", "E", "F", "G", "H");
	
	my %position2sheet = %{$sheet2batch{"$btype"}};
	
	my @sample_rows  = Spreadsheet::Read::rows($batch_wkbkRef->[$position2sheet{"sample ids"}]);
	my @comment_rows = Spreadsheet::Read::rows($batch_wkbkRef->[$position2sheet{"comments"}]);
	my @backref_rows = Spreadsheet::Read::rows($batch_wkbkRef->[$position2sheet{"backrefs"}]) if defined $position2sheet{"backrefs"};
	my @volume_rows  = Spreadsheet::Read::rows($batch_wkbkRef->[$position2sheet{"volume overrides"}]) if defined $position2sheet{"volume overrides"};
	my @storage_rows = Spreadsheet::Read::rows($batch_wkbkRef->[$position2sheet{"storage ids"}]) if defined $position2sheet{"storage ids"};
	
	# all batch sheets except the Metadata sheet are based off sample_rows (ie the "Sample IDs" sheet), so we only need to loop over that
	for (my $i=1; $i < scalar(@sample_rows); $i++) {
		for (my $j=1; $j < scalar(@{$sample_rows[$i]}); $j++) {
			if (!defined $sample_rows[$i][$j]) {
				warn "****WARNING**** $btype batch has no sample ID at row $i, column $j.";
				next;
			} else {
				my $sample_id = "$sample_rows[$i][$j]";
				$sample_id =~ s/^\s+//;
				$sample_id =~ s/\s+$//;
				
				my $batch_loc = "$row2name[$i]" . sprintf("%02d", $j);
				
				unless ("" eq "$sample_id" or "buffer" eq lc("$sample_id")) {
					$batchmetaRef->{"wells"}->{"$batch_loc"} = {"${btype}_location_in_batch" => "$batch_loc", "${btype}_batch_id" => "$bid"};
					# record the sample id
					$batchmetaRef->{"wells"}->{"$batch_loc"}->{"sample_id"} = "$sample_id";

					# make an id for this entry in the batch, based on the sample id
					my $lead = substr $btype, 0, 1;
					$sample2idRef->{"$sample_id"} = {} unless defined $sample2idRef->{"$sample_id"};
					my $num = 1;
					my $uid = "$sample_id.${lead}${num}";
					$sample2idRef->{"$sample_id"}->{"${btype}_id"} = "$uid";

					if (defined $comment_rows[0]) {
						# record any comment on the sample in this batch
						if (defined $comment_rows[$i][$j] and "$comment_rows[$i][$j]" ne "") {
							$batchmetaRef->{"wells"}->{"$batch_loc"}->{"${btype}_comments"} = "$comment_rows[$i][$j]";
						}
					}
					
					if (defined $volume_rows[0]) {
						# Record the input volume
						if (defined $volume_rows[$i][$j] and "$volume_rows[$i][$j]" ne "") {
							$batchmetaRef->{"wells"}->{"$batch_loc"}->{"${btype}_input_ul"} = "$volume_rows[$i][$j]";
						} else {
							$batchmetaRef->{"wells"}->{"$batch_loc"}->{"${btype}_input_ul"} = $batchmetaRef->{"${btype}_input_ul"};
						}
					}
					
					if (defined $storage_rows[0]) {
						# Record the storage location
						if (defined $storage_rows[$i][$j] and "$storage_rows[$i][$j]" ne "") {
							$batchmetaRef->{"wells"}->{"$batch_loc"}->{"${btype}_location_in_storage"} = "$storage_rows[$i][$j]";
						} else {
							$batchmetaRef->{"wells"}->{"$batch_loc"}->{"${btype}_input_ul"} = $batchmetaRef->{"${btype}_input_ul"};
						}
					}

					if (defined $backref_rows[0]) {
						# Record the back reference.
						# For extractions the source is the concentration. For assays the source is the extraction.
						my $source      = $ref2source{"$btype"};
						my $source_code = substr "$source", 0, 1;
						my $num         = 1;
						$num = "$backref_rows[$i][$j]" if defined $backref_rows[$i][$j] and "$backref_rows[$i][$j]" ne "";
						$batchmetaRef->{"wells"}->{"$batch_loc"}->{"${source}_id"} = "${sample_id}.${source_code}${num}";
					}

				}
			}
		}	# j loop
	}	# i loop
}

