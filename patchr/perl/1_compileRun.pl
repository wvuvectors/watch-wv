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
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=   "Compile WaTCH lab run files into database update tables.\n";
$usage   .=   "RUNDIR must contain exactly one ABATCH Excel file and one ddPCR results file, in csv format.\n";
$usage   .=   "RUNDIR can contain zero, one, or multiple CBATCH, EBATCH, and RBATCH files as needed.\n";
$usage   .=   "A new directory RUNDIR/updates will be created to hold the database update files.\n";
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

opendir(my $dirHa, $rundir) or die "Unable to open directory $rundir: $!";
my @rund_files  = grep { (/\.csv$/) && (!/^~/) && -f "$rundir/$_" } readdir($dirHa);
closedir $dirHa;

opendir(my $dirHb, $rundir) or die "Unable to open directory $rundir: $!";
my @batch_files = grep { (/\.xlsx$/) && (!/^~/) && -f "$rundir/$_" } readdir($dirHb);
closedir $dirHb;


my %batches       = ("assay" => {}, "concentration" => {}, "extraction" => {}, "archive" => {});
my %batchcols     = ("assay" => {}, "concentration" => {}, "extraction" => {}, "archive" => {});
my %control_wells = ();
my %sample2id     = ();

my %updatecols    = (
	"assay"         => ["assay_id", "extraction_id", "assay_batch_id", "assay_location_in_batch", "assay_input_ul", "assay_target", "assay_target_genetic_locus", "assay_target_macromolecule", "assay_target_fluorophore", "assay_accepted_droplets", "assay_target_predicted_copies_per_ul_reaction", "assay_target_copies_per_ul_reaction", "assay_comment"], 
	"concentration" => ["concentration_id", "sample_id", "concentration_batch_id", "concentration_location_in_batch", "concentration_comment"], 
	"control"       => ["control_id", "assay_batch_id", "control_type", "control_location_in_batch", "control_input_ul", "control_template", "control_macromolecule", "control_fluorophore", "control_accepted_droplets", "control_predicted_copies_per_ul_reaction", "control_copies_per_ul_reaction", "control_comment"], 
	"control_k"     => ["control_id", "assay_batch_id", "control_type", "assay_location_in_batch", "assay_input_ul", "control_template", "control_macromolecule", "assay_target_fluorophore", "assay_accepted_droplets", "control_predicted_copies_per_ul_reaction", "assay_target_copies_per_ul_reaction", "assay_comment"], 
	"extraction"    => ["extraction_id", "concentration_id", "extraction_batch_id", "extraction_location_in_batch", "extraction_location_in_storage", "extraction_comment"],
	"archive"       => ["archive_id", "sample_id", "archive_batch_id", "archive_location_in_batch", "archive_location_in_storage", "archive_comment"]
);

my %fluorcols = ();

# id link map
my %ref2source = ("assay" => "extraction", "extraction" => "concentration", "concentration" => "sample");

# Batch sheet orders
my %sheet2batch = ("assay"         => {"metadata" => 1, "sample ids" => 2, "backrefs" => 3, "comments" => 4, "volume overrides" => 5}, # backrefs are extraction ids, defaults to 1
									 "concentration" => {"metadata" => 1, "sample ids" => 2, "comments" => 3}, 
									 "extraction"    => {"metadata" => 1, "sample ids" => 2, "backrefs" => 3, "comments" => 4, "storage ids" => 5}, # backrefs are concentration ids, defaults to 1
									 "archive"       => {"metadata" => 1, "sample ids" => 2, "comments" => 3, "storage ids" => 4});


# Read in the run data file(s) (csv)
# There many be multiple csv files in the run directory so we want to test each one
# ATM there can be only one data file per run; that may change in the future?
#
#Well,Sample description 1,Sample description 2,Sample description 3,Sample description 4,Target,Conc(copies/uL),Status,Experiment,SampleType,TargetType,Supermix,DyeName(s),Accepted Droplets,Positives,Negatives
#my @plate2wells = ();
my %well2data = ();
foreach my $rundF (@rund_files) {
	my $csvA  = Text::CSV->new({auto_diag => 4, binary => 1});
	open (my $afh, "<", "$rundir/$rundF") or die "Unable to open $rundir/$rundF for reading: $!\n";
	#my %well2data = ();
	my %fieldHash = ();
	my $lineCount = 0;
	while (my $line = $csvA->getline($afh)) {
		my @cols = @$line;
		if ($lineCount == 0) {
			if ("$cols[0]" eq "Well") {
				# record the cell headers so we can access the data by hash
				my $n = 0;
				foreach (@cols) {
					$fieldHash{"$_"} = $n;
					$n++;
				}
			} else {
				# Not a data file so close the fh and skip to the next
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
				$well2data{"$well"} = {} unless defined $well2data{"$well"};
				$well2data{"$well"}->{"$dye"} = {
					"assay_target_copies_per_ul_reaction" => $t_conc, 
					"assay_accepted_droplets"             => $axdrop
				};
			} else {
				print "***** ERROR : Malformed line in run data. *****\n";
				print "File: $rundF.\n";
				print "Line: $lineCount.\n";
				print "Data: " . join(", ", @$line) . "\n";
				print "***** ACTION: This line will be ignored. *****\n";
			}
		}
		$lineCount++;
	}
	close $afh;
	#push @plate2wells, \%well2data;
}



foreach (@batch_files) {
	my $batch_wkbk = ReadData("$rundir/$_", dtfmt => "mm/dd/yy");
	
	# Parse data based on batch type, which is recorded in cell B1 (page 1, row 1, column 2) of the Excel file
	my @mdrows  = Spreadsheet::Read::rows($batch_wkbk->[1]);
	my @typerow = Spreadsheet::Read::row($batch_wkbk->[1], 1);
	my $btype   = lc $typerow[1];
	
	next unless defined $batches{"$btype"};
	
	# Read the batch metadata (on sheet 1 of the Excel file)
	my $batch_id = read_batch_metadata("$btype", \@mdrows, $batches{"$btype"}, $batchcols{"$btype"});
	
	# Now read the samples in each plate
	read_batch_data("$btype", $batch_id, $batch_wkbk, $batches{"$btype"}->{"$batch_id"}, \%sample2id);
}

#print Dumper(\%batches);
#print Dumper(\%control_wells);
#print Dumper(\%sample2id);
#print Dumper(\%batchcols);
#die;


# Write batch metadata to update files
`rm -r "$rundir/updates/"` if -d "$rundir/updates/";
`mkdir "$rundir/updates/"`;

my %batchf = ("assay" => "abatch", "concentration" => "cbatch", "extraction" => "ebatch", "archive" => "rbatch");
foreach my $btype (keys %batches) {
	my $bname = $batchf{"$btype"};
	open (my $pfh, ">", "$rundir/updates/update.$bname.txt") or die "Unable to open $rundir/updates/update.$bname.txt for writing: $!";
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
	foreach my $bid (keys %{$batches{"$btype"}}) {
		for (my $i=0; $i<scalar(@cols); $i++) {
			my $col = $cols[$i];
			next if "assay_input_ul" eq "$col";
			print $pfh "\t" unless $i == 0;
			print $pfh "$batches{$btype}->{$bid}->{$col}";
		}
		print $pfh "\n";
	}
	close $pfh;
}


# write one-to-one update files (conc, extract, archive)
foreach my $btype (keys %batches) {
	next if "$btype" eq "assay";
	#my $lead = substr $btype, 0, 1;
	open (my $ufh, ">", "$rundir/updates/update.$btype.txt") or die "Unable to open $rundir/updates/update.$btype.txt for writing: $!";
	print $ufh join("\t", @{$updatecols{"$btype"}}) . "\n";
	foreach my $bid (keys %{$batches{"$btype"}}) {
		foreach my $well (keys %{$batches{"$btype"}->{"$bid"}->{"wells"}}) {
			my $sample_id = $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"sample_id"};
			my $uid = $sample2id{"$sample_id"}->{"${btype}_id"};
			print $ufh "$uid";
			for (my $i=1; $i<scalar(@{$updatecols{"$btype"}}); $i++) {
				my $key = $updatecols{"$btype"}->[$i];
				if (defined $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"}) {
					print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"};
				} elsif (defined $batches{"$btype"}->{"$bid"}->{"$key"}) {
					print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"$key"};
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

# Simplified fluorophore codes, mostly to append to the control IDs to make them unique
my %fabbrevs = ("FAM" => "FA",
								"HEX" => "HE",
								"Cy5" => "CY",
								"Cy5.5" => "CZ",
								"ROX" => "RO", 
								"ATTO 590" => "AT",
								"FAM/HEX" => "FH",
								"Cy5/Cy5.5" => "CC",
								"ROX/ATTO 590" => "RA");


# write assay and control update files, a many-to-one set
my $btype = "assay";
open (my $ufh, ">", "$rundir/updates/update.$btype.txt") or die "Unable to open $rundir/updates/update.$btype.txt for writing: $!";
print $ufh join("\t", @{$updatecols{"$btype"}}) . "\n";
my $lead = substr $btype, 0, 1;
open (my $xfh, ">", "$rundir/updates/update.control.txt") or die "Unable to open $rundir/updates/update.control.txt for writing: $!";
print $xfh join("\t", @{$updatecols{"control"}}) . "\n";

foreach my $bid (keys %{$batches{"$btype"}}) {
	foreach my $well (keys %{$batches{"$btype"}->{"$bid"}->{"wells"}}) {

		if (defined $control_wells{$well}) {
			#die "$well\n";
			foreach my $fluor (keys %{$batches{"$btype"}->{"$bid"}->{"fluorophores"}}) {
				next unless defined $control_wells{$well}->{"$fluor"};
				my $rep = $control_wells{$well}->{"$fluor"}->{"rep"};
				my $sid = $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"sample_id"};
				$sid =~ s/\s//gi;
				my $fabbrev = $fabbrevs{"$fluor"};
				my $uid = "$bid.$sid.$fabbrev$rep";
				print $xfh "$uid";
				for (my $i=1; $i<scalar(@{$updatecols{"control_k"}}); $i++) {
					my $key = $updatecols{"control_k"}->[$i];
					my $fkey = "${key}_${rep}";
					if (defined $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"}) {
						print $xfh "\t" . $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"};
					} elsif (defined $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$key"}) {
						print $xfh "\t" . $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$key"};
					} elsif (defined $well2data{"$well"}->{"$fluor"}->{"$key"}) {
						print $xfh "\t" . $well2data{"$well"}->{"$fluor"}->{"$key"};
					} elsif (defined $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$fkey"}) {
						print $xfh "\t" . $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$fkey"};
					} elsif (defined $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$fkey"}) {
						print $xfh "\t" . $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$fkey"};
					} elsif (defined $well2data{"$well"}->{"$fluor"}->{"$fkey"}) {
						print $xfh "\t" . $well2data{"$well"}->{"$fluor"}->{"$fkey"};
					} else {
						print $xfh "\t";
					}
				}
				print $xfh "\n";
			}

		} else {
			# handle assay data
			my $sample_id = $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"sample_id"};
			my $count = 1;
			foreach my $fluor (keys %{$batches{"$btype"}->{"$bid"}->{"fluorophores"}}) {
				my $uid = "$sample_id.${lead}$count";
				$count++;
				print $ufh "$uid";
				for (my $i=1; $i<scalar(@{$updatecols{"$btype"}}); $i++) {
					my $key = $updatecols{"$btype"}->[$i];
					if (defined $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"}) {
						print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"wells"}->{"$well"}->{"$key"};
					} elsif (defined $batches{"$btype"}->{"$bid"}->{"$key"}) {
						print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"$key"};
					} elsif (defined $sample2id{"$sample_id"}->{"$key"}) {
						print $ufh "\t" . $sample2id{"$sample_id"}->{"$key"};
					} elsif (defined $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$key"}) {
						print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$key"};
					} elsif (defined $well2data{"$well"}->{"$fluor"}->{"$key"}) {
						print $ufh "\t" . $well2data{"$well"}->{"$fluor"}->{"$key"};
					} else {
						print $ufh "\t";
					}
				}
				print $ufh "\n";
			}
		}
		
	}
}
close $ufh;
close $xfh;



sub read_batch_metadata {
	my $btype = shift;
	my $sheet = shift;
	my $batch = shift;
	my $bcols = shift;
	
	$btype = lc $btype;
	
	my $bid;
	my @dyes = ();
	
	my %local_batch  = ();
	my %local_fluoro = ();
	
	for (my $i=1; $i < scalar(@$sheet); $i++) {
		next if "" eq "$sheet->[$i][0]";
		my $key = lc "$sheet->[$i][0]";
		my $val = "$sheet->[$i][1]";

		push @dyes, "$val" if "target fluorophore" eq "$key";

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
	
#	if (scalar keys %local_fluoro > 0) {
#		print Dumper(\%local_fluoro);
#		die;
#	}

	die "FATAL: No batch id recognized: $btype\n" unless defined $bid;
	
	$batch->{"$bid"} = {} unless defined $batch->{"$bid"};
	
	# placeholder to hold sample data (keyed per well)	
	$batch->{"$bid"}->{"wells"} = {};
	
	foreach my $key (keys %local_batch) {
		$batch->{"$bid"}->{"$key"} = "$local_batch{$key}";
		$bcols->{"$key"} = 1;
	}
	if (scalar keys %local_fluoro > 0) {
		$batch->{"$bid"}->{"fluorophores"} = {};
		foreach my $dye (keys %local_fluoro) {
			$batch->{"$bid"}->{"fluorophores"}->{$dye} = {};
			foreach my $key (keys %{$local_fluoro{$dye}}) {
				if ("$key" =~ /^assay_control_(\d+)_(.+)$/) {
					my $rep   = $1;
					my $field = $2;
					my $key_new = "control_${field}_${rep}";
					if ("well" eq "$field") {
						my $ctl_well = "$local_fluoro{$dye}->{$key}";
						if ("$ctl_well" =~ /([A-Za-z])(\d+)/) {
							my ($row, $col) = ($1, $2);
							$row = uc $row;
							$col = "0$col" if $col < 10;
							$ctl_well = "$row$col";
							$control_wells{"$ctl_well"} = {} unless defined $control_wells{"$ctl_well"};
							$control_wells{"$ctl_well"}->{"$dye"} = {"rep" => $rep};
						}
					}
					$batch->{"$bid"}->{"fluorophores"}->{$dye}->{"$key_new"} = "$local_fluoro{$dye}->{$key}";
				} else {
					$batch->{"$bid"}->{"fluorophores"}->{$dye}->{"$key"} = "$local_fluoro{$dye}->{$key}";
				}
			}
		}
	}
	
	return $bid;
}	



sub read_batch_data {
	my $btype = shift;
	my $bid   = shift;
	my $wkbk  = shift;
	my $batch = shift;
	my $s2id  = shift;
	
	$btype = lc $btype;
	
	my @row2name = ("", "A", "B", "C", "D", "E", "F", "G", "H");
	
	my %position2sheet = %{$sheet2batch{"$btype"}};
	
	my @sample_rows  = Spreadsheet::Read::rows($wkbk->[$position2sheet{"sample ids"}]);
	my @comment_rows = Spreadsheet::Read::rows($wkbk->[$position2sheet{"comments"}]);
	my @backref_rows = Spreadsheet::Read::rows($wkbk->[$position2sheet{"backrefs"}]) if defined $position2sheet{"backrefs"};
	my @volume_rows  = Spreadsheet::Read::rows($wkbk->[$position2sheet{"volume overrides"}]) if defined $position2sheet{"volume overrides"};
	my @storage_rows = Spreadsheet::Read::rows($wkbk->[$position2sheet{"storage ids"}]) if defined $position2sheet{"storage ids"};
	
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
					$batch->{"wells"}->{"$batch_loc"} = {"${btype}_location_in_batch" => "$batch_loc", "${btype}_batch_id" => "$bid"};
					# record the sample id
					$batch->{"wells"}->{"$batch_loc"}->{"sample_id"} = "$sample_id";

					# make an id for this entry in the batch, based on the sample id
					my $lead = substr $btype, 0, 1;
					$s2id->{"$sample_id"} = {} unless defined $s2id->{"$sample_id"};
					my $num = 1;
					my $uid = "$sample_id.${lead}${num}";
					$s2id->{"$sample_id"}->{"${btype}_id"} = "$uid";

					if (defined $comment_rows[0]) {
						# record any comment on the sample in this batch
						if (defined $comment_rows[$i][$j] and "$comment_rows[$i][$j]" ne "") {
							$batch->{"wells"}->{"$batch_loc"}->{"${btype}_comments"} = "$comment_rows[$i][$j]";
						}
					}
					
					if (defined $volume_rows[0]) {
						# Record the input volume
						if (defined $volume_rows[$i][$j] and "$volume_rows[$i][$j]" ne "") {
							$batch->{"wells"}->{"$batch_loc"}->{"${btype}_input_ul"} = "$volume_rows[$i][$j]";
						} else {
							$batch->{"wells"}->{"$batch_loc"}->{"${btype}_input_ul"} = $batch->{"${btype}_input_ul"};
						}
					}
					
					if (defined $storage_rows[0]) {
						# Record the storage location
						if (defined $storage_rows[$i][$j] and "$storage_rows[$i][$j]" ne "") {
							$batch->{"wells"}->{"$batch_loc"}->{"${btype}_location_in_storage"} = "$storage_rows[$i][$j]";
						} else {
							$batch->{"wells"}->{"$batch_loc"}->{"${btype}_input_ul"} = $batch->{"${btype}_input_ul"};
						}
					}

					if (defined $backref_rows[0]) {
						# Record the back reference.
						# For extractions the source is the concentration. For assays the source is the extraction.
						my $source      = $ref2source{"$btype"};
						my $source_code = substr "$source", 0, 1;
						my $num         = 1;
						$num = "$backref_rows[$i][$j]" if defined $backref_rows[$i][$j] and "$backref_rows[$i][$j]" ne "";
						$batch->{"wells"}->{"$batch_loc"}->{"${source}_id"} = "${sample_id}.${source_code}${num}";
					}

				}
			}
		}	# j loop
	}	# i loop
}


exit 0;

