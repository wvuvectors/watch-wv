#! /usr/bin/env perl

use strict;
use warnings;
use Text::CSV;
#use DateTime qw( );
use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use Data::Dumper;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=   "Merge WaTCH testing laboratory Excel files into a single table.\n";
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

my $rundf = "$rundir/run_data.csv";

opendir(my $dirH, $rundir) or die "Unable to open directory $rundir: $!";
my @batch_files = grep { (/\.xlsx$/) && (!/^~/) && -f "$rundir/$_" } readdir($dirH);
closedir $dirH;

my %batches       = ("assay" => {}, "concentration" => {}, "extraction" => {}, "archive" => {});
my %batchcols    = ("assay" => {}, "concentration" => {}, "extraction" => {}, "archive" => {});
my %sample2id    = ();

my %updatecols   = (
	"assay"         => ["assay_id", "extraction_id", "assay_batch_id", "assay_location_in_batch", "assay_input_ul", "assay_target", "assay_target_category", "assay_target_genetic_locus", "assay_target_macromolecule", "assay_target_fluorophore", "assay_accepted_droplets", "assay_target_calculated_copies_per_ul_reaction", "assay_target_copies_per_ul_reaction", "assay_comment"], 
	"concentration" => ["concentration_id", "sample_id", "concentration_batch_id", "concentration_location_in_batch", "concentration_comment"], 
#	"concentration" => ["concentration_id", "sample_id", "concentration_batch_id", "concentration_location_in_batch", "concentration_input_ml", "concentration_output_ml", "concentration_comment"], 
	"control"         => ["control_id", "assay_batch_id", "control_type", "control_location_in_batch", "control_template", "control_target_macromolecule", "control_target_fluorophore", "control_accepted_droplets", "control_target_calculated_copies_per_ul_reaction", "control_target_copies_per_ul_reaction", "control_comment"], 
#	"extraction"    => ["extraction_id", "concentration_id", "extraction_batch_id", "extraction_location_in_batch", "extraction_location_in_storage", "extraction_input_ul", "extraction_output_ul", "extraction_comment"],
	"extraction"    => ["extraction_id", "concentration_id", "extraction_batch_id", "extraction_location_in_batch", "extraction_location_in_storage", "extraction_comment"],
	"archive"       => ["archive_id", "sample_id", "archive_batch_id", "archive_location_in_batch", "archive_location_in_storage", "archive_comment"]
);

my %fluorcols = ();


# read in assay data file, run_data.csv
# Well,Sample description 1,Sample description 2,Sample description 3,Sample description 4,Target,Conc(copies/uL),Status,Experiment,SampleType,TargetType,Supermix,DyeName(s),Accepted Droplets,Positives,Negatives
my %cell2data = ();
my $csvA = Text::CSV->new({auto_diag => 4, binary => 1});
open (my $afh, "<", "$rundf") or die "Unable to open $rundf for reading: $!\n";
my $count = 0;
while (my $line = $csvA->getline($afh)) {
	my @cols = @$line;
	$count++;
	if (scalar @cols < 15) {
		warn "WARN  : Line $count in $rundf has " . scalar(@cols) . " columns but requires *exactly* 15. Skipping.\n";
		next;
	}
	my ($cell, $result, $dye, $accepted_droplets) = ($cols[0], $cols[6], $cols[12], $cols[13]);
	$cell2data{"$cell"} = {} unless defined $cell2data{"$cell"};
	$cell2data{"$cell"}->{"$dye"} = {"assay_target_copies_per_ul_reaction" => $result, "assay_accepted_droplets" => $accepted_droplets};
}
close $afh;


foreach (@batch_files) {
	my $batch_wkbk = ReadData("$rundir/$_", dtfmt => "mm/dd/yy");
	
	# Parse data based on batch type, which is recorded in cell B1 (page 1, row 1, column 2) of the Excel file
	my @mdrows = Spreadsheet::Read::rows($batch_wkbk->[1]);
	my @typerow = Spreadsheet::Read::row($batch_wkbk->[1], 1);
	my $btype = lc $typerow[1];
	
	next unless defined $batches{"$btype"};
	
	# Read the batch metadata (on sheet 1 of the Excel file)
	my $batch_id = read_batch_metadata("$btype", \@mdrows, $batches{"$btype"}, $batchcols{"$btype"});
	
	# Now read the samples in each plate
	read_batch_data("$btype", $batch_id, $batch_wkbk, $batches{"$btype"}->{"$batch_id"}, \%sample2id);
}

#print Dumper(\%cell2data);
#print Dumper(\%batches);
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
		foreach my $cell (keys %{$batches{"$btype"}->{"$bid"}->{"cells"}}) {
			my $sample_id = $batches{"$btype"}->{"$bid"}->{"cells"}->{"$cell"}->{"sample_id"};
			my $uid = $sample2id{"$sample_id"}->{"${btype}_id"};
			print $ufh "$uid";
			for (my $i=1; $i<scalar(@{$updatecols{"$btype"}}); $i++) {
				my $key = $updatecols{"$btype"}->[$i];
				if (defined $batches{"$btype"}->{"$bid"}->{"cells"}->{"$cell"}->{"$key"}) {
					print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"cells"}->{"$cell"}->{"$key"};
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

# write assay update file, a many-to-one set
# assay_id
# extraction_id
# assay_plate_id
# assay_plate_cell
# assay_input_ml
# assay_type
# assay_target
# assay_target_category
# assay_target_genetic_locus
# assay_target_macromolecule
# assay_target_fluorophore
# assay_accepted_droplets
# assay_calculated_copies_per_reaction
# assay_copies_per_reaction
# assay_result
# assay_comments
my $btype = "assay";
open (my $ufh, ">", "$rundir/updates/update.$btype.txt") or die "Unable to open $rundir/updates/update.$btype.txt for writing: $!";
print $ufh join("\t", @{$updatecols{"$btype"}}) . "\n";
my $lead = substr $btype, 0, 1;
foreach my $bid (keys %{$batches{"$btype"}}) {
	foreach my $cell (keys %{$batches{"$btype"}->{"$bid"}->{"cells"}}) {
		my $sample_id = $batches{"$btype"}->{"$bid"}->{"cells"}->{"$cell"}->{"sample_id"};
		
		# REMOVE CONTROLS
		next if "$sample_id" eq "PCR NC" or "$sample_id" eq "PCR PC";
		
		my $count = 1;
		foreach my $fluor (keys %{$batches{"$btype"}->{"$bid"}->{"fluorophores"}}) {
			my $uid = "$sample_id.${lead}$count";
			$count++;
			print $ufh "$uid";
			for (my $i=1; $i<scalar(@{$updatecols{"$btype"}}); $i++) {
				my $key = $updatecols{"$btype"}->[$i];
				if (defined $batches{"$btype"}->{"$bid"}->{"cells"}->{"$cell"}->{"$key"}) {
					print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"cells"}->{"$cell"}->{"$key"};
				} elsif (defined $batches{"$btype"}->{"$bid"}->{"$key"}) {
					print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"$key"};
				} elsif (defined $sample2id{"$sample_id"}->{"$key"}) {
					print $ufh "\t" . $sample2id{"$sample_id"}->{"$key"};
				} elsif (defined $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$key"}) {
					print $ufh "\t" . $batches{"$btype"}->{"$bid"}->{"fluorophores"}->{"$fluor"}->{"$key"};
				} elsif (defined $cell2data{"$cell"}->{"$fluor"}->{"$key"}) {
					print $ufh "\t" . $cell2data{"$cell"}->{"$fluor"}->{"$key"};
				} else {
					print $ufh "\t";
				}
			}
			print $ufh "\n";
		}
	}
}
close $ufh;


=cut
# write batch assay controls to update file
open (my $cfh, ">", "$rundir/update.controls.txt") or die "Unable to open $rundir/update.controls.txt for writing: $!";
my @cols = sort keys %fluorcols;
print $cfh join("\t", @cols) . "\n";
foreach my $bid (keys %{$batches{"assay"}}) {
	foreach my $fid (keys %{$batches{"assay"}->{"$bid"}->{"fluorophores"}}) {
		next unless 
		for (my $i=0; $i<scalar(@cols); $i++) {
			my $col = $cols[$i];
			print $cfh "\t" unless $i == 0;
			print $cfh "$batches{'assay'}->{$bid}->{'fluorophores'}->{$bid}->{$col}";
		}
		print $cfh "\n";
	}
}
close $cfh;
=cut


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
	$batch->{"$bid"}->{"cells"} = {};
	
	foreach my $key (keys %local_batch) {
		$batch->{"$bid"}->{"$key"} = "$local_batch{$key}";
		$bcols->{"$key"} = 1;
	}
	if (scalar keys %local_fluoro > 0) {
		$batch->{"$bid"}->{"fluorophores"} = {};
		foreach my $dye (keys %local_fluoro) {
			$batch->{"$bid"}->{"fluorophores"}->{$dye} = {};
			foreach my $key (keys %{$local_fluoro{$dye}}) {
				$batch->{"$bid"}->{"fluorophores"}->{$dye}->{"$key"} = "$local_fluoro{$dye}->{$key}";
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

	my (@sample_rows, @volstor_rows, @comment_rows);
	@sample_rows  = Spreadsheet::Read::rows($wkbk->[2]);
	@comment_rows = Spreadsheet::Read::rows($wkbk->[3]);
	@volstor_rows = Spreadsheet::Read::rows($wkbk->[4]) if defined Spreadsheet::Read::rows($wkbk->[4]);
	
	# all batch sheets except the Metadata sheet are based off sample_rows (ie the "Plate" sheet), so we only need to loop over that
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
				
				unless ("" eq "$sample_id" or "buffer" eq lc($sample_id)) {
					$batch->{"cells"}->{"$batch_loc"} = {"${btype}_location_in_batch" => "$batch_loc", "${btype}_batch_id" => "$bid"};
					# record the sample id
					$batch->{"cells"}->{"$batch_loc"}->{"sample_id"} = "$sample_id";

					# record the sample to id mapping
					my $lead = substr $btype, 0, 1;
					$s2id->{"$sample_id"} = {} unless defined $s2id->{"$sample_id"};
					my $num = 1;
					my $uid = "$sample_id.${lead}${num}";
					$s2id->{"$sample_id"}->{"${btype}_id"} = "$uid";

					# record any comment on the sample in this batch
					$batch->{"cells"}->{"$batch_loc"}->{"${btype}_comment"} = "$comment_rows[$i][$j]" if defined $comment_rows[$i][$j];
					if ("assay" eq "$btype") {
						# for assay batch, record any custom input volume
						if (defined $volstor_rows[$i][$j] and "$volstor_rows[$i][$j]" ne "") {
							$batch->{"cells"}->{"$batch_loc"}->{"${btype}_input_ul"} = "$volstor_rows[$i][$j]";
						} else {
							$batch->{"cells"}->{"$batch_loc"}->{"${btype}_input_ul"} = $batch->{"assay_input_ul"};
						}
					} elsif ("extraction" eq "$btype" or "archive" eq "$btype") {
						# for extraction and archive batches, record the storage location
						$batch->{"cells"}->{"$batch_loc"}->{"${btype}_location_in_storage"} = "$volstor_rows[$i][$j]" if defined $volstor_rows[$i][$j];
					}
				}
			}
		}	# j loop
	}	# i loop
}

exit 0;

