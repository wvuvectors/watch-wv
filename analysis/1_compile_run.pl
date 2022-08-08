#! /usr/bin/env perl

use strict;
use warnings;
use Text::CSV;
#use DateTime qw( );
use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use Data::Dumper;


# read in plate files from directory and convert to single table

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=   "Merge WaTCH lab excel files into a single table.\n";
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

opendir(my $dirH, $rundir) or die "Unable to open directory $rundir: $!";
my @plate_files = grep { (/\.xlsx$/) && (!/^~/) && -f "$rundir/$_" } readdir($dirH);
closedir $dirH;

my %plates       = ("assay" => {}, "concentration" => {}, "extraction" => {});
my %platecols    = ("assay" => {}, "concentration" => {}, "extraction" => {});
my %sample2id    = ();

my %updatecols   = (
	"assay" => ["assay_id", "extraction_id", "assay_plate_id", "assay_plate_cell", "assay_input_ul", "assay_type", "assay_target", "assay_target_category", "assay_target_genetic_locus", "assay_target_macromolecule", "assay_fluorophore", "assay_accepted_droplets", "assay_calculated_copies_per_ul_reaction", "assay_copies_per_ul_reaction", "assay_comments"], 
	"concentration" => ["concentration_id", "sample_id", "concentration_plate_id", "concentration_plate_cell", "concentration_storage_location", "concentration_comments"], 
	"extraction" => ["extraction_id", "concentration_id", "extraction_plate_id", "extraction_plate_cell", "extraction_storage_location", "extraction_comments"]
);

my %fluorcols = ();


# read in assay data file, assay_plate.csv
# Well,Sample description 1,Sample description 2,Sample description 3,Sample description 4,Target,Conc(copies/uL),Status,Experiment,SampleType,TargetType,Supermix,DyeName(s),Accepted Droplets,Positives,Negatives
my %cell2data = ();
my $csvA = Text::CSV->new({auto_diag => 4, binary => 1});
open (my $afh, "<", "$rundir/assay_plate.csv") or die "Unable to open $rundir/assay_plate.csv for reading: $!\n";
my $count = 0;
while (my $line = $csvA->getline($afh)) {
	my @cols = @$line;
	$count++;
	if (scalar @cols < 15) {
		warn "WARN  : Line $count in assay_plate.csv has " . scalar(@cols) . " but requires exactly 15. Skipping.\n";
		next;
	}
	my ($cell, $result, $dye, $accepted_droplets) = ($cols[0], $cols[6], $cols[12], $cols[13]);
	$cell2data{"$cell"} = {} unless defined $cell2data{"$cell"};
	$cell2data{"$cell"}->{"$dye"} = {"assay_copies_per_ul_reaction" => $result, "assay_accepted_droplets" => $accepted_droplets};
}
close $afh;


foreach (@plate_files) {
	my $plate_wkbk = ReadData("$rundir/$_", dtfmt => "mm/dd/yy");
	
	# parse data based on plate type, which is recorded in cell B1 (page 1, row 1, column 2)
	my @mdrows = Spreadsheet::Read::rows($plate_wkbk->[1]);
	my @typerow = Spreadsheet::Read::row($plate_wkbk->[1], 1);
	my $ptype = lc $typerow[1];

	my $plate_id = read_plate_metadata("$ptype", \@mdrows, $plates{"$ptype"}, $platecols{"$ptype"});
	
	# Now read the samples in each plate
	read_plate_data("$ptype", $plate_id, $plate_wkbk, $plates{"$ptype"}->{"$plate_id"}, \%sample2id);
}

#print Dumper(\%cell2data);
#print Dumper(\%plates);
#print Dumper(\%sample2id);
#die;


# Write plate metadata to update files
my %platef = ("assay" => "aplate", "concentration" => "cplate", "extraction" => "eplate");
foreach my $ptype (keys %plates) {
	my $pname = $platef{"$ptype"};
	open (my $pfh, ">", "$rundir/update.$pname.txt") or die "Unable to open $rundir/update.$pname.txt for writing: $!";
	my @cols = sort keys %{$platecols{"$ptype"}};
	print $pfh join("\t", @cols) . "\n";
	foreach my $pid (keys %{$plates{"$ptype"}}) {
		for (my $i=0; $i<scalar(@cols); $i++) {
			my $col = $cols[$i];
			print $pfh "\t" unless $i == 0;
			print $pfh "$plates{$ptype}->{$pid}->{$col}";
		}
		print $pfh "\n";
	}
	close $pfh;
}

# write one-to-one update files (conc, extract)
foreach my $ptype (keys %plates) {
	next if "$ptype" eq "assay";
	#my $lead = substr $ptype, 0, 1;
	open (my $ufh, ">", "$rundir/update.$ptype.txt") or die "Unable to open $rundir/update.$ptype.txt for writing: $!";
	print $ufh join("\t", @{$updatecols{"$ptype"}}) . "\n";
	foreach my $pid (keys %{$plates{"$ptype"}}) {
		foreach my $cell (keys %{$plates{"$ptype"}->{"$pid"}->{"cells"}}) {
			my $sample_id = $plates{"$ptype"}->{"$pid"}->{"cells"}->{"$cell"}->{"sample_id"};
			my $uid = $sample2id{"$sample_id"}->{"${ptype}_id"};
			print $ufh "$uid";
			for (my $i=1; $i<scalar(@{$updatecols{"$ptype"}}); $i++) {
				my $key = $updatecols{"$ptype"}->[$i];
				if (defined $plates{"$ptype"}->{"$pid"}->{"cells"}->{"$cell"}->{"$key"}) {
					print $ufh "\t" . $plates{"$ptype"}->{"$pid"}->{"cells"}->{"$cell"}->{"$key"};
				} elsif (defined $plates{"$ptype"}->{"$pid"}->{"$key"}) {
					print $ufh "\t" . $plates{"$ptype"}->{"$pid"}->{"$key"};
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
# assay_fluorophore
# assay_accepted_droplets
# assay_calculated_copies_per_reaction
# assay_copies_per_reaction
# assay_result
# assay_comments
my $ptype = "assay";
open (my $ufh, ">", "$rundir/update.$ptype.txt") or die "Unable to open $rundir/update.$ptype.txt for writing: $!";
print $ufh join("\t", @{$updatecols{"$ptype"}}) . "\n";
my $lead = substr $ptype, 0, 1;
foreach my $pid (keys %{$plates{"$ptype"}}) {
	foreach my $cell (keys %{$plates{"$ptype"}->{"$pid"}->{"cells"}}) {
		my $sample_id = $plates{"$ptype"}->{"$pid"}->{"cells"}->{"$cell"}->{"sample_id"};
		
		# REMOVE CONTROLS
		next if "$sample_id" eq "PCR NC" or "$sample_id" eq "PCR PC";
		
		my $count = 1;
		foreach my $fluor (keys %{$plates{"$ptype"}->{"$pid"}->{"fluorophores"}}) {
			my $uid = "$sample_id.${lead}$count";
			$count++;
			print $ufh "$uid";
			for (my $i=1; $i<scalar(@{$updatecols{"$ptype"}}); $i++) {
				my $key = $updatecols{"$ptype"}->[$i];
				if (defined $plates{"$ptype"}->{"$pid"}->{"cells"}->{"$cell"}->{"$key"}) {
					print $ufh "\t" . $plates{"$ptype"}->{"$pid"}->{"cells"}->{"$cell"}->{"$key"};
				} elsif (defined $plates{"$ptype"}->{"$pid"}->{"$key"}) {
					print $ufh "\t" . $plates{"$ptype"}->{"$pid"}->{"$key"};
				} elsif (defined $sample2id{"$sample_id"}->{"$key"}) {
					print $ufh "\t" . $sample2id{"$sample_id"}->{"$key"};
				} elsif (defined $plates{"$ptype"}->{"$pid"}->{"fluorophores"}->{"$fluor"}->{"$key"}) {
					print $ufh "\t" . $plates{"$ptype"}->{"$pid"}->{"fluorophores"}->{"$fluor"}->{"$key"};
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
# write plate assay controls to update file
open (my $cfh, ">", "$rundir/update.controls.txt") or die "Unable to open $rundir/update.controls.txt for writing: $!";
my @cols = sort keys %fluorcols;
print $cfh join("\t", @cols) . "\n";
foreach my $pid (keys %{$plates{"assay"}}) {
	foreach my $fid (keys %{$plates{"assay"}->{"$pid"}->{"fluorophores"}}) {
		next unless 
		for (my $i=0; $i<scalar(@cols); $i++) {
			my $col = $cols[$i];
			print $cfh "\t" unless $i == 0;
			print $cfh "$plates{'assay'}->{$pid}->{'fluorophores'}->{$pid}->{$col}";
		}
		print $cfh "\n";
	}
}
close $cfh;
=cut


sub read_plate_metadata {
	my $ptype = shift;
	my $sheet = shift;
	my $plate = shift;
	my $pcols = shift;
	
	$ptype = lc $ptype;
	
	my $pid;
	my @dyes = ();
	
	my %local_plate  = ();
	my %local_fluoro = ();
	
	for (my $i=1; $i < scalar(@$sheet); $i++) {
		next if "" eq "$sheet->[$i][0]";
		my $key = lc "$sheet->[$i][0]";
		my $val = "$sheet->[$i][1]";

		push @dyes, "$val" if "fluorophore" eq "$key";

		if ("$key" eq "plate id") {
			$pid = "$val";
		} elsif ("$key" eq "date") {
			my $datetime = DateTime::Format::Excel->parse_datetime($val);
			$val = $datetime->ymd('-');
		}
		$key =~ s/ /_/gi;
		$key = "${ptype}_$key";
		if (scalar @dyes > 0) {
			$local_fluoro{$dyes[-1]} = {} unless defined $local_fluoro{$dyes[-1]};
			$local_fluoro{$dyes[-1]}->{"$key"} = "$val";
		} else {
			$local_plate{"$key"} = "$val";
		}
	}
	
#	if (scalar keys %local_fluoro > 0) {
#		print Dumper(\%local_fluoro);
#		die;
#	}

	die "FATAL: No plate id recognized: $ptype\n" unless defined $pid;
	
	$plate->{"$pid"} = {} unless defined $plate->{"$pid"};
	
	# placeholder to hold sample data (keyed per well)	
	$plate->{"$pid"}->{"cells"} = {};
	
	foreach my $key (keys %local_plate) {
		$plate->{"$pid"}->{"$key"} = "$local_plate{$key}";
		$pcols->{"$key"} = 1;
	}
	if (scalar keys %local_fluoro > 0) {
		$plate->{"$pid"}->{"fluorophores"} = {};
		foreach my $dye (keys %local_fluoro) {
			$plate->{"$pid"}->{"fluorophores"}->{$dye} = {};
			foreach my $key (keys %{$local_fluoro{$dye}}) {
				$plate->{"$pid"}->{"fluorophores"}->{$dye}->{"$key"} = "$local_fluoro{$dye}->{$key}";
			}
		}
	}
	
	return $pid;
}	

sub read_plate_data {
	my $ptype = shift;
	my $pid   = shift;
	my $wkbk  = shift;
	my $plate = shift;
	my $s2id  = shift;
	
	$ptype = lc $ptype;
	
	my @row2name = ("", "A", "B", "C", "D", "E", "F", "G", "H");

	my (@sample_rows, @volstor_rows, @comment_rows);
	@sample_rows  = Spreadsheet::Read::rows($wkbk->[2]);
	@volstor_rows = Spreadsheet::Read::rows($wkbk->[3]);
	@comment_rows = Spreadsheet::Read::rows($wkbk->[4]);
	
	# all other plate sheets are based off the sample_rows sheet, so we only need to loop over that
	for (my $i=1; $i < scalar(@sample_rows); $i++) {
		for (my $j=1; $j < scalar(@{$sample_rows[$i]}); $j++) {
			if (!defined $sample_rows[$i][$j]) {
				warn "****WARNING**** $ptype plate has no sample ID at row $i, column $j.";
				next;
			} else {
				my $sample_id = "$sample_rows[$i][$j]";
				$sample_id =~ s/^\s+//;
				$sample_id =~ s/\s+$//;
				
				my $plate_cell = "$row2name[$i]" . sprintf("%02d", $j);
				
				unless ("" eq "$sample_id" or "buffer" eq lc($sample_id)) {
					$plate->{"cells"}->{"$plate_cell"} = {"${ptype}_plate_cell" => "$plate_cell", "${ptype}_plate_id" => "$pid"};
					# record the sample id
					$plate->{"cells"}->{"$plate_cell"}->{"sample_id"} = "$sample_id";

					# record the sample to id mapping
					my $lead = substr $ptype, 0, 1;
					$s2id->{"$sample_id"} = {} unless defined $s2id->{"$sample_id"};
					my $num = 1;
					my $uid = "$sample_id.${lead}${num}";
					$s2id->{"$sample_id"}->{"${ptype}_id"} = "$uid";

					# record any comment on the sample in this plate
					$plate->{"cells"}->{"$plate_cell"}->{"${ptype}_comments"} = "$comment_rows[$i][$j]" if defined $comment_rows[$i][$j];
					if ("assay" eq "$ptype") {
						# for assay plates, store any custom input volume
						if (defined $volstor_rows[$i][$j] and "$volstor_rows[$i][$j]" ne "") {
							$plate->{"cells"}->{"$plate_cell"}->{"${ptype}_input_ml"} = "$volstor_rows[$i][$j]";
						} else {
							$plate->{"cells"}->{"$plate_cell"}->{"${ptype}_input_ul"} = $plate->{"assay_input_ul"};
						}
					} else {
						# for all other plates, store the storage location
						$plate->{"cells"}->{"$plate_cell"}->{"${ptype}_storage_location"} = "$volstor_rows[$i][$j]" if defined $volstor_rows[$i][$j];
					}
				}
			}
		}	# j loop
	}	# i loop
}

exit 0;

