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
$usage   .= "Usage: $progname [options] -i RUNDIR -o DBDIR\n";
$usage   .=   "Query RUNDIR for any unprocessed batch files in subfolders '0 SAMPLES', '1 CB_PLATES', '2 EB_PLATES', '3 AB_PLATES,' '4 AB_RESULTS,' and '5 RB_PLATES.'\n";
$usage   .=   "Any batch files in these subfolders that do not appear in the latest database will be printed as file paths to STDOUT.\n";
$usage   .=   "The list of batch IDs that already exist is in DBDIR/latest/watchdb.completed_batches.txt.\n";
$usage   .=   "\n";

my $rundir;
my $outdir;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } elsif ($arg eq "-i") {
		$rundir = shift;
	} elsif ($arg eq "-o") {
		$outdir = shift;
	}
}

die "FATAL: $progname requires a valid run directory.\n$usage\n" unless defined $rundir and -d "$rundir";

# First read in ids for the batches we've already completed.

my $barchive_f = "$outdir/latest/watchdb.completed_batches.txt";
my %completed    = ();

if (-f "$barchive_f") {

	open(my $achFH, "<", "$barchive_f") or die "Unable to read list of completed batches from $barchive_f: $!";
	my $count = 1;
	while (my $line = <$achFH>) {
		chomp $line;
		$completed{"$line"} = 1;
	}
	close $achFH;
} else {
	warn "\nNo file found listing completed batches! Assuming there have been no batches processed for this run.\n";
}

#print Dumper(\%completed);
#die;

# Now identify batch files from the "run directory" sub-folders that we have not completed.
my @files2proc = ();
#my @subfolders = ("3 AB_PLATES");
my @subfolders = ("1 CB_PLATES", "2 EB_PLATES", "3 AB_PLATES", "5 RB_PLATES");
foreach my $subf (@subfolders) {
	#print "$rundir/$subf\n";
	next unless -d "$rundir/$subf/";
	opendir(my $dirBatch, "$rundir/$subf/") or die "Although I found it, I am unable to open $rundir/$subf/: $!";
	my @batch_files  = grep { (/\.xlsx$/) && (!/^~/) && -f "$rundir/$subf/$_" } readdir($dirBatch);
	closedir $dirBatch;
	#print Dumper(\@batch_files);
	foreach my $f (@batch_files) {
		#print "$f\n";
		my $batch_wkbk = ReadData("$rundir/$subf/$f", dtfmt => "mm/dd/yy");

		# Get the batch type and id, which are recorded in cells B1 & B6 (sheet 1, rows 1 & 6, zero-based column 1) of each Excel file.
		#my @mdrows  = Spreadsheet::Read::rows($batch_wkbk->[1]);
		my @typerow = Spreadsheet::Read::row($batch_wkbk->[1], 1);
		my $btype   = lc $typerow[1];
		my @idrow = Spreadsheet::Read::row($batch_wkbk->[1], 6);
		my $bid   = $idrow[1];
		#print "$btype\t$bid\n";
		unless (defined $completed{"$bid"}) {
			#print "$btype\t$bid\n";
			if ("assay" eq "$btype") {
				# For assay batches, we first confirm the assay data file exists. 
				# These are in "4 AB_RESULTS" subfolder and begin with the assay batch id and end in csv.
				# If it doesn't exist we want to skip this assay batch file for now.
				my $resultsf = "4 AB_RESULTS";
				opendir(my $dirResults, "$rundir/$resultsf/") or die "Although I found it, I am unable to open $rundir/$resultsf/: $!";
				my @result_files  = grep { (/^${bid}_.+\.csv$/) && (!/^~/) && -f "$rundir/$resultsf/$_" } readdir($dirResults);
				if (scalar(@result_files) > 0) {
					push(@files2proc, "$rundir/$subf/$f");
					foreach my $rf (@result_files) {
						push(@files2proc, "$rundir/$resultsf/$rf");
					}
				} else {
					die "No results files found in $rundir/$resultsf/";
				}
			} else {
				push(@files2proc, "$rundir/$subf/$f");
			}
		}
	}
}
#print Dumper(\@files2proc);
#die;

# Handle the Asset Tiger file that contains the sample metadata.
# Only do this if there is at least one batch file to process.
unless (scalar(@files2proc) == 0) {
	my $subf = "0 SAMPLES";
	my $f = "AssetTagReport.csv";
	opendir(my $dirResults, "$rundir/$subf/") or die "Although I found it, I am unable to open $rundir/$subf/: $!";
	my @at_files  = grep { (/^AssetTagReport\.csv$/) && (!/^~/) && -f "$rundir/$subf/$_" } readdir($dirResults);
	if (scalar(@at_files) > 0) {
		push(@files2proc, "$rundir/$subf/$f");
	}
}


# Print the list of files to STDOUT.
foreach my $uf (@files2proc) {
	print "$uf\n";
}

exit 0;

