#! /usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

print "******\n";
print "Running $progname.\n";
print "******\n";


my $usage = "\n";
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=  "Validate the compiled run in RUNDIR vs the main watch database. Particularly look for duplication of table IDs.\n";
$usage   .=   "\n";

my $dbdir  = "data/watchdb/LATEST";
my $buffer = "             "; # 13 spaces to buffer short table names in print statements
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

# WaTCH key columns
my %table2key = ("abatch"        => "assay_batch_id", 
								 "archive"       => "archive_id", 
								 "assay"         => "assay_id", 
								 "cbatch"        => "concentration_batch_id", 
								 "concentration" => "concentration_id", 
								 "control"       => "control_id", 
								 "ebatch"        => "extraction_batch_id", 
								 "extraction"    => "extraction_id", 
								 "rbatch"        => "archive_batch_id");

# Existing WaTCH ids keyed by WaTCH table
my %watch2id = ("abatch"        => {}, 
								"archive"       => {}, 
								"assay"         => {}, 
								"cbatch"        => {}, 
								"concentration" => {}, 
								"control"       => {}, 
								"ebatch"        => {}, 
								"extraction"    => {}, 
								"rbatch"        => {});

# Count lines in each WaTCH table file
my %watch2lines = ("abatch"        => 0, 
									 "archive"       => 0, 
									 "assay"         => 0, 
									 "cbatch"        => 0, 
									 "concentration" => 0, 
									 "control"       => 0, 
									 "ebatch"        => 0, 
									 "extraction"    => 0, 
									 "rbatch"        => 0);

# Track duplicated IDs in the WaTCH database
my %watchdup = ("abatch"        => {}, 
								"archive"       => {}, 
								"assay"         => {}, 
								"cbatch"        => {}, 
								"concentration" => {}, 
								"control"       => {}, 
								"ebatch"        => {}, 
								"extraction"    => {}, 
								"rbatch"        => {});


# Run update IDs keyed by WaTCH table
my %runup2id = ("abatch"        => {}, 
								"archive"       => {}, 
								"assay"         => {}, 
								"cbatch"        => {}, 
								"concentration" => {}, 
								"control"       => {}, 
								"ebatch"        => {}, 
								"extraction"    => {}, 
								"rbatch"        => {});

# Count lines in each update file
my %runup2lines = ("abatch"        => 0, 
									 "archive"       => 0, 
									 "assay"         => 0, 
									 "cbatch"        => 0, 
									 "concentration" => 0, 
									 "control"       => 0, 
									 "ebatch"        => 0, 
									 "extraction"    => 0, 
									 "rbatch"        => 0);

# Track duplicated IDs in the update tables
my %runupdup = ("abatch"        => {}, 
								"archive"       => {}, 
								"assay"         => {}, 
								"cbatch"        => {}, 
								"concentration" => {}, 
								"control"       => {}, 
								"ebatch"        => {}, 
								"extraction"    => {}, 
								"rbatch"        => {});


# Track ID collisions between update and the WaTCH database
my %collisions = ("abatch"        => {}, 
									"archive"       => {}, 
									"assay"         => {}, 
									"cbatch"        => {}, 
									"concentration" => {}, 
									"control"       => {}, 
									"ebatch"        => {}, 
									"extraction"    => {}, 
									"rbatch"        => {});


print "Validating run data from $rundir/updates/ folder against WaTCH database tables in $dbdir.\n";

# Fetch ids from the current WaTCH database tables.
#
my $idtotal = 0;
foreach my $table (keys %watch2id) {
	print "Retrieving ids from WaTCH table $table\n";
	my $keyname = $table2key{"$table"};
	my $keycol  = -1;
	my $linenum = 0;
	open (my $dbFH, "<", "$dbdir/watchdb.${table}.txt") or die "Unable to open $dbdir/watchdb.${table}.txt for reading: $!\n";
	while (my $line = <$dbFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		my @cols = split "\t", "$line", -1;
		if ($linenum == 0) {
			# First line of the file contains the column names
			# Determine which column matches the key and contains the table ids
			FIND:
			foreach (my $i=0; $i<scalar(@cols); $i++) {
				if ("$keyname" eq "$cols[$i]") {
					$keycol = $i;
					last FIND;
				}
			}
		} else {
			# Extract the id for this row
			my $thisId = "$cols[$keycol]";
			if (defined $watch2id{"$table"}->{"$thisId"}) {
				# If the id already exists in the id hash, mark it as a duplicate
				$watchdup{"$table"}->{"$thisId"} = 0 unless defined $watchdup{"$table"}->{"$thisId"};
				$watchdup{"$table"}->{"$thisId"} = $watchdup{"$table"}->{"$thisId"} + 1;
			} else {
				# If it is a new id, add it to the main id hash
				$watch2id{"$table"}->{"$cols[$keycol]"} = 1;
			}
		}
		$linenum++;
	}
	$watch2lines{"$table"} = $linenum;
	close $dbFH;
}
print "Finished reading " . scalar(keys %watch2id) . " WaTCH tables.\n";

print "-------------------------------\n";
print "WaTCH internal validation results ($dbdir):\n";
print sprintf('%-13s', "TABLE") . "\t" . sprintf('%-13s', "WaTCH IDS") . "\t" . sprintf('%-13s', "DUPLICATES") . "\n";
foreach my $table (keys %watch2id) {
	print sprintf('%-13s', "$table") . "\t" . sprintf('%-13s', scalar(keys %{$watch2id{"$table"}})) . "\t" . sprintf('%-13s', scalar(keys %{$watchdup{"$table"}})) . "\n";
}
print "-------------------------------\n\n";


#print Dumper(\%table2id);
#die;


# Scan update tables for duplicate IDs and flag.
#
foreach my $table (keys %runup2id) {
	print "Validating ids from update table $table\n";
	my $keyname = $table2key{"$table"};
	my $keycol  = -1;
	my $linenum = 0;
	open (my $dbFH, "<", "$rundir/updates/update.${table}.txt") or die "Unable to open $rundir/updates/update.${table}.txt for reading: $!\n";
	while (my $line = <$dbFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		my @cols = split "\t", "$line", -1;
		if ($linenum == 0) {
			# First line of the file contains the column names
			# Determine which column matches the key and contains the table ids
			FIND:
			foreach (my $i=0; $i<scalar(@cols); $i++) {
				if ("$keyname" eq "$cols[$i]") {
					$keycol = $i;
					last FIND;
				}
			}
		} else {
			# Extract the id for this row
			my $thisId = "$cols[$keycol]";
			if (defined $runup2id{"$table"}->{"$thisId"}) {
				# If the id already exists in the id hash, mark it as a duplicate
				# No need to compare against WaTCH since that was done on the first instance
				$runupdup{"$table"}->{"$thisId"} = 0 unless defined $runupdup{"$table"}->{"$thisId"};
				$runupdup{"$table"}->{"$thisId"} = $runupdup{"$table"}->{"$thisId"} + 1;
			} else {
				# If it is a new id, add it to the main id hash
				$runup2id{"$table"}->{"$cols[$keycol]"} = 1;
				if (defined $watch2id{"$table"}->{"$thisId"}) {
					# If the id already exists in the watch id hash, mark it as a collision
					$collisions{"$table"}->{"$thisId"} = 1;
				}
			}
		}
		$linenum++;
	}
	$runup2lines{"$table"} = $linenum;
	close $dbFH;
}

print "-------------------------------\n";
print "Update validation results ($rundir):\n";
print sprintf('%-13s', "TABLE") . "\t" . sprintf('%-13s', "UPDATE IDS") . "\t" . sprintf('%-13s', "DUPLICATES") . "\t" . sprintf('%-13s', "COLLISIONS") . "\n";
foreach my $table (keys %runup2id) {
	print sprintf('%-13s', "$table") . "\t" . sprintf('%-13s', scalar(keys %{$runup2id{"$table"}})) . "\t" . sprintf('%-13s', scalar(keys %{$runupdup{"$table"}})) . "\t" . sprintf('%-13s', scalar(keys %{$collisions{"$table"}})) . "\n";
}
print "-------------------------------\n\n";

print "******\n";
print "Finished $progname.\n";
print "******\n";


my $status = 0;
my ($num_dups_watch, $num_dups_runup, $num_collisions) = (0,0,0);
foreach my $table (keys %table2key) {
	$num_dups_watch = $num_dups_watch + scalar(keys %{$watchdup{"$table"}});
	$num_dups_runup = $num_dups_runup + scalar(keys %{$runupdup{"$table"}});
	$num_collisions = $num_collisions + scalar(keys %{$collisions{"$table"}});
}
$status = $num_dups_watch + $num_dups_runup + $num_collisions;
if ($status > 0) {
	# print dup ids and collisions!
}
exit $status;


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
=cut


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


