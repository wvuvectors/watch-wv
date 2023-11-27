#! /usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

#print "******\n";
#print "Running $progname.\n";
#print "******\n";


my $usage = "\n";
$usage   .= "Usage: $progname [options] RUNDIR\n";
$usage   .=  "Validate the compiled run in RUNDIR vs the main watch database. Particularly look for duplication of table IDs.\n";
$usage   .=   "\n";

my $dbdir  = "data/latest";
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



# Fetch ids from the current WaTCH database tables.
#
print "Validating ids from existing WaTCH database tables...\n";
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
#print "Finished reading " . scalar(keys %watch2id) . " WaTCH tables.\n";
print "Done.\n\n";

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
print "Validating ids from run update tables...\n";
foreach my $table (keys %runup2id) {
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
print "Done.\n\n";

print "-------------------------------\n";
print "Run validation results ($rundir):\n";
print sprintf('%-13s', "TABLE") . "\t" . sprintf('%-13s', "UPDATE IDS") . "\t" . sprintf('%-13s', "DUPLICATES") . "\t" . sprintf('%-13s', "COLLISIONS") . "\n";
foreach my $table (keys %runup2id) {
	print sprintf('%-13s', "$table") . "\t" . sprintf('%-13s', scalar(keys %{$runup2id{"$table"}})) . "\t" . sprintf('%-13s', scalar(keys %{$runupdup{"$table"}})) . "\t" . sprintf('%-13s', scalar(keys %{$collisions{"$table"}})) . "\n";
}
print "-------------------------------\n\n";


my $status = 0;
my ($num_dups_watch, $num_dups_runup, $num_collisions) = (0,0,0);
foreach my $table (keys %table2key) {
	$num_dups_watch = $num_dups_watch + scalar(keys %{$watchdup{"$table"}});
	$num_dups_runup = $num_dups_runup + scalar(keys %{$runupdup{"$table"}});
	$num_collisions = $num_collisions + scalar(keys %{$collisions{"$table"}});
}
$status = $num_dups_watch + $num_dups_runup + $num_collisions;
#if ($status > 0) {
	open (my $v1FH, ">", "$rundir/updates/_collisions.txt");
	print $v1FH "The following IDs are present in both the existing WaTCH database and the current run.\n\nTABLE\tID\n";
	foreach my $table (keys %collisions) {
		foreach my $id (keys %{$collisions{"$table"}}) {
			print $v1FH "$table\t$id\n";
		}
	}
	close $v1FH;
	open (my $v2FH, ">", "$rundir/updates/_rundups.txt");
	print $v2FH "The following IDs are duplicated in the current run data.\n\nTABLE\tID\n";
	foreach my $table (keys %runupdup) {
		foreach my $id (keys %{$runupdup{"$table"}}) {
			print $v2FH "$table\t$id\n";
		}
	}
	close $v2FH;
	open (my $v3FH, ">", "$rundir/updates/_watchdups.txt");
	print $v3FH "The following IDs are duplicated in the WaTCH database (before appending any data from the current run).\n\nTABLE\tID\n";
	foreach my $table (keys %watchdup) {
		foreach my $id (keys %{$watchdup{"$table"}}) {
			print $v3FH "$table\t$id\n";
		}
	}
	close $v3FH;
#}

#print "******\n";
#print "Finished $progname.\n";
#print "******\n";


exit $status;


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


