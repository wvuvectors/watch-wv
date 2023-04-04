#! /usr/bin/env perl


use strict;
use warnings;
use DateTime::Format::Excel;
use Data::Dumper;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;


my $usage = "\n";
$usage   .= "Usage: $progname\n";
$usage   .=   "Convert old format WaTCH db file to new format. Pass the old file on STDIN.\n";
$usage   .=   "A new directory CONVERSION will be created to hold the converted files.\n";
$usage   .=   "\n";

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
	}
}

# Define id keys for each table
my %table2key = (
	"sample" => "sample_id",
	"abatch" => "assay_batch_id",
	"archive" => "archive_id",
	"assay" => "assay_id",
	"cbatch" => "concentration_batch_id",
	"concentration" => "concentration_id",
	"control" => "control_id",
	"ebatch" => "extraction_batch_id",
	"extraction" => "extraction_id",
	"rbatch" => "archive_batch_id"
);

# Read in the conversion map
# ID map from old to new format, by table
my %new2old = ();
# Set up some constant vlaues to use in the new tables
my %predefs = (
	"sample_event"											 => "Routine Surveillance",
	"assay_batch_record_version" 				 => "1.0",
	"assay_quantification_type"  				 => "Direct Quantification",
	"assay_target_predicted_copies_per_ul_reaction" => 0, 
	"assay_analysis_software_version"		 => "1.2",
	"assay_reaction_ul"          				 => "20",
	"assay_target_macromolecule" 				 => "RNA",
	"concentration_batch_record_version" => "1.0",
	"control_input_ul"									 => "5.0",
	"extraction_batch_record_version"		 => "1.0",
	"extraction_eluant"									 => "Kit Buffer",
	"archive_batch_record_version"			 => "1.0",
	"archive_eluant"										 => "Zymo DNA/RNA Shield"
);
# Location for converted data, by table
my %tables = ();

open(my $fh, "<", "data/conversion_map.txt") or die "Unable to open map file: $!";
while (my $line = <$fh>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	
	my @cols = split "\t", "$line", -1;
	$new2old{$cols[1]} = {} unless defined $new2old{$cols[1]};
	$new2old{$cols[1]}->{$cols[2]} = "NA";
	$new2old{$cols[1]}->{$cols[2]} = "$cols[0]" unless "$cols[0]" eq "";

	$tables{"$cols[1]"} = {} unless defined $tables{"$cols[1]"};
}
close $fh;
#print Dumper(\%new2old);
#die;

# Read in the old WaTCHdb data
# Old data keyed by sample ID and replicate ID
my %watchdb  = ();
my $linenum  = 0;
my @colnames = ();

while (my $line = <>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	
	my @cols = split "\t", "$line", -1;
	if ($linenum == 0) {
		# First line of the file contains the column names
		# Store them and determine which column matches the key
		foreach (my $i=0; $i<scalar(@cols); $i++) {
			push @colnames, "$cols[$i]";
		}
	} else {
		my $thisId  = "$cols[0]";
		my $thisRep = "$cols[1]";
		$watchdb{"$thisId"} = {} unless defined $watchdb{"$thisId"};
		$watchdb{"$thisId"}->{"$thisRep"} = {} unless defined $watchdb{"$thisId"}->{"$thisRep"};
		for (my $i=2; $i<scalar(@cols); $i++) {
			$watchdb{"$thisId"}->{"$thisRep"}->{"$colnames[$i]"} = "$cols[$i]";
		}
	}
	$linenum++;
}
#print Dumper(\%watchdb);
#die;

my %table2newid = ("concentration" => {"code" => "c",  "incr" => 1},
									 "extraction"    => {"code" => "e",  "incr" => 1}, 
									 "assay"         => {"code" => "a",  "incr" => 1},
									 "control"       => {"code" => "t",  "incr" => 1},
									 "archive"       => {"code" => "r",  "incr" => 1});

# Fill in the new table hashes
foreach my $sample_id (keys %watchdb) {
	foreach my $rep_id (keys %{$watchdb{"$sample_id"}}) {
		my %entry = %{$watchdb{"$sample_id"}->{"$rep_id"}};

		# sample table is easiest
		$tables{"sample"} = {} unless defined $tables{"sample"};
		$tables{"sample"}->{"$sample_id"} = {} unless defined $tables{"sample"}->{"$sample_id"};
		foreach my $new_id (keys %{$new2old{"sample"}}) {
			if (defined $predefs{"$new_id"}) {
				$tables{"sample"}->{"$sample_id"}->{"$new_id"} = $predefs{"$new_id"};
			} elsif ("$new_id" eq "sample_id") {
				$tables{"sample"}->{"$sample_id"}->{"$new_id"} = "$sample_id";
			} else {
				my $old_id = $new2old{"sample"}->{"$new_id"};
				$tables{"sample"}->{"$sample_id"}->{"$new_id"} = $entry{"$old_id"};
			}
		}
		
		# Handle the batches
		foreach my $table_name (keys %tables) {
			next unless "$table_name" =~ /batch/;
			
			my $id_col_new = $table2key{"$table_name"};
			my $id_col_old = $new2old{"$table_name"}->{"$id_col_new"};

			my $uid = $entry{"$id_col_old"};
			die "FATAL: No $table_name id exists for $sample_id!\n" unless defined $uid;

			# can skip if the batch already exists
			next if $table_name =~ /batch/ and defined $tables{"$table_name"}->{"$uid"};
			
#			if (defined $entry{"$id_col_old"}) {
#				$uid = $entry{"$id_col_old"};
#			} else {
#				$uid = "${sample_id}." . "$table2newid{$table_name}->{'code'}" . "$table2newid{$table_name}->{'incr'}";
#				$table2newid{$table_name}->{"incr"} = $table2newid{$table_name}->{"incr"} + 1;
#			}
			
			$tables{"$table_name"}->{"$uid"} = {};
			foreach my $new_id (keys %{$new2old{"$table_name"}}) {
				if (defined $predefs{"$new_id"}) {
					$tables{"$table_name"}->{"$uid"}->{"$new_id"} = $predefs{"$new_id"};
				} else {
					my $old_id = $new2old{"$table_name"}->{"$new_id"};
					$tables{"$table_name"}->{"$uid"}->{"$new_id"} = $entry{"$old_id"};
				}
			}
		}
	}
}
print Dumper(\%tables);
die;

# Write the tables to file
mkdir("data/converted/") unless -d "data/converted/";

foreach my $table_name (keys %new2old) {
	my @hdr = sort {$a cmp $b} keys %{$new2old{"$table_name"}};
	my $id_col = $table2key{"$table_name"};
	
	open(my $ofh, ">", "data/converted/watchdb.$table_name.txt") or die "Unable to open $table_name file: $!";
	print $ofh "$id_col";
	for (my $i=0; $i<scalar(@hdr); $i++) {
		print $ofh "\t$hdr[$i]" unless "$hdr[$i]" eq "$id_col";
	}
	print $ofh "\n";

	foreach my $uid (keys %{$tables{"$table_name"}}) {
		print $ofh "$uid";
		for (my $i=0; $i<scalar(@hdr); $i++) {
			print $ofh "\t$tables{$table_name}->{$uid}->{$hdr[$i]}" unless "$hdr[$i]" eq "$id_col";
		}
		print $ofh "\n";
	}
	close $ofh;
}


exit 0;

