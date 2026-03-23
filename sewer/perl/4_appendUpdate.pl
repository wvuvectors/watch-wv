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
$usage   .= "Usage: $progname [options] DBDIR UPDATEDIR\n";
$usage   .=  "Append the update in UPDATEDIR to the database in DBDIR.\n";
$usage   .=   "\n";

my ($dbdir, $updir) = ($ARGV[0], $ARGV[1]);

=cut
while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } elsif (defined $dbdir) {
		$updir = $arg;
	} else {
		$dbdir = $arg;
	}
}
=cut

die "FATAL: $progname requires a valid database directory ($dbdir).\n$usage\n" unless defined $dbdir and -d "$dbdir";
die "FATAL: $progname requires a valid update directory ($updir).\n$usage\n" unless defined $updir and -d "$updir";

my %orderedcols = (
	"abatch"        => ["assay_batch_id", "assay_date", "assay_reaction_ul", "assay_machine", "assay_amplification_method", "assay_amplification_method_lot_id", "assay_quantification_method", "assay_quantification_type", "assay_batch_record_version", "assay_method", "assay_method_lot_id", "assay_qx_manager_version", "assay_run_by", "assay_batch_comment"],
	"archive"       => ["archive_id", "sample_id", "archive_batch_id", "archive_location_in_batch", "archive_location_in_storage", "archive_comment"],
	"assay"         => ["assay_id", "extraction_id", "sample_id", "assay_batch_id", "assay_location_in_batch", "assay_input_ul", "assay_class", "assay_type", "assay_target", "assay_target_genetic_locus", "assay_template", "assay_target_macromolecule", "assay_target_fluorophore", "assay_accepted_droplets", "assay_target_predicted_copies_per_ul_reaction", "assay_target_copies_per_ul_reaction", "assay_comment"], 
	"cbatch"        => ["concentration_batch_id", "concentration_date", "concentration_input_ml", "concentration_machine", "concentration_method", "concentration_method_lot_id", "concentration_output_ml", "concentration_run_by", "concentration_batch_record_version", "concentration_batch_comment"],
	"concentration" => ["concentration_id", "sample_id", "concentration_batch_id", "concentration_location_in_batch", "concentration_comment"], 
	"ebatch"        => ["extraction_batch_id", "extraction_date", "extraction_input_ul", "extraction_eluant", "extraction_machine", "extraction_method", "extraction_method_lot_id", "extraction_output_ul", "extraction_batch_record_version", "extraction_run_by", "extraction_batch_comment"],
	"extraction"    => ["extraction_id", "concentration_id", "extraction_batch_id", "extraction_location_in_batch", "extraction_location_in_storage", "extraction_comment"],
	"rbatch"        => ["archive_batch_id", "archive_date", "archive_input_ml", "archive_eluant", "archive_machine", "archive_method", "archive_method_lot_id", "archive_output_ml", "archive_run_by", "archive_batch_record_version", "archive_batch_comment"],
	"sample"        => ["sample_id", "sample_status", "location_id", "sample_event", "sample_qc", "sample_collection_start_datetime", "sample_collection_end_datetime", "sample_recovered_datetime", "sample_collection_by", "sample_flow", "sample_received_by", "sample_received_date", "sample_ph_lab", "sample_comment"]
);

my %watch = ();

# Put the current data (before update) into %watch, keyed by table.
#
opendir(my $dbdFH, "$dbdir") or die "Although I found it, I am unable to open $dbdir: $!";
my @dbfiles  = grep { (/^watchdb\..+\.txt$/) && (!/^~/) && -f "$dbdir/$_" } readdir($dbdFH);
for my $f (@dbfiles) {
	#print "$_\n";
	my $table = "$f";
	$table =~ s/watchdb\.(.+?)\.txt/$1/gi;
	next unless defined $table and defined $orderedcols{"$table"};

	$watch{"$table"} = {} unless defined $watch{"$table"};
	#print "$table\n";
	my @colheaders = ();
	my $linenum = 1;
	my $idcol;

	open(my $fh, "<", "$dbdir/$f") or die "FATAL: Unable to open $dbdir/$f for reading: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		my @cols = split /\t/, "$line", -1;
		if ($linenum == 1) {
			$linenum++;
			for (my $i=0; $i<scalar(@cols); $i++) {
				my $val = $cols[$i];
				push @colheaders, "$val";
				$idcol = $i if "$val" eq "$orderedcols{$table}->[0]";
			}
		} else {
			my $uid = "$cols[$idcol]";
			$watch{"$table"}->{"$uid"} = {};
			for (my $i=0; $i<scalar(@cols); $i++) {
				my $val = "$cols[$i]";
				my $key = "$colheaders[$i]";
				$watch{"$table"}->{"$uid"}->{"$key"} = "$val";
			}
		}
	}
	close $fh;
}


# Now add the update data to the tables in %watch.
#
opendir(my $updFH, "$updir") or die "Although I found it, I am unable to open $updir: $!";
my @upfiles  = grep { (/^update\..+\.txt$/) && (!/^~/) && -f "$updir/$_" } readdir($updFH);
for my $f (@upfiles) {
	#print "$_\n";
	my $table = "$f";
	$table =~ s/update\.(.+?)\.txt/$1/gi;
	next unless defined $table and defined $orderedcols{"$table"};

	$watch{"$table"} = {} unless defined $watch{"$table"};
	#print "$table\n";
	my @colheaders = ();
	my $linenum = 1;
	my $idcol;

	open(my $fh, "<", "$updir/$f") or die "FATAL: Unable to open $updir/$f for reading: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		my @cols = split /\t/, "$line", -1;
		if ($linenum == 1) {
			$linenum++;
			for (my $i=0; $i<scalar(@cols); $i++) {
				my $val = $cols[$i];
				push @colheaders, "$val";
				$idcol = $i if "$val" eq "$orderedcols{$table}->[0]";
			}
		} else {
			my $uid = "$cols[$idcol]";
			# The validation step should have caught ids in any table duplicated between the update 
			# and the existing data, except for the sample table. 
			#
			# A fraction of samples, including but not just the newest, are stored in AT so there 
			# will be some overlap with the watchDB. We don't need to load these dups into %watch 
			# but they should not throw an error.
			#
			if (defined $watch{"$table"}->{"$uid"}) {
				unless ("$table" eq "sample") {
					print "INFO : $table id $uid already exists in either the update or the main watchdb (before update).\n";
					print "INFO : $uid will not be added again!\n";
					print "INFO : It is unclear how this passed the validation stage, so you may want to look into it.\n";
				}
				next;
			}
			$watch{"$table"}->{"$uid"} = {};
			for (my $i=0; $i<scalar(@cols); $i++) {
				my $val = "$cols[$i]";
				my $key = "$colheaders[$i]";
				$watch{"$table"}->{"$uid"}->{"$key"} = "$val";
			}
		}
	}
	close $fh;
}

#print Dumper(\%watch);
#die;

foreach my $table (keys %watch) {
	open(my $fh, ">", "$dbdir/watchdb.$table.txt.tmp") or die "Unable to open $dbdir/watchdb.$table.txt.tmp for writing: $!\n";
	print $fh join("\t", @{$orderedcols{"$table"}}) . "\n";
	foreach my $uid (sort keys %{$watch{"$table"}}) {
		print $fh "$uid";
		for (my $i=1; $i<scalar(@{$orderedcols{"$table"}}); $i++) {
			my $colheader = $orderedcols{"$table"}->[$i];
			my $val = "";
			$val = $watch{"$table"}->{"$uid"}->{"$colheader"} if defined $watch{"$table"}->{"$uid"}->{"$colheader"};
			print $fh "\t$val";
		}
		print $fh "\n";
	}
}

foreach my $table (keys %watch) {
	`mv $dbdir/watchdb.$table.txt.tmp $dbdir/watchdb.$table.txt`;
}


exit 0;


