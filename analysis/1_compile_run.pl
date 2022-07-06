#! /usr/bin/env perl

use strict;
use warnings;
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


my %sample2metadata = ();
my %outhead = ();

foreach (@plate_files) {
	my $plate_wkbk = ReadData("$rundir/$_", dtfmt => "mm/dd/yy");
	
	# parse data based on plate type
	my @row = Spreadsheet::Read::row($plate_wkbk->[1], 1);
	parse_plate("$row[1]", $plate_wkbk, \%sample2metadata, \%outhead);
}

#print Dumper(\%sample2metadata);
#die;

# Merge with field data
#
=cut
my %asset2metadata = ();
my @keys = ();
open (my $assetFH, "<", "$rundir/sample_assets.csv") or die "Unable to open $rundir/sample_assets.csv for reading: $!\n";
my $count = 0;
while (my $line = <$assetFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @fields = split /,/, "$line", -1;
	if ($count == 0) {
		for (my $j=0; $j < scalar(@fields); $j++) {
			$asset2metadata{$fields[$j]} = $j;
			$keys[$j] = $fields[$j];
		}
	} else {
		for (my $j=0; $j < scalar(@fields); $j++) {
			$asset2metadata{$keys[$j]} = $fields[$j];
		}
	}
	$count++;
}
close $assetFH;

print Dumper(\%asset2metadata);
die;

# If any samples in the current run are from previous plates, need to get concentration and extraction data from 
# the watchDB file
#
my %sample2watchDB = ();
my @dbkeys = ();
open (my $dbFH, "<", "watchDB.csv") or die "Unable to open watchDB.csv for reading: $!\n";
my $count = 0;
my $sampleIDcol = 0;
while (my $line = <$dbFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @fields = split /,/, "$line", -1;
	if ($count == 0) {
		for (my $j=0; $j < scalar(@fields); $j++) {
			$sample2watchDB{$fields[$j]} = $j;
			$dbkeys[$j] = $fields[$j];
			$sampleIDcol = $j if "$fields[$j]" eq "sampleID";
		}
	} else {
		$sample2watchDB{$fields[$sampleIDcol]} = {};
		for (my $j=0; $j < scalar(@fields); $j++) {
			next if $j == $sampleIDcol;
			$sample2watchDB{$fields[$sampleIDcol]}->{$dbkeys[$j]} = $fields[$j];
		}
	}
	$count++;
}
close $dbFH;
=cut

#print Dumper(\%sample2watchDB);
#die;


open (my $outFH, ">", "$rundir/run_metadata.csv") or die "Unable to open $rundir/run_metadata.csv for writing: $!\n";
print $outFH "\"SampleID\",\"" . join("\"\,\"", sort keys %outhead) . "\"\n";

foreach my $sampleID (keys %sample2metadata) {
	print $outFH "\"$sampleID\"";
	for my $k (sort keys %outhead) {
		my $val = "";
		if (defined $sample2metadata{"$sampleID"}->{"$k"}) {
			$val = $sample2metadata{"$sampleID"}->{"$k"};
		}
		print $outFH ",\"$val\"";
	}
	print $outFH "\n";
}
close $outFH;



sub parse_plate {
	my $src = shift;
	my $wkbk = shift;
	my $s2m = shift;
	my $oh = shift;
	
	my @row2name = ("", "A", "B", "C", "D", "E", "F", "G", "H");
	my %plate_md = ();
	my @mdrows = Spreadsheet::Read::rows($wkbk->[1]);
	for (my $i=1; $i < scalar(@mdrows); $i++) {
		if (index($mdrows[$i][0], $src) == 0) {
			unless ("" eq "$mdrows[$i][0]") {
				my $key = "$mdrows[$i][0]";
				$key =~ s/ (Type)//gi;
				my $val = "$mdrows[$i][1]";
				if ("$key" =~ /Date/gi) {
					my $datetime = DateTime::Format::Excel->parse_datetime($val);
					$val = $datetime->ymd('-');
				}
				$plate_md{"$key"} = "$val";
			}
			$oh->{"$mdrows[$i][0]"} = 1;
		} else {
			unless ("" eq "$mdrows[$i][0]") {
				my $key = "$mdrows[$i][0]";
				$key =~ s/ (Type)//gi;
				my $val = "$mdrows[$i][1]";
				if ("$key" =~ /Date/gi) {
					my $datetime = DateTime::Format::Excel->parse_datetime($val);
					$val = $datetime->ymd('-');
				}
				$plate_md{"$src $key"} = "$val";
			}
			$oh->{"$src $mdrows[$i][0]"} = 1;
		}
	}
	
	$oh ->{"${src} Plate Location"} = 1;
	$oh ->{"${src} Plate Comment"} = 1;
	
	my $row_sheet = 2;
	$row_sheet = 3 if "$src" eq "Assay";
	my $comment_sheet = 3;
	$comment_sheet = 4 if "$src" eq "Assay";
	
	# Fix to ignore entries outside of the identified plate boundaries!
	#
	my @plate_rows = Spreadsheet::Read::rows($wkbk->[$row_sheet]);
	#print Dumper(\@plate_rows);
	#die;
	my @plate_rows_comment = Spreadsheet::Read::rows($wkbk->[$comment_sheet]);
#	print "$src Plate:\n";
	for (my $i=1; $i < scalar(@plate_rows); $i++) {
		#print "i: $i" if $src eq "Assay";
		for (my $j=1; $j < scalar(@{$plate_rows[$i]}); $j++) {
			#print ", j: $j\n" if $src eq "Assay";
			if (!defined $plate_rows[$i][$j]) {
				warn "****WARNING**** $src Plate has no sample ID at row $i, column $j";
				next;
#			} elsif ("" ne "$plate_rows[$i][$j]") {
#				print "  Sample at $i, $j is $plate_rows[$i][$j]\n";
			}
			my $sample_id = "$plate_rows[$i][$j]";
			$sample_id =~ s/^\s+//;
			$sample_id =~ s/\s+$//;
			unless ("" eq "$sample_id" or "buffer" eq lc($sample_id)) {
				#print "  $row2name[$i]" . sprintf("%02d", $j) . "\n" if $src eq "Assay";
				$s2m->{"$sample_id"} = {} unless defined $s2m->{"$sample_id"};
				$s2m->{"$sample_id"}->{"${src} Plate Location"} = "$row2name[$i]" . sprintf("%02d", $j);
				$s2m->{"$sample_id"}->{"${src} Plate Comment"} = "$plate_rows_comment[$i][$j]" if defined $plate_rows_comment[$i][$j];
				foreach my $k (keys %plate_md) {
					$s2m->{"$sample_id"}->{"$k"} = "$plate_md{$k}";
				}
			}
		}
	}
}

=cut
sub sort_samples {
	my ($county_a) = $a =~ /([A-Z]+)/;
	my ($county_b) = $b =~ /([A-Z]+)/;

	my ($site_a) = $a =~ /(\d+)\-/;
	my ($site_b) = $b =~ /(\d+)\-/;

	my ($sample_a) = $a =~ /\-(\d+)/;
	my ($sample_b) = $b =~ /\-(\d+)/;

	# Compare and return
	return $county_a cmp $county_b and $site_a <=> $site_b and $sample_a <=> $sample_b;
}
=cut


exit 0;

