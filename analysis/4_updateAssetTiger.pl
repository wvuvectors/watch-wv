#! /usr/bin/env perl

use strict;
use warnings;
use DateTime qw( );

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options]\n";
$usage   .=   "Build update table for AssetTiger database.\n";
$usage   .=   "\n";

my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');

my $WATCHFILE      = "updates/watchdb.LATEST.txt";
my $FIELDSFILE     = "resources/fields_atgr.txt";
my $ATFILE_UPDATE  = "updates/asset_tiger.$NOW.txt";


# keyed master hash for all compiled input data
my %asset2data = ();

#
# Input from latest WaTCH database.
# Use this to populate the keys of the master hash %asset2data
#
my @keysW = ();
my $count = 0;

if (-f "$WATCHFILE") {
	open (my $watchInFH, "<", "$WATCHFILE") or die "Unable to open WATCHFILE for reading: $!\n";
	while (my $line = <$watchInFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @fields = split "\t", "$line", -1;

	if ($count == 0) {
		# populate the keys array with values on the first line
		for (my $j=0; $j < scalar(@fields); $j++) {
			$keysW[$j] = $fields[$j];
		}
	} else {
			# first field is the asset ID
			my $asset_id = $fields[0];
			# populate the hash entry for this asset with the remaining values
			$asset2data{$asset_id} = {} unless defined $asset2data{$asset_id};
			for (my $j=1; $j < scalar(@fields); $j++) {
				$fields[$j] = "" if "$fields[$j]" eq "NaN" or "$fields[$j]" eq "-" or "$fields[$j]" eq "none";
				$asset2data{"$asset_id"}->{"$keysW[$j]"} = $fields[$j];
				if ("Description" eq "$keysW[$j]" and "" eq "$fields[$j]") {
					$asset2data{"$asset_id"}->{"$keysW[$j]"} = "Sample";
				}
			}
		}
		$count++;
	}
	close $watchInFH;
}


# read in the fields for the AssetTiger database
my @atdb_fields = ();
open (my $atFieldsFH, "<", "$FIELDSFILE") or die "Unable to open $FIELDSFILE for reading: $!";
while (<$atFieldsFH>) {
	chomp;
	push @atdb_fields, "$_";
}
close $atFieldsFH;


#
# Write the AssetTiger database update to an Excel file
# AT require Excel for data import :-@
#
open (my $atOutFH, ">", "$ATFILE_UPDATE") or die "Unable to open $ATFILE_UPDATE for writing: $!";
print $atOutFH join("\t", @atdb_fields) . "\n";
foreach my $sample_id (keys %asset2data) {
	print $atOutFH "$sample_id";
	for (my $i=1; $i < scalar(@atdb_fields); $i++) {
		my $field = $atdb_fields[$i];
		my $value = "";
		$value = trim($asset2data{$sample_id}->{"$field"}) if defined $asset2data{$sample_id}->{"$field"};
		print $atOutFH "\t$value";
	}
	print $atOutFH "\n";
}
close $atOutFH;


exit 0;



sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


