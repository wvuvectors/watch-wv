#! /usr/bin/env perl


use strict;
use warnings;
use DateTime qw( );

use DateTime::Format::Excel;
use Data::Dumper;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options]\n";
$usage   .=   "\n";

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
	}
}

my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');


# keyed master hash for all compiled input data
my %asset2data = ();
my %flu2data   = ();

# Paths for input files
my $WATCHFILE_MAIN = "updates/watchdb.LATEST.txt";
my $FLUFILE_MAIN   = "updates/watchdb_flu.LATEST.txt";

# Paths for output files
my $FLUFILE_INCR   = "updates/fludb/watchdb_flu.$NOW.txt";

#
# Fetch input from latest WaTCH database.
# Use this to populate the keys of the master hash %asset2data
# Skip if there is no database file
#
my @w_field_names = ();
my $count = 0;

if (-f "$WATCHFILE_MAIN") {
	open (my $watchInFH, "<", "$WATCHFILE_MAIN") or die "Unable to open $WATCHFILE_MAIN for reading: $!\n";
	while (my $line = <$watchInFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @values = split "\t", "$line", -1;

	if ($count == 0) {
		# populate the keys array with values on the first line
		for (my $i=0; $i < scalar(@values); $i++) {
			$w_field_names[$i] = $values[$i];
		}
	} else {
			# first field is the asset ID. second is the replicate id.
			my $asset_id = $values[0];
			my $rep_id = $values[1];
			
			# populate the hash entry for this asset with the remaining values
			$asset2data{"$asset_id"} = "NA" unless defined $asset2data{"$asset_id"};
			for (my $i=0; $i < scalar(@values); $i++) {
				$asset2data{"$asset_id"} = "$values[$i]" if "$w_field_names[$i]" eq "Location";
			}
		}
		$count++;
	}
	close $watchInFH;
}

#print Dumper(\%asset2data);
#die;


#
# Fetch input from latest Flu database.
# Skip if there is no database file
#
my @f_field_names = ();
$count = 0;

if (-f "$FLUFILE_MAIN") {
	open (my $watchInFH, "<", "$FLUFILE_MAIN") or die "Unable to open $FLUFILE_MAIN for reading: $!\n";
	while (my $line = <$watchInFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @values = split "\t", "$line", -1;

	if ($count == 0) {
		# populate the keys array with values on the first line
		for (my $i=0; $i < scalar(@values); $i++) {
			$f_field_names[$i] = $values[$i];
		}
	} else {
			# first field is the asset ID. second is the replicate id.
			my $asset_id = $values[0];
			
			# populate the hash entry for this asset with the remaining values
			$flu2data{"$asset_id"} = {};
			for (my $i=0; $i < scalar(@values); $i++) {
				$flu2data{"$asset_id"}->{"$f_field_names[$i]"} = $values[$i];
			}
		}
		$count++;
	}
	close $watchInFH;
}

#print Dumper(\%flu2data);
#die;



my @flu_outfields = ("Sample ID", "Sample Composite Start", "Sample Composite End", "Sample Received Date", "Sample Flow (MGD)", 
										 "Location", "Assay Target 1", "Assay Target 1 Result (CN/L)", "Assay Target 2", "Assay Target 2 Result (CN/L)");

`cp $FLUFILE_MAIN $FLUFILE_MAIN.OLD`;

open (my $watchFH, ">", "$FLUFILE_INCR") or die "Unable to open $FLUFILE_INCR for writing: $!";

print $watchFH join("\t", @flu_outfields) . "\n";
foreach my $asset_id (keys %flu2data) {
	my $fluset = $flu2data{"$asset_id"};
	print $watchFH "$asset_id";
	foreach my $field (@flu_outfields) {
		next if "$field" eq "Sample ID";
		my $value = "NA";
		if ("$field" eq "Location") {
			$value = $asset2data{"$asset_id"};
		} elsif (defined $fluset->{"$field"}) {
			$value = "$fluset->{$field}";
		}
		print $watchFH "\t$value";
	}
	print $watchFH "\n";
}
close $watchFH;

`cp $FLUFILE_INCR $FLUFILE_MAIN`;


exit 0;



