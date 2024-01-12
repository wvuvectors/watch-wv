#! /usr/bin/env perl

use strict;
use warnings;

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use DateTime qw( );

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options] [DAYS]\n";
$usage   .=   "Build update table for AssetTiger database. Pass an integer DAYS to determine ";
$usage   .=   "how many days of samples to include in the update (default = 90). Writes to STDOUT.\n";
$usage   .=   "\n";

my $status = 0;
my $days = 90;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
  } else {
		$days = $arg;
		$days =~ s/\D//gi;
	}
}


#
# Get time stamp for now
#
my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');

my $feed         = "NWSS";
my @feed_headers = ();
my %feed_cvmap   = ();

# Get the required feed fields (column headers), in order.
# These are located in the resources/ folder.
open (my $fFH, "<", "resources/${feed}_fields.txt") or die "FATAL: Unable to open resources/${feed}_fields.txt for reading: $!.";
while (my $line = <$fFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	push @feed_headers, "$line";
}
close $fFH;


# Get the controlled vocabulary map for this feed.
# This is located in the resources/ folder.
my $linenum = 0;
open (my $vFH, "<", "resources/${feed}_cvm.txt") or die "FATAL: Unable to open resources/${feed}_cvm.txt for reading: $!.";
while (my $line = <$vFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	#my @cols = split /\t/, "$line", -1;
	my ($extref_field, $extref_value, $watch_table, $watch_field, $watch_value, $reference) = split /\t/, "$line", -1;
	if ($linenum == 0) {
		$linenum++;
	} else {
		$feed_cvmap{"$extref_field"} = {} unless defined $feed_cvmap{"$extref_field"};
		$feed_cvmap{"$extref_field"}->{"$watch_field"} = {} unless defined $feed_cvmap{"$extref_field"}->{"$watch_field"};
		$feed_cvmap{"$extref_field"}->{"$watch_field"}->{"$watch_value"} = {"value" => "$extref_value", "reference" => "$reference"};
	}
}
close $vFH;


my %tables = (
	"abatch"				=> {},
	"assay" 				=> {},
	"cbatch" 				=> {},
	"concentration" => {},
	"ebatch" 				=> {},
	"extraction" 		=> {},
	"result" 				=> {},
	"sample"				=> {}
);

my %resources = (
	"location" => {},
	"wwtp"		 => {},
	"county"	 => {},
	"lab"			 => {}
);


#
# Read WaTCH database tables from dbdir.
#

foreach my $table (keys %tables) {
	$linenum = 0;
	my @headers = ();
	open (my $tFH, "<", "$dbdir/watchdb.$table.txt") or die "Unable to open $dbdir/watchdb.$table.txt for reading: $!\n";
	while (my $line = <$tFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @cols = split "\t", "$line", -1;

		if ($linenum == 0) {
			$linenum++;
			for (@cols) {
				push @headers, "$_";
			}
		} else {
			# first field is always the uid regardless of table
			my $uid = $cols[0];

			$tables{"$table"}->{"$uid"} = {};
			for (my $i=0; $i < scalar(@cols); $i++) {
				$tables{"$table"}->{"$uid"}->{"$headers[$i]"} = trim("$cols[$i]");
			}
		}
	}
	close $tFH;
}

#print Dumper(\%tables);
#die;


#
# Read WaTCHdb resources from resources/ dir.
#
# Read resource tables into hash
my $resource_wkbk = ReadData("resources/watchdb.all_tables.xlsx", dtfmt => "mm/dd/yy");
#print Dumper($resource_wkbk);
#die;

foreach my $sheet_name (keys %{$resource_wkbk->[0]->{"sheet"}}) {
#	print "$sheet_name\n";
	my @field_names = ();

	$resources{"$sheet_name"} = {} unless defined $resources{"$sheet_name"};
	my $pos = $resource_wkbk->[0]->{"sheet"}->{"$sheet_name"};
	my $sheetRef = $resource_wkbk->[$pos];
	
	# The field labels are in row A (or, index 1 of each array in 'cell' ). Store these in an array.
	for (my $i=1; $i < scalar(@{$sheetRef->{"cell"}}); $i++) {
		#print "$sheetRef->{cell}->[$i]->[1]\n";
		next unless defined $sheetRef->{"cell"}->[$i]->[1] and "$sheetRef->{cell}->[$i]->[1]" ne "";
		my $col_name = "$sheetRef->{cell}->[$i]->[1]";
		push @field_names, "$col_name";
	}
	#print Dumper(\@field_names);
	#last;

	# The keys for each resource sheet are in column A (or, array 1 of 'cell'). 
	# Loop over these to add entries to the resources hash.
	for (my $i=2; $i < scalar(@{$sheetRef->{"cell"}->[1]}); $i++) {
		#print "$sheetRef->{cell}->[1]->[$i]\n";
		next unless defined $sheetRef->{"cell"}->[1]->[$i] and "$sheetRef->{cell}->[1]->[$i]" ne "";
		my $key = "$sheetRef->{cell}->[1]->[$i]";
		$resources{"$sheet_name"}->{"$key"} = {};
	}
	#last;
	
	# To get the actual values, loop over the 'cell' array.
	# Get the key from array 1.
	# Loop over the field_names array to get the values.
	for (my $i=0; $i < scalar(@{$sheetRef->{"cell"}})-1; $i++) {
		next unless defined $sheetRef->{"cell"}->[$i+1] and defined $field_names[$i] and "$field_names[$i]" ne "";
		my $colRef = $sheetRef->{"cell"}->[$i+1];
		my $field = "$field_names[$i]";
		my $resourceRef = $resources{"$sheet_name"};
		# Loop over the values in this column. Add each to the resources hash.
		for (my $j=2; $j < scalar(@{$colRef}); $j++) {
			next unless defined $sheetRef->{"cell"}->[1]->[$j];
			my $key = "$sheetRef->{cell}->[1]->[$j]";
			my $value = "";
			$value = "$colRef->[$j]" if defined $colRef->[$j];
			$resources{"$sheet_name"}->{"$key"}->{"$field"} = "$value";
		}
	}
}

#print Dumper(\%resources);
#die;




exit $status;



sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


