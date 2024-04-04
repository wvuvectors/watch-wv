#! /usr/bin/env perl

use strict;
use warnings;

use Date::Format;
use DateTime qw( );
#use DateTime::Format::Strptime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname OUTDIR\n";
$usage   .=  "Build resource tables for WVWD using resource Excel file.\n";
$usage   .=   "\n";


my $status = 0;
my $resdir  = "resources";
my $outdir;

while (@ARGV) {
  my $arg = shift;
  if ("$arg" eq "-h") {
		die "$usage";
	} else {
		$outdir = "$arg";
	}
}

die "An output directory is required, but none was provided.\n$usage\n" unless defined $outdir;
die "The provided directory $outdir is not a directory.\n$usage\n" unless -d "$outdir";


#
# Get time stamp for now
#
my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');


my %resources = (
	"location" => {},
	"wwtp"		 => {},
	"county"	 => {},
	"lab"			 => {}
);


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

	my $outf = "$outdir/watchdb_${sheet_name}.txt";
	open(my $fh, ">", "$outf") or die "Unable to open $outf for writing: $!";
	print $fh join("\t", @field_names) . "\n";
	foreach my $uid (keys %{$resources{"$sheet_name"}}) {
		my $count = 0;
		foreach my $field (@field_names) {
			print $fh "\t" unless $count == 0;
			print $fh "$resources{$sheet_name}->{$uid}->{$field}";
			$count = 1;
		}
		print $fh "\n";
	}
	close $fh;

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


sub formatDT {
	# Format the date and time appropriately. I hate date hacking.
	# DateTime::Format::Strptime doesn't work although I have not given up hope.
	#
	my $dtstring = shift;
	
	my %dtHash = ("date" => "", "time" => "");
	my ($mon, $day, $yr, $hrs, $min);
	
	$dtstring =~ s/ AM//gi;
	if (($mon, $day, $yr, $hrs, $min) = $dtstring =~ /(.+?)\/(.+?)\/(.+?) (.+?):(.+?)/) {
		$yr = "20$yr" if length $yr == 2;
		$mon = "0$mon" if length $mon == 1;
		$day = "0$day" if length $day == 1;
		$dtHash{"date"} = "$yr-$mon-$day";

		$hrs = "0$hrs" if length $hrs == 1;
		$min = "0$min" if length $min == 1;
		$dtHash{"time"} = "$hrs:$min";

	} elsif (($mon, $day, $yr) = $dtstring =~ /(.+?)\/(.+?)\/(.+)/) {
		$yr = "20$yr" if length $yr == 2;
		$mon = "0$mon" if length $mon == 1;
		$day = "0$day" if length $day == 1;
		$dtHash{"date"} = "$yr-$mon-$day";
	}
	
	return \%dtHash;
}



