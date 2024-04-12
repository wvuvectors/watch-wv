#! /usr/bin/env perl

use strict;
use warnings;

use Text::CSV qw(csv);

use List::Util qw(sum);

use Date::Format;
use DateTime qw( );
use DateTime::Duration;
use DateTime::Format::Strptime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname\n";
$usage   .=  "Pre-calculates delta and freshness for faster dashboard performance.\n";
$usage   .=   "\n";


my $status = 0;

while (@ARGV) {
  my $arg = shift;
  if ("$arg" eq "-h") {
		die "$usage";
	}
}

#die "A file path must be provided to $progname!\n\n$usage\n" unless defined $muf;
#die "$muf is not a readable file!\n\n$usage\n" unless -f "$muf";

my $dur_mo = DateTime::Duration->new(
	months       => 1,
	end_of_month => 'wrap'
);

my $dur_yr = DateTime::Duration->new(
	years				 => 1,
	end_of_month => 'wrap'
);

my $today = DateTime->today(time_zone => 'local');
$today->set_time_zone('UTC');

# write a tab-delim file deltas.txt
# region_name
# region_geolevel (county or facility)
# region_lab (zoowvu, muidsl, or other)
# annual mean value for the region, for each disease
# last month's mean value for the region, for each disease
# delta from annual mean, for each disease
# data freshness (days since last update), for each disease

my $pdate = "collection_end_datetime";

my %targ2disease = (
	"Influenza Virus A (FluA)" => "FLUA", 
	"Influenza Virus B (FluB)" => "FLUB", 
	"Human Norovirus GII (HuNoV-GII)" => "NoV", 
	"SARS-CoV-2" => "COVID", 
	"Respiratory Syncitial Virus, Human (RSV)" => "RSV"
);

my %infiles = ("../dashboard/data/watchdb.result.txt" => "result",
							 "../dashboard/data/mu.result.txt" 			=> "result",
							 "../dashboard/data/watchdb.sample.txt" => "sample", 
							 "../dashboard/data/mu.sample.txt" 			=> "sample"); 

my %tables  = ("result" => {}, 
							 "sample" => {});


# Read table data into tables hash
foreach my $f (keys %infiles) {
	my $tableRef = $tables{"$infiles{$f}"};
	my @colnames = ();
	my $linenum  = 0;
	open (my $fh, "<", "$f") or die "Unable to open $f for reading: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		next if "$line" =~ /^\s*$/;
		my @cols = split "\t", "$line", -1;
		if ($linenum == 0) {
			# First line of the file contains the column names
			# Each row is a hash keyed by the uid
			foreach (my $i=0; $i<scalar(@cols); $i++) {
				push @colnames, "$cols[$i]";
			}
		} else {
			# Extract the data for this row, keyed by the ID
			my $thisId = "$cols[0]";
			$tableRef->{"$thisId"} = {};
			for (my $i=0; $i<scalar(@cols); $i++) {
				if (!defined $colnames[$i]) {
					print "DEBUG:: $thisId col name ($i) does not exist!\n";
					print "$line\n";
					print Dumper(@colnames);
					exit -1;
				}
				$tableRef->{"$thisId"}->{"$colnames[$i]"} = "$cols[$i]";
			}
		}
		$linenum++;
	}
	close $fh;
}
#print Dumper(\%tables);
#die;


# WaTCH resource tables read into their own hash
my %resources = ();

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
			next unless defined $sheetRef->{"cell"}->[1]->[$j] and "$sheetRef->{cell}->[1]->[$j]" ne "";
			my $key = "$sheetRef->{cell}->[1]->[$j]";
			my $value = "";
			$value = "$colRef->[$j]" if defined $colRef->[$j];
			$resources{"$sheet_name"}->{"$key"}->{"$field"} = "$value";
		}
	}
}
#print Dumper(\%resources);
#die;


my %printable = ();

# Output file columns, metadata:
# region_name	region_geolevel	region_lab
foreach my $locid (keys %{$resources{"location"}}) {
	my $loc_name = $resources{"location"}->{"$locid"}->{"location_common_name"};
	my $geo = "facility";
	$geo = "upstream" unless "$resources{location}->{$locid}->{location_category}" eq "wwtp";
	my $county = $resources{"location"}->{"$locid"}->{"location_counties_served"};
	my $lab = $resources{"location"}->{"$locid"}->{"location_primary_lab"};
	$printable{"$locid"} = {
		"region" 					=> "$loc_name", 
		"region_geolevel" => "$geo",
		"region_lab"			=> "$lab"
	};
	foreach my $target (keys %targ2disease) {
		$printable{"$locid"}->{"$target"} = {
				"vals_yr" 				=> [],
				"vals_mo" 				=> [],
				"included_dates"  => []
		};
	}
	$printable{"$county"} = {
		"region" 					=> "$county", 
		"region_geolevel" => "county",
		"region_lab"			=> "$lab"
	};
	foreach my $target (keys %targ2disease) {
		$printable{"$county"}->{"$target"} = {
				"vals_yr" 				=> [],
				"vals_mo" 				=> [],
				"included_dates"  => []
		};
	}
}

# Output file columns, per disease:
# disease1_fresh	disease1_delta	disease1_annual_mean	disease1_lastmo_mean	...
foreach my $uid (keys %{$tables{"result"}}) {
	my $resultRef = $tables{"result"}->{"$uid"};
	my $locid     = $resultRef->{"location_id"};
	my $target    = $resultRef->{"target"};
	my $dtObj     = makeDT($resultRef->{"$pdate"});
	my $key       = "target_copies_fn_per_cap";
	$key = "target_copies_per_ld" if isEmpty($resultRef->{"sample_flow"});
	
	$dtObj->set_time_zone('UTC');
	my $date_diff = $today - $dtObj;
	#print Dumper($dtObj->ymd());
	#print Dumper($diff_mo);
	#print Dumper($diff_mo->delta_months());
	#die;
	if ($date_diff->delta_months() < 12 or ($date_diff->delta_months() == 12 and $date_diff->delta_days() == 0)) {
		push @{$printable{"$locid"}->{"$target"}->{"vals_yr"}}, $resultRef->{"$key"};
		#push @{$printable{"$locid"}->{"$target"}->{"included_dates"}}, $dtObj;
	}
	if ($date_diff->delta_months() < 1 or ($date_diff->delta_months() == 1 and $date_diff->delta_days() == 0)) {
		push @{$printable{"$locid"}->{"$target"}->{"vals_mo"}}, $resultRef->{"$key"};
	}
}
#print Dumper(\%printable);
#die;

my @targets_sorted = sort keys %targ2disease;
print "region_name\tregion_geolevel\tregion_lab";
foreach my $t (@targets_sorted) {
	my $d = $targ2disease{"$t"};
	print "\tmonthly_mean_$d\tannual_mean_$d\tlast_val_$d\tmonthly_delta_$d\tannual_delta_$d";
}
print "\n";

foreach my $rkey (keys %printable) {
	my $regRef = $printable{"$rkey"};
	my $name = $regRef->{"region"};
	my $geo  = $regRef->{"region_geolevel"};
	my $lab  = $regRef->{"region_lab"};
	
	print "$name\t$geo\t$lab";
	# Calc mean of all val_yr & val_mo entries for each target a given region, and store both the values and the % changes
	foreach my $t (@targets_sorted) {
		my $d = $targ2disease{"$t"};
		my $m = "NA";
		unless (scalar(@{$regRef->{"$t"}->{"vals_mo"}}) == 0) {
			$m = sum(@{$regRef->{"$t"}->{"vals_mo"}}) / scalar(@{$regRef->{"$t"}->{"vals_mo"}});
		}
		my $a = "NA";
		unless (scalar(@{$regRef->{"$t"}->{"vals_yr"}}) == 0) {
			$a = sum(@{$regRef->{"$t"}->{"vals_yr"}}) / scalar(@{$regRef->{"$t"}->{"vals_yr"}});
		}
		print "\t$m\t$a\tNA\tNA\tNA";
	}
	print "\n";
}


exit $status;



sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


sub isEmpty {
	my $val = shift;
	
	my $is_empty = 0;
	$is_empty = 1 if "$val" eq "NaN" or "$val" eq "-" or "$val" eq "none" or "$val" eq "" or "$val" eq "TBD" or "$val" eq "NA";
	
	return $is_empty;
}


sub makeDT {
	my $dtstr = shift;
	
	my $dformat1 = DateTime::Format::Strptime->new(
		pattern   => '%Y-%m-%d %H:%M',
		time_zone => 'local',
		on_error  => 'croak',
	);
	my $dformat2 = DateTime::Format::Strptime->new(
		pattern   => '%m/%d/%Y %H:%M',
		time_zone => 'local',
		on_error  => 'croak',
	);
	
	#print "$dtstr\n";
	my $dtobj;
	if ("$dtstr" =~ /-/) {
		$dtobj = $dformat1->parse_datetime("$dtstr");
	} elsif ("$dtstr" =~ /\//) {
		$dtobj = $dformat2->parse_datetime("$dtstr");
	} else {
		$dtobj = "";
	}
	
	return $dtobj;
}


