#! /usr/bin/env perl

use strict;
use warnings;

use POSIX 'ceil';

use Text::CSV qw(csv);

use List::Util qw(sum);

use Date::Format;
use DateTime qw( );
use DateTime::Duration;
use DateTime::Format::Strptime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use Time::Piece;
use Time::Local;


use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname\n";
$usage   .=  "Pre-calculates alert levels for faster dashboard performance and writes them ";
$usage   .=  "in tab-delimited format to STDOUT.\n";
$usage   .=   "\n";

my $status = 0;

while (@ARGV) {
  my $arg = shift;
  if ("$arg" eq "-h") {
		die "$usage";
	}
}

my $ABUND_ALERT_WINDOW = 12;	# Number of weeks to use in determining abundance alert level.
my $TREND_ALERT_WINDOW = 4;		# Number of consecutive samples to use in trend calculations.
my $SPIKE_THRESHOLD 	 = 5;		# Trigger (X) for identifying spikes.

my %targ2disease = (
	"Influenza Virus A (FluA)" => "FLUA", 
	"Influenza Virus B (FluB)" => "FLUB", 
	"Human Norovirus GII (HuNoV-GII)" => "NoV", 
	"SARS-CoV-2" => "COVID", 
	"Respiratory Syncitial Virus, Human (RSV)" => "RSV"
);

my %targ2loci = (
	"Influenza Virus A (FluA)" => {"M" => 1}, 
	"Influenza Virus B (FluB)" => {"NEP/NS1" => 1}, 
	"Human Norovirus GII (HuNoV-GII)" => {"ORF1_ORF2" => 1}, 
	"SARS-CoV-2" => {"SC2" => 1, "N2" => 1}, 
	"Respiratory Syncitial Virus, Human (RSV)" => {"G" => 1}
);

my %infiles = ("../dashboard/data/watchdb.result.txt" => 1,
							 "../dashboard/data/mu.result.txt" 			=> 1);


# Read resource tables into hash
my %resources = ();
my $resource_wkbk = ReadData("resources/watchdb.all_tables.xlsx", dtfmt => "mm/dd/yy");
#print Dumper($resource_wkbk);
#die;

foreach my $sheet_name (keys %{$resource_wkbk->[0]->{"sheet"}}) {
#	print "$sheet_name\n";
	next unless "$sheet_name" eq "location";
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

my %do_not_rollup = ();

my %loc2county = ();
foreach my $locid (keys %{$resources{"location"}}) {
	next unless $resources{"location"}->{"$locid"}->{"location_status"} eq "active";
	my $county = $resources{"location"}->{"$locid"}->{"location_counties_served"};
	$loc2county{"$locid"} = "$county";
	$do_not_rollup{"$locid"} = 1 if $resources{"location"}->{"$locid"}->{"location_category"} ne "wwtp";
}

#print Dumper(\%loc2county);
#print Dumper(\%resources);
#die;

# Used to pre-fill the dated_abundance hash.
#
my @epiweek_bounds = (2050.1, 2019.11);

# Read result tables into hash, keyed by county, location_id, target, and sample id
my %results = ();

foreach my $f (keys %infiles) {
	my @colnames = ();
	my $linenum  = 0;

	my $loc_id_col = -1;
	my $uid_col = -1;
	my $targ_col = -1;
	my $locus_col = -1;
	my $dt_col = -1;
	my $val1_col = -1;
	my $val2_col = -1;
	my $basis_col = -1;
	my $sample_valid_col = -1;
	my $result_valid_col = -1;
	my $extreme_col = -1;

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
				$uid_col = $i if "$cols[$i]" eq "assay_id";
				$loc_id_col = $i if "$cols[$i]" eq "location_id";
				$targ_col = $i if "$cols[$i]" eq "target";
				$locus_col = $i if "$cols[$i]" eq "target_genetic_locus";
				$sample_valid_col = $i if "$cols[$i]" eq "sample_qc";
				$result_valid_col = $i if "$cols[$i]" eq "target_result_validated";
				$dt_col = $i if "$cols[$i]" eq "collection_end_datetime";
				$val1_col = $i if "$cols[$i]" eq "target_copies_fn_per_cap";
				$val2_col = $i if "$cols[$i]" eq "target_copies_per_ld";
				$basis_col = $i if "$cols[$i]" eq "target_per_capita_basis";
				$extreme_col = $i if "$cols[$i]" eq "is.extreme";
			}
		} else {
			# Extract the data for this row, keyed by the ID
			my ($thisId, $thisLoc, $thisTarg, $thisLocus, $thisDateTime) = ("$cols[$uid_col]", "$cols[$loc_id_col]", "$cols[$targ_col]", "$cols[$locus_col]", "$cols[$dt_col]");
			my $val = $cols["$val1_col"];
			if ("$val" eq "NA") {
				$val = $cols["$val2_col"];
			} else {
				$val = $val / $cols["$basis_col"];
			}
			
			# ignore results that are not validated (sample fail, NTC too high, etc.)
			next unless lc("$cols[$sample_valid_col]") eq "pass";
			next if lc("$cols[$result_valid_col]") eq "ntc above threshold";
			next unless defined $targ2disease{"$thisTarg"};
			next unless defined $loc2county{"$thisLoc"};
			
			# skip any extreme values
			next if lc("$cols[$extreme_col]") eq lc("TRUE");

			# skip any results that don't match one of the accepted gene loci
			next unless defined $targ2loci{"$thisTarg"}->{"$thisLocus"};
			
			# skip any NA abundance values
			next if lc("$cols[$val1_col]") eq lc("NA");
			next if lc("$cols[$val2_col]") eq lc("NA");

			# date format hack area to make sure all times are in 24-hour format.
			# Also calculate a collection week as year.week (eg, 2023.1) to enable roll-up to 
			# county and state.
			$thisDateTime =~ s/ AM$//i;
			my ($m, $d, $y, $h, $min) = ("", "", "", "", "");
			if ("$thisDateTime" =~ /(\d+)\/(\d+)\/(\d+) (\d+):(\d+) PM$/) {
				($m, $d, $y, $h, $min) = ($1, $2, $3, $4, $5);
				$h += 12;
				$thisDateTime = "$m/$d/$y $h:$min";
			} else {
				my ($date, $time) = split " ", "$thisDateTime";
				($m, $d, $y) = split "/", "$date";
			}					
			$y = "20$y" unless ($y =~ /^20/);
			my $epoch = timelocal( 0, 0, 0, $d, $m - 1, $y - 1900 );
			my $week  = strftime( "%U", @{localtime($epoch)} );
			$week = "0$week" unless ($week =~ /^\d{2}$/);
			my $epi_week = "$y.$week";
			
			# Also update the oldest and newest epi weeks.
			#
			$epiweek_bounds[0] = $epi_week unless $epiweek_bounds[0] < $epi_week;
			$epiweek_bounds[1] = $epi_week unless $epiweek_bounds[1] > $epi_week;
			
			# Add this datum to the results hash, keyed on epi week.
			# 
			my $thisCounty = $loc2county{"$thisLoc"};
			$results{$epi_week} = {} unless defined $results{$epi_week};
			$results{$epi_week}->{"$thisCounty"} = {} unless defined $results{$epi_week}->{"$thisCounty"};
			$results{$epi_week}->{"$thisCounty"}->{"$thisLoc"} = {} unless defined $results{$epi_week}->{"$thisCounty"}->{"$thisLoc"};
			$results{$epi_week}->{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"} = [] unless defined $results{$epi_week}->{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"};
			push @{$results{$epi_week}->{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}}, $val;
			#push @{$results{$epi_week}->{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}}, "$thisId";

		}
		$linenum++;
	}
	close $fh;
}
#print Dumper(\%results);
#print Dumper(\@epiweek_bounds);
#die;

my %dated_abundances = ();
 
foreach my $epiweek (sort {$b <=> $a} keys %results) {
	$dated_abundances{"WV"} = {} unless defined $dated_abundances{"WV"};
	# Set up array to hold the values for all locations in the state
	my %state_vals = ();
	foreach my $county (keys %{$results{$epiweek}}) {
		$dated_abundances{"$county"} = {} unless defined $dated_abundances{"$county"};
		# Set up array to hold the values for all locations in the county
		my %county_vals = ();
		foreach my $loc (keys %{$results{$epiweek}->{"$county"}}) {
			$dated_abundances{"$loc"} = {} unless defined $dated_abundances{"$loc"};
			foreach my $targ (keys %{$results{$epiweek}->{"$county"}->{"$loc"}}) {
				# Add this target to the state and county rollup arrays unless it already exists
				$state_vals{"$targ"} = [] unless defined $state_vals{"$targ"};
				$county_vals{"$targ"} = [] unless defined $county_vals{"$targ"};
				
				# Add the target to the dated_abundances hashes
				$dated_abundances{"WV"}->{"$targ"} = {} unless defined $dated_abundances{"WV"}->{"$targ"};
				$dated_abundances{"$county"}->{"$targ"} = {} unless defined $dated_abundances{"$county"}->{"$targ"};
				$dated_abundances{"$loc"}->{"$targ"} = {} unless defined $dated_abundances{"$loc"}->{"$targ"};

				# In case there are multiple samples for this epi week, make sure to take the average
				my $mean = calcSimpleMean($results{$epiweek}->{"$county"}->{"$loc"}->{"$targ"});
				# Store the average (or value) in the dated_abundances hash for this location, target, and epi week
				$dated_abundances{"$loc"}->{"$targ"}->{$epiweek} = $mean;
				# Temporarily put the average (or value) in the arrays set up to hold county and state rollups
				push @{$county_vals{"$targ"}}, $mean unless defined $do_not_rollup{"$loc"};
				push @{$state_vals{"$targ"}}, $mean unless defined $do_not_rollup{"$loc"};
			}
		}
		# county rollup
		foreach my $targ (keys %county_vals) {
			my $county_mean = calcSimpleMean($county_vals{"$targ"});
			$dated_abundances{"$county"}->{"$targ"}->{$epiweek} = $county_mean;
		}
	}
	# state rollup
	foreach my $targ (keys %state_vals) {
		my $state_mean = calcSimpleMean($state_vals{"$targ"});
		$dated_abundances{"WV"}->{"$targ"}->{$epiweek} = $state_mean;
	}
}
#print Dumper(\%dated_abundances);
#die;


print "region_name\tregion_geolevel\ttarget\tepi_year\tepi_week\tabundance_change\tabundance_level\ttrend_slope\ttrend_level\n";
foreach my $region (sort keys %dated_abundances) {

	my $geolevel = "county";
	if (defined $loc2county{"$region"}) {
		$geolevel = "facility";
	} elsif ("$region" eq "WV") {
		$geolevel = "state";
	}
	
	foreach my $targ (sort keys %{$dated_abundances{"$region"}}) {
		# Traverse the epi weeks from oldest to newest
		my @sorted_eweeks = sort {$a <=> $b} keys %{$dated_abundances{"$region"}->{"$targ"}};
		#print Dumper(\@sorted_eweeks);
		foreach (my $i=0; $i<scalar(@sorted_eweeks); $i++) {
			my $epiweek = $sorted_eweeks[$i];
			my ($y, $w) = split /\./, $epiweek, -1;
			my $prnt_str = "$region\t$geolevel\t$targ\t$y\t$w\t";
			
			#
			# Calculate the abundance level. This is written as an integer so the exact text 
			# can be tweaked in the dashboard. In the simplest case it refers to an assay index.
			#
			
			# Get the abundance for the most recent available epi week.
			my $this_val = $dated_abundances{"$region"}->{"$targ"}->{$epiweek};
			# Extract array of values for determining abundance level of the current value.
			my $al_compare_arr = getComparisonTimeWindow($epiweek, $dated_abundances{"$region"}->{"$targ"}, $ABUND_ALERT_WINDOW);
			# Calculate the abundance level text to use for this epiweek
			my $alevel = getAbundanceLevel($this_val, $al_compare_arr);
			
			$prnt_str .= "$alevel->{change}\t$alevel->{level}\t";

			#
			# Calculate the trend level. This is written as an integer so the exact text 
			# can be tweaked in the dashboard. In the simplest case it refers to an assay index.
			#

			my $tl_compare_arr = getComparisonSampleWindow($epiweek, $dated_abundances{"$region"}->{"$targ"}, $TREND_ALERT_WINDOW);
			# Calculate the trend level text to use for this epiweek
			my $tlevel = getTrendLevel($this_val, $tl_compare_arr);
			
			$prnt_str .= "$tlevel->{slope}\t$tlevel->{level}";
			
			print "$prnt_str\n";
		}
	}
}

exit $status;



sub getComparisonTimeWindow {
	my $this_wk    = shift;
	my $value_hash = shift;
	my $win_width  = shift;
	
	# Working backward in time from $this_wk, not including $this_wk, retrieve the appropriate 
	# values from the keyed values in $value_hash. Since this is time-based, not sample-based, 
	# we can simply step through the hash $win_width times.
	# Store the values in a new array @comp.
	#
	my @comp = ();
	my ($y, $w) = split /\./, $this_wk;
	# Step backwards through the epi weeks until we reach $win_width time points in the past.
	for (my $j=0; $j < $win_width; $j++) {
		$w = $w - 1;
		if ($w < 0 or ($w == 0 and !defined $value_hash->{"$y.$w"})) {
			$y--;
			$w = 52;
		}
		push(@comp, $value_hash->{"$y.$w"}) if defined $value_hash->{"$y.$w"};
	}
	#print Dumper(\@comp);
	return \@comp;
}


sub getComparisonSampleWindow {
	my $this_wk = shift;
	my $value_hash = shift;
	my $win_width = shift;
	
	# Working backward in time from $this_wk, not including $this_wk, retrieve the appropriate 
	# values from the keyed values in $values_hash. Since this is sample-based, we need to track 
	# how many samples we actually find.
	# Store the values in a new array @comp.
	#
	my @comp = ();
	my ($y, $w) = split /\./, $this_wk;
	
	my $stopper = $win_width * 2;
	# Step backwards through the epi weeks until we extract $win_width samples.
	while (scalar(@comp) <= $win_width and $stopper > 0) {
		$w = $w - 1;
		if ($w < 0 or ($w == 0 and !defined $value_hash->{"$y.$w"})) {
			$y--;
			$w = 52;
		}
		push(@comp, $value_hash->{"$y.$w"}) if defined $value_hash->{"$y.$w"};
		$stopper--;
	}
	#print Dumper(\@comp);
	return \@comp;
}


sub getAbundanceLevel {
	my $val = shift;
	my $arr = shift;
	
	my %result = ("change" => "NA", "level" => 1);
	
	my $mean = calcSimpleMean($arr);
	# Handle edge case: no values in the array returns NA.
	return \%result if "$mean" eq "NA";	# Not enough data
	
	if ($mean > 0) {
		$result{"change"} = ($val-$mean)/$mean;
		if ($result{"change"} <= -0.5) {
			$result{"level"} = 3;	# LOW
		} elsif ($result{"change"} <= 0) {
			$result{"level"} = 4;	# MODERATE
		} elsif ($result{"change"} <= 0.5) {
			$result{"level"} = 5;	#HIGH
		} elsif ($result{"change"} > 1.5) {
			$result{"level"} = 6;	#VERY HIGH
		}
	} else {
		# If the mean is actually 0, we need a special case.
		$result{"change"} = "NA";	
		$result{"level"} = 2;	
	}

	return \%result;
}



sub getTrendLevel {
	my $val = shift;
	my $arr = shift;
	
	my %result = ("slope" => "NA", "level" => 1);
	return \%result if scalar @$arr < 2;	# Not enough data.
			
	# Add the current value $val to the front of the comparison array $arr.
	# This makes the length of $arr equal to $TREND_LEVEL_WINDOW.
	unshift(@$arr, $val);
	
	# Convert the input data into x,y pairs (points) where x is the date index and y is 
	# the abundance value. This will allow us to calculate the slope of a linear trend line 
	# through the points and determine the direction and magnitude (slope) of the change.
	#
	my @points = ();
	for (my $i=0; $i<scalar(@$arr); $i++) {
		my $x = $i+1;
		my $y = $arr->[$i];
		push(@points, [$x, $y]);
	}
		
	# Look for a spiking trend, where the slope of the line connecting the most recent pair 
	# of points exceeds the SPIKE_THRESHOLD. If the actual values are less than 1, this step 
	# is skipped.
	#
 	if ($points[1]->[0] > 1) {
 		my @latest2 = ($points[0], $points[1]);
 		my $sline = calcLSM(\@latest2);
		if ($sline->{"slope"} > 8) {
			$result{"slope"} = $sline->{"slope"};
			$result{"level"} = 2;
		}
 	}
 	
 	unless ($result{"level"} == 2) {

		my $line = calcLSM(\@points);
		my $slope = $line->{"slope"};
		$result{"slope"} = $slope;
		
		if ($slope <= -5) {
			$result{"level"} = 3; # RAPIDLY DECREASING
		} elsif ($slope <= -0.5) {
			$result{"level"} = 4; # DECREASING
		} elsif ($slope <= 0.5) {
			$result{"level"} = 5; # STABLE
		} elsif ($slope <= 5) {
			$result{"level"} = 6; # INCREASING
		} elsif ($slope > 5) {
			$result{"level"} = 7; # RAPIDLY INCREASING
		}
	}
	
	return \%result;
}


sub calcSimpleMean {
	my $a = shift;
	
	return "NA" if scalar @{$a} == 0;
	
	my $sum = 0;
	foreach my $val (@{$a}) {
		$sum += $val;
	}
	my $mean = $sum / scalar @{$a};

	return $mean;
}


sub calcLSM {
	
	my $points = shift; # array of arrays of paired points.
	
	my %lsm = ("slope" => "NA", "intercept" => "NA");
	
	my $n 		= scalar @$points;
	my $sumX  = 0;
	my $sumY  = 0;
	my $sumX2 = 0;
	my $sumXY = 0;
	
	
	foreach my $point (@$points) {
		my ($x, $y) = @$point;
		$sumX  += $x;
		$sumY  += $y;
		$sumX2 += $x * $x;
		$sumXY += $x * $y;
	}
	
	# Calculate the slope (m)
	my $m = ($n * $sumXY - $sumX * $sumY) / ($n * $sumX2 - $sumX * $sumX);
	
	# Calculate the y-intercept (b)
	my $b = ($sumY - $m * $sumX) / $n;
	
	$lsm{"slope"} = $m;
	$lsm{"intercept"} = $b;
	
	return \%lsm;
}



