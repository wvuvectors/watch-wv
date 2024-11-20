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

#die "A file path must be provided to $progname!\n\n$usage\n" unless defined $muf;
#die "$muf is not a readable file!\n\n$usage\n" unless -f "$muf";

# my $dur_mo = DateTime::Duration->new(
# 	months       => 1,
# 	end_of_month => 'wrap'
# );
# 
# my $dur_yr = DateTime::Duration->new(
# 	years				 => 1,
# 	end_of_month => 'wrap'
# );
# 
# my $today = DateTime->today(time_zone => 'local');
# $today->set_time_zone('UTC');
# 

my $WINDOW = 3;	# Number of samples to use in trend calculations.
my $SPIKE_THRESHOLD = 500; # Trigger for identifying spikes/despikes; given as percent change.
my $TREND_THRESHOLD = 25; # Trigger for identifying consistent trend; given as percent change.


my %targ2disease = (
	"Influenza Virus A (FluA)" => "FLUA", 
	"Influenza Virus B (FluB)" => "FLUB", 
	"Human Norovirus GII (HuNoV-GII)" => "NoV", 
	"SARS-CoV-2" => "COVID", 
	"Respiratory Syncitial Virus, Human (RSV)" => "RSV"
);

my %infiles = ("../dashboard/data/watchdb.result.txt" => 1,
							 "../dashboard/data/mu.result.txt" 			=> 1);

my %keepers = (
		"sample_id" => 1, 
		"location_id" => 1, 
		"collection_end_datetime" => 1, 
		"target" => 1, 
		"target_copies_per_ld" => 1, 
		"target_copies_fn_per_cap" => 1,
		"target_per_capita_basis" => 1,
		"sample_flow" => 1,  
		"target_result_validated" => 1,
		"sample_qc" => 1
);

my %sentinels = (
	"PrincetonWWTP-01" => 12, 
	"StarCityWWTP-01" => 15, 
	"WheelingWWTP-01" => 12
);

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


my %loc2county = ();
foreach my $locid (keys %{$resources{"location"}}) {
	next unless $resources{"location"}->{"$locid"}->{"location_status"} eq "active";
	my $county = $resources{"location"}->{"$locid"}->{"location_counties_served"};
	$loc2county{"$locid"} = "$county";
}

#print Dumper(\%loc2county);
#print Dumper(\%resources);
#die;


# Read result tables into hash, keyed by county, location_id, target, and sample id
my %results = ();

foreach my $f (keys %infiles) {
	my @colnames = ();
	my $linenum  = 0;
	my $loc_id_col = -1;
	my $uid_col = -1;
	my $targ_col = -1;

	my $sample_valid_col = -1;
	my $result_valid_col = -1;

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
				$sample_valid_col = $i if "$cols[$i]" eq "sample_qc";
				$result_valid_col = $i if "$cols[$i]" eq "target_result_validated";
			}
		} else {
			# Extract the data for this row, keyed by the ID
			my ($thisId, $thisLoc, $thisTarg) = ("$cols[$uid_col]", "$cols[$loc_id_col]", "$cols[$targ_col]");
			# ignore results that are not validated (sample fail, NTC too high, etc.)
			next unless lc("$cols[$sample_valid_col]") eq "pass";
			next unless lc("$cols[$result_valid_col]") eq "yes";
			next unless defined $targ2disease{"$thisTarg"};
			next unless defined $loc2county{"$thisLoc"};
			my $thisCounty = $loc2county{"$thisLoc"};
			$results{"$thisCounty"} = {} unless defined $results{"$thisCounty"};
			$results{"$thisCounty"}->{"$thisLoc"} = {} unless defined $results{"$thisCounty"}->{"$thisLoc"};
			$results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"} = {} unless defined $results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"};
			for (my $i=0; $i<scalar(@cols); $i++) {
				if (!defined $colnames[$i]) {
					print "DEBUG:: $thisId col name ($i) does not exist!\n";
					print "$line\n";
					print Dumper(@colnames);
					exit -1;
				}
				next unless defined $keepers{"$colnames[$i]"};
				$results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}->{"$thisId"} = {} unless defined $results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}->{"$thisId"};

				# date format hack area to make sure all times are in 24-hour format.
				# Also calculate a collection week as year.week (eg, 2023.1) to enable roll-up to 
				# county and state.
				if ("$colnames[$i]" eq "collection_end_datetime") {
				
					$cols[$i] =~ s/ AM$//i;
					my ($m, $d, $y, $h, $min) = ("", "", "", "", "");
					if ("$cols[$i]" =~ /(\d+)\/(\d+)\/(\d+) (\d+):(\d+) PM$/) {
						($m, $d, $y, $h, $min) = ($1, $2, $3, $4, $5);
						$h += 12;
						$cols[$i] = "$m/$d/$y $h:$min";
					} else {
						my ($date, $time) = split " ", "$cols[$i]";
						($m, $d, $y) = split "/", "$date";
					}					
					my $epoch = timelocal( 0, 0, 0, $d, $m - 1, $y - 1900 );
					my $week  = strftime( "%U", @{localtime($epoch)} );
					$y = "20$y" unless ($y =~ /^20/);
					$results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}->{"$thisId"}->{"collection_week"} = "$y.$week";
				}
				$results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}->{"$thisId"}->{"$colnames[$i]"} = "$cols[$i]";
			}
		}
		$linenum++;
	}
	close $fh;
}
#print Dumper(\%results);
#die;


# Build a calculation hash keyed on county, location_id, and target.
# for each location, include its county and the array of most recent abundance values to use in temporal order
# 
my %ordered_uids = ();
foreach my $county (keys %results) {
	$ordered_uids{"$county"} = {};
	foreach my $locid (keys %{$results{"$county"}}) {
		$ordered_uids{"$county"}->{"$locid"} = {};
		foreach my $targ (keys %{$results{"$county"}->{"$locid"}}) {
			$ordered_uids{"$county"}->{"$locid"}->{"$targ"} = [];
			my %result = %{$results{"$county"}->{"$locid"}->{"$targ"}};
			my @a = sort {Time::Piece->strptime($result{$b}->{"collection_end_datetime"}, '%m/%d/%Y %H:%M')->epoch <=> Time::Piece->strptime($result{$a}->{"collection_end_datetime"}, '%m/%d/%Y %H:%M')->epoch} keys %result;
			$ordered_uids{"$county"}->{"$locid"}->{"$targ"} = \@a;
	#			# Schwartzian transform makes this sort more efficient, but I can't get it to work.
	#			# It is not particularly slow anyway. :-)
	# 		my @sorted_ids = map { $_->[1] }
	# 				sort { $a->{"collection_end_datetime"} <=> $b->{"collection_end_datetime"} }
	# 				map {
	# 						[ Time::Piece->strptime( $resultRef->{$_}->{"collection_end_datetime"}, '%m/%d/%Y %H:%M' )->epoch, $_ ]
	# 				}
	# 				keys %$resultRef;
		}
	}
}

#print Dumper(\%ordered_uids);
#die;

my %calc = ();
foreach my $county (keys %results) {
	$calc{"$county"} = {};
	foreach my $locid (keys %{$results{"$county"}}) {
		$calc{"$county"}->{"$locid"} = {};
		foreach my $targ (keys %{$results{"$county"}->{"$locid"}}) {
			$calc{"$county"}->{"$locid"}->{"$targ"} = [];
			my $loop_count = $WINDOW;
			if (defined $sentinels{"$locid"}) {
				$loop_count = $sentinels{"$locid"};
			}
			$loop_count = scalar @{$ordered_uids{"$county"}->{"$locid"}->{"$targ"}} if $loop_count > scalar @{$ordered_uids{"$county"}->{"$locid"}->{"$targ"}};
			for (my $i=0; $i<$loop_count; $i++) {
				if (defined $ordered_uids{"$county"}->{"$locid"}->{"$targ"}->[$i]) {
					my $uid = $ordered_uids{"$county"}->{"$locid"}->{"$targ"}->[$i];
					push @{$calc{"$county"}->{"$locid"}->{"$targ"}}, $results{"$county"}->{"$locid"}->{"$targ"}->{"$uid"};
					# Normalize to per person
					unless ($results{"$county"}->{"$locid"}->{"$targ"}->{"$uid"}->{"target_copies_fn_per_cap"} eq "NA") {
						$results{"$county"}->{"$locid"}->{"$targ"}->{"$uid"}->{"target_copies_fn_per_cap"} = $results{"$county"}->{"$locid"}->{"$targ"}->{"$uid"}->{"target_copies_fn_per_cap"} / $results{"$county"}->{"$locid"}->{"$targ"}->{"$uid"}->{"target_per_capita_basis"};
					}

				} else {
					print Dumper($ordered_uids{"$county"}->{"$locid"}->{"$targ"});
					die;
				}
			}
		}
	}
}
#print Dumper(\%calc);
#die;

my (%county_rollup, %state_rollup) = ((), ());

foreach my $targ (keys %targ2disease) {
	$state_rollup{"$targ"} = {};
}

print "region_name\ttarget\tnormalization\tlatest_abundance\trecent_trend\n";

foreach my $county (keys %calc) {
	$county_rollup{"$county"} = {};
	foreach my $locid (keys %{$calc{"$county"}}) {
		foreach my $targ (keys %{$calc{"$county"}->{"$locid"}}) {
			#next unless "$targ" eq "SARS-CoV-2";
			
			$county_rollup{"$county"}->{"$targ"} = {};
			# No abundance values for this location & target.
			if (scalar @{$calc{"$county"}->{"$locid"}->{"$targ"}} == 0) {
				print "$locid\t$targ\tnone\tNA\tNA\n";
#				push @{$county_rollup{"$county"}->{"$targ"}->{"abundance"}}, "NA";
#				push @{$county_rollup{"$county"}->{"$targ"}->{"trend"}}, "NA";
#				push @{$county_rollup{"$county"}->{"$targ"}->{"normalization"}}, "NA";
				next;
			}
			
			my $valkey = "target_copies_fn_per_cap";	# Default value is the normalized value.
			my $norm = "time, flow, population";
			foreach my $n (@{$calc{"$county"}->{"$locid"}->{"$targ"}}) {
				# If any normalized values are NA, use raw data instead.
				if ("$n->{target_copies_fn_per_cap}" eq "NA") {
					$valkey = "target_copies_per_ld";
					$norm = "time";
					last;
				}
			}
			
			# Only one abundance value for this location & target.
			if (scalar @{$calc{"$county"}->{"$locid"}->{"$targ"}} == 1) {
				my $abund = $calc{"$county"}->{"$locid"}->{"$targ"}->[0]->{"$valkey"};
				my $week  = $calc{"$county"}->{"$locid"}->{"$targ"}->[0]->{"collection_week"};
				print "$locid\t$targ\t$norm\t$abund\tNA\n";
				if ($resources{"location"}->{"$locid"}->{"location_category"} eq "wwtp") {
					$county_rollup{"$county"}->{"$targ"}->{$week} = [] unless defined $county_rollup{"$county"}->{"$targ"}->{$week};
					push @{$county_rollup{"$county"}->{"$targ"}->{$week}}, $abund;
					$state_rollup{"$targ"}->{$week} = [] unless defined $state_rollup{"$targ"}->{$week};
					push @{$state_rollup{"$targ"}->{$week}}, $abund;
				}
				next;
			}
			
			# More than one abundance value, look for spiking or despiking trend.
		 	# SPIKING: most recent sample is >SPIKE_THRESHOLD higher than the previous sample.
			# DESPIKING: most recent sample is >SPIKE_THRESHOLD lower than the previous sample.
			# Actual values less than 1 are not eligible for spiking/despiking.
			#
			my $ultim_abund  = $calc{"$county"}->{"$locid"}->{"$targ"}->[0]->{"$valkey"};
			my $penult_abund = $calc{"$county"}->{"$locid"}->{"$targ"}->[1]->{"$valkey"};
			my $spike_result = checkForSpike($ultim_abund, $penult_abund);	# returns 1 for spike, -1 for despike, 0 for neither
			
			my $this_sample = $calc{"$county"}->{"$locid"}->{"$targ"}->[0];
			my $this_abund = $this_sample->{"$valkey"};
			my $week  = $this_sample->{"collection_week"};
			if ($resources{"location"}->{"$locid"}->{"location_category"} eq "wwtp") {
				$county_rollup{"$county"}->{"$targ"}->{$week} = [] unless defined $county_rollup{"$county"}->{"$targ"}->{$week};
				push @{$county_rollup{"$county"}->{"$targ"}->{$week}}, $this_abund;
				$state_rollup{"$targ"}->{$week} = [] unless defined $state_rollup{"$targ"}->{$week};
				push @{$state_rollup{"$targ"}->{$week}}, $this_abund;
			}

			my @pct_diffs = ();
			for (my $i=1; $i<scalar(@{$calc{"$county"}->{"$locid"}->{"$targ"}}); $i++) {
				my $prev_sample = $calc{"$county"}->{"$locid"}->{"$targ"}->[$i];
				my $prev_abund = $prev_sample->{"$valkey"};
				my $week  = $prev_sample->{"collection_week"};

				push @pct_diffs, (100 * ($this_abund - $prev_abund) / ($prev_abund+1));
				if ($resources{"location"}->{"$locid"}->{"location_category"} eq "wwtp") {
					$county_rollup{"$county"}->{"$targ"}->{$week} = [] unless defined $county_rollup{"$county"}->{"$targ"}->{$week};
					push @{$county_rollup{"$county"}->{"$targ"}->{$week}}, $prev_abund;
					$state_rollup{"$targ"}->{$week} = [] unless defined $state_rollup{"$targ"}->{$week};
					push @{$state_rollup{"$targ"}->{$week}}, $prev_abund;
				}
				$this_abund = $prev_abund;
			}
			
			# Reset to most recent abundance
			$this_abund = $this_sample->{"$valkey"};
			my $trend = "NA";
			if ($spike_result == 1) {
				$trend = "SPIKING";
			} elsif ($spike_result == -1) {
				$trend = "DESPIKING";
			} else {
				$trend = getTrend(\@pct_diffs);
# 				if ("$locid" eq "WheelingWWTP-01") {
# 					print Dumper(\@pct_diffs);
# 					print "$trend\n";
# 					die;
# 				}
			}
			
			print "$locid\t$targ\t$norm\t$this_abund\t$trend\n";
		}
	}
}
#print Dumper(\%state_rollup);
#die;
#my @temp = sort {$b <=> $a} keys %{$state_rollup{"SARS-CoV-2"}};
#print Dumper(\@temp);
#die;


# Print out county rollups.
foreach my $county (keys %county_rollup) {
	foreach my $targ (keys %{$county_rollup{"$county"}}) {
		my $a = getMeanAbund($county_rollup{"$county"}->{"$targ"});
		my $t = getMeanTrend($county_rollup{"$county"}->{"$targ"});
		print "$county\t$targ\tNA\t$a\t$t\n";
	}
}

# Print the statewide rollup.
foreach my $targ (keys %targ2disease) {
	my $a = getMeanAbund($state_rollup{"$targ"});
	my $t = getMeanTrend($state_rollup{"$targ"});
	print "WV\t$targ\tNA\t$a\t$t\n";
}
exit $status;



sub checkForSpike {
	my $ult = shift;
	my $penult = shift;
	return 0 unless $ult > 0 and $penult > 0;
	
	my $val  = 100 * ($ult - $penult) / ($penult+1);

	if ($val >= $SPIKE_THRESHOLD) {
		return 1;
	} elsif ($val <= (-1 * $SPIKE_THRESHOLD)) {
		return -1;
	} else {
		return 0;
	}
	
}

sub getTrend {
	# Ordered by date, first element is % diff between most recent and second most recent.
	my $pct_diff_arr = shift;

# 	INCREASING: at least a TREND_THRESHOLD increase over each of the most recent WINDOW samples.
# 	DECREASING: at least a TREND_THRESHOLD decrease over each of the most recent WINDOW samples.
# 	VARIABLE: at least WINDOW-1 samples >TREND_THRESHOLD change in either direction.
# 	STABLE: at most WINDOW-2 samples >TREND_THRESHOLD change in either direction.
		
		my $total_count = scalar @{$pct_diff_arr};
		my $total_change = 0;
		my ($count_incr, $count_decr) = (0, 0);

		for (my $i=0; $i<scalar(@{$pct_diff_arr}); $i++) {
			$total_change += $pct_diff_arr->[$i];
			$count_incr++ if $pct_diff_arr->[$i] >= $TREND_THRESHOLD;
			$count_decr++ if $pct_diff_arr->[$i] <= (-1 * $TREND_THRESHOLD);
		}
		my $mean_change = $total_change / scalar(@{$pct_diff_arr});
		
		my $t = "INDETERMINATE";
		if ($total_change >= $TREND_THRESHOLD) {
			$t = "INCREASING";
		} elsif ($total_change <= (-1 * $TREND_THRESHOLD)) {
			$t = "DECREASING";
		} elsif ($count_incr + $count_decr == $total_count) {
			$t = "VARIABLE";
		} else {
			$t = "STABLE";
		}

# 		my $t = "INDETERMINATE";
# 		if ($count_decr == 0 and $count_incr > 0) {
# 			$t = "INCREASING";
# 		} elsif ($count_incr == 0 and $count_decr > 0) {
# 			$t = "DECREASING";
# 		} elsif ($count_incr + $count_decr == $total) {
# 			$t = "VARIABLE";
# 		} else {
# 			$t = "STABLE";
# 		}
		
	return "$t";
	
}

sub getMeanAbund {
	my $week2vals = shift;
	
	return 0 if scalar keys %{$week2vals} == 0;
	
	my @weeks_sorted = sort {$b <=> $a} keys %{$week2vals};
	my $ultim_week = $weeks_sorted[0];
	return 0 if scalar @{$week2vals->{$ultim_week}} == 0;
		
	my $sum = 0;
	foreach my $val (@{$week2vals->{$ultim_week}}) {
		$sum += $val;
	}
	my $mean = $sum / scalar @{$week2vals->{$ultim_week}};

	return $mean;
}

sub getMeanTrend {
	my $week2vals = shift;
	
	return "NA" if scalar keys %{$week2vals} == 0;
	
	my @weekly_means = ();
	my @weeks_sorted = sort {$b <=> $a} keys %{$week2vals};

	foreach my $week (@weeks_sorted) {
		my $sum = 0;
		for (my $i=0; $i < scalar(@{$week2vals->{$week}}); $i++) {
			$sum += $week2vals->{$week}->[$i];
		}
		my $mean = $sum / scalar @{$week2vals->{$week}};
		push @weekly_means, $mean;
	}
	
	my @pct_diffs = ();
	my $this_mean = $weekly_means[0];
	
	my $limit = $WINDOW;
	$limit = scalar(@weekly_means) if scalar(@weekly_means) < $WINDOW;
	
	for (my $i=1; $i < $limit; $i++) {
		my $prev_mean = $weekly_means[$i];
		push @pct_diffs, (100 * ($this_mean - $prev_mean) / ($prev_mean+1));
		$this_mean = $prev_mean;
	}
	
	my $trend = getTrend(\@pct_diffs);
	return "$trend";
}


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


# sub getTrendComplicated {
# 	# Ordered by date, first element is % diff between most recent and second most recent.
# 	my $pct_diff_arr = shift;
# 
# # 	INCREASING: at least a TREND_THRESHOLD increase over each of the most recent WINDOW samples.
# # 	DECREASING: at least a TREND_THRESHOLD decrease over each of the most recent WINDOW samples.
# # 	VARIABLE: at least WINDOW-1 samples >TREND_THRESHOLD change in either direction.
# # 	STABLE: at most WINDOW-2 samples >TREND_THRESHOLD change in either direction.
# 
# 		my ($pct_incr, $pct_decr) = (0, 0);
# 		my $base_trigger = ceil(0.6 * scalar @{$pct_diff_arr});
# 		
# 		for (my $i=0; $i<scalar(@{$pct_diff_arr}); $i++) {
# 			$pct_incr++ if $pct_diff_arr->[$i] >= $TREND_THRESHOLD;
# 			$pct_decr++ if $pct_diff_arr->[$i] <= (-1 * $TREND_THRESHOLD);
# 		}
# 		my $incr_trigger = $base_trigger;		# increasing across base trigger (number of gaps)
# 		my $decr_trigger = $base_trigger;		# decreasing across base trigger (number of gaps)
# 		my $var_trigger  = $base_trigger-1;	# either increasing or decreasing across base trigger save 1
# 		
# 		my $t = "INDETERMINATE";
# 		if ($pct_incr >= $incr_trigger) {
# 			$t = "INCREASING";
# 		} elsif ($pct_decr >= $decr_trigger) {
# 			$t = "DECREASING";
# 		} elsif ($pct_incr + $pct_decr >= $var_trigger) {
# 			$t = "VARIABLE";
# 		} else {
# 			$t = "STABLE";
# 		}
# 		
# 	return "$t";
# 	
# }

