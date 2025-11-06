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


my $WINDOW = 4;	# Number of weeks to use in trend calculations.
my $SPIKE_THRESHOLD = 500; # Trigger (percent change) for identifying spikes/despikes.
my $TREND_THRESHOLD = 20; # Trigger (percent change) for identifying consistent trend; given as percent change.

my %targ2disease = (
	"Influenza Virus A (FluA)" => "FLUA", 
	"Influenza Virus B (FluB)" => "FLUB", 
	"Human Norovirus GII (HuNoV-GII)" => "NoV", 
	"SARS-CoV-2" => "COVID", 
	"Respiratory Syncitial Virus, Human (RSV)" => "RSV"
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
			
			# skip any legacy SARS-CoV-2 N1 results
			next if "$thisTarg" eq "SARS-CoV-2" and "$thisLocus" eq "N1";
			
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
#die;

my %ordered_means = ();
 
foreach my $epiweek (sort {$b <=> $a} keys %results) {
	$ordered_means{"WV"} = {} unless defined $ordered_means{"WV"};
	my %state_vals = ();
	foreach my $county (keys %{$results{$epiweek}}) {
		$ordered_means{"$county"} = {} unless defined $ordered_means{"$county"};
		my %county_vals = ();
		foreach my $loc (keys %{$results{$epiweek}->{"$county"}}) {
			$ordered_means{"$loc"} = {} unless defined $ordered_means{"$loc"};
			foreach my $targ (keys %{$results{$epiweek}->{"$county"}->{"$loc"}}) {
				$state_vals{"$targ"} = [] unless defined $state_vals{"$targ"};
				$county_vals{"$targ"} = [] unless defined $county_vals{"$targ"};
				$ordered_means{"$county"}->{"$targ"} = [] unless defined $ordered_means{"$county"}->{"$targ"};
				$ordered_means{"$loc"}->{"$targ"} = [] unless defined $ordered_means{"$loc"}->{"$targ"};
				my $mean = calcSimpleMean($results{$epiweek}->{"$county"}->{"$loc"}->{"$targ"});
				push @{$ordered_means{"$loc"}->{"$targ"}}, $mean;
				push @{$county_vals{"$targ"}}, $mean unless defined $do_not_rollup{"$loc"};
				push @{$state_vals{"$targ"}}, $mean unless defined $do_not_rollup{"$loc"};
			}
		}
		# county rollup
		foreach my $targ (keys %county_vals) {
			my $county_mean = calcSimpleMean($county_vals{"$targ"});
			push @{$ordered_means{"$county"}->{"$targ"}}, $county_mean;
		}
	}
	# state rollup
	foreach my $targ (keys %state_vals) {
		my $state_mean = calcSimpleMean($state_vals{"$targ"});
		push @{$ordered_means{"WV"}->{"$targ"}}, $state_mean;
	}
}
#print Dumper(\%ordered_means);
#die;


print "region_name\ttarget\tabundance_pct_change\ttrend\n";
foreach my $region (keys %ordered_means) {
	foreach my $targ (keys %{$ordered_means{"$region"}}) {
		my $this_arr = $ordered_means{"$region"}->{"$targ"};
		my $total = scalar @{$this_arr};
		# No abundance values for this region.
		if ($total <= 1) {
			print "$region\t$targ\tNA\tNA\n";
		} else {
			# More than one abundance value: try to calculate risk.
			my $comp = 13;
			$comp = scalar @{$this_arr} if scalar @{$this_arr} < 13;
			my $pct_ch = "NA";
			if ($comp > 3) {
				my $sum = 0;
				for (my $i=0; $i < $comp; $i++) {
					$sum += $this_arr->[$i];
				}
				my $mean_comp = $sum / $comp;
				if ($mean_comp == 0) {
					warn "The sum of the ordered means for target '$targ' from region '$region' is 0!";
					print "$region\t$targ\tNA\tNA\n";
					next;
				} else {
					$pct_ch = 100 * ($this_arr->[0] - $mean_comp)/$mean_comp;
				}
			}
			
			# More than one abundance value: look for trend.
			my $trend = "NA";
			my $prnt_str = "";
			
			# First look for spiking or despiking trend.
		 	# SPIKING: most recent sample is >SPIKE_THRESHOLD higher than the previous sample.
			# DESPIKING: most recent sample is >SPIKE_THRESHOLD lower than the previous sample.
			# Actual values less than 1 are not eligible for spiking/despiking.
			#
			my $ultim_ab  = $this_arr->[0];
			my $penult_ab = $this_arr->[1];
			my $spike_result = checkForSpike($ultim_ab, $penult_ab);	# returns 1 for spike, -1 for despike, 0 for neither
			
			if ($spike_result == 1) {
				$trend = "SPIKING";
				$prnt_str .= "$ultim_ab,$penult_ab";
			} elsif ($spike_result == -1) {
				$trend = "DESPIKING";
				$prnt_str .= "$ultim_ab,$penult_ab";
			} else {
				# Look at the pct diff between values for the trend.
				my @abundances = ();
				for (my $i=0; $i<$WINDOW; $i++) {
					last if $i >= scalar(@{$this_arr});
					push @abundances, $this_arr->[$i];
					$prnt_str .= "$this_arr->[$i],";
				}
				$trend = getTrend(\@abundances);
			}
#			print "$region\t$targ\t$prnt_str\t$trend\n";
			print "$region\t$targ\t$pct_ch\t$trend\n";
		}
	}
}

exit $status;



sub checkForSpike {
	my $ult = shift;
	my $penult = shift;
	return 0 unless $ult > 1 and $penult > 1;
	
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
	# Array of abundance values, ordered chronologically, element 0 is most recent.
	my $a = shift;

# 	INCREASING: at least a TREND_THRESHOLD increase over each of the most recent WINDOW samples.
# 	DECREASING: at least a TREND_THRESHOLD decrease over each of the most recent WINDOW samples.
# 	VARIABLE: at least WINDOW-1 samples >TREND_THRESHOLD change in either direction.
# 	STABLE: at most WINDOW-2 samples >TREND_THRESHOLD change in either direction.
		
		my ($incr, $decr, $stable) = (0,0,0);
		
		for (my $i=scalar(@{$a})-1; $i>0; $i--) {
			my ($a_now, $a_next) = ($a->[$i], $a->[$i-1]);
			my $pct_diff = (100 * ($a_next - $a_now) / ($a_now+1));
			if ($pct_diff >= $TREND_THRESHOLD) {
				$incr++;
			} elsif ($pct_diff <= (-1 * $TREND_THRESHOLD)) {
				$decr++;
			} else {
				$stable++;
			}
		}
		
		my $t = "INDETERMINATE";
		if ($incr == 0 and $decr == 0) {
			$t = "STABLE";
		} elsif ($incr - $decr > 0) {
			$t = "INCREASING";
		} elsif ($decr - $incr > 0) {
			$t = "DECREASING";
		} else {
			$t = "VARIABLE";
		}

	return "$t";
	
}


sub calcSimpleMean {
	my $a = shift;
	
	return 0 if scalar @{$a} == 0;
	
	my $sum = 0;
	foreach my $val (@{$a}) {
		$sum += $val;
	}
	my $mean = $sum / scalar @{$a};

	return $mean;
}




