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
use Time::Piece;

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

# my @abundance = ("low", "moderate", "high", "very high", "unknown");
# my @trend = ("stable", "increasing", "decreasing", "variable", "spiking", "unknown");
# my @risk_level = ("low", "moderate", "high", "very high", "unknown");


my $WINDOW = 3;

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

my %out_table = ();
# Columns in the output table (tab-delim), keyed on region_name:
# region_name (either location_id or county)
# target_id
# abundance
# trend
# risk_level
# dominant_variant
# Calculate these based on the most recent values in the results

my %loc2county = ();
foreach my $locid (keys %{$resources{"location"}}) {
	next unless $resources{"location"}->{"$locid"}->{"location_status"} eq "active";
	my $county = $resources{"location"}->{"$locid"}->{"location_counties_served"};
	$out_table{"$locid"}  = {};
	$out_table{"$county"} = {};
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
				# date format hack area
				if ("$colnames[$i]" eq "collection_end_datetime") {
					$cols[$i] =~ s/ AM$//i;
					if ("$cols[$i]" =~ /(\d+)\/(\d+)\/(\d+) (\d+):(\d+) PM$/) {
						my ($m, $d, $y, $h, $min) = ($1, $2, $3, $4, $5);
						$h += 12;
						$cols[$i] = "$m/$d/$y $h:$min";
					}
				}
				$results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}->{"$thisId"} = {} unless defined $results{"$thisCounty"}->{"$thisLoc"}->{"$thisTarg"}->{"$thisId"};
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
# To calculate trend, need most recent $WINDOW samples.
# 	INCREASING: at least a 5% increase over the most recent $WINDOW samples.
# 	DECREASING: at least a 5% decrease over the most recent $WINDOW samples.
# 	VARIABLE: at least 2 samples >5% away from the mean value over the most recent $WINDOW samples.
# 	STABLE: no samples >5% away from the mean value over the most recent $WINDOW samples.
# 	SPIKING: most recent sample is >200% higher than the previous sample.
# 	DESPIKING: most recent sample is >200% lower than the previous sample, which was a spike.
#
# To calculate abundance, need the most recent 1 sample.
# To calculate dominant variant, need the most recent sample and all variant proportions.
# To calculate risk level, need...?
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
			$loop_count = scalar @{$ordered_uids{"$county"}->{"$locid"}->{"$targ"}} if $loop_count > scalar @{$ordered_uids{"$county"}->{"$locid"}->{"$targ"}};
			for (my $i=0; $i<$loop_count; $i++) {
				if (defined $ordered_uids{"$county"}->{"$locid"}->{"$targ"}->[$i]) {
					my $uid = $ordered_uids{"$county"}->{"$locid"}->{"$targ"}->[$i];
					push @{$calc{"$county"}->{"$locid"}->{"$targ"}}, $results{"$county"}->{"$locid"}->{"$targ"}->{"$uid"};
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

# Columns in the output table (tab-delim), keyed on region_name:
# region_name (either location_id or county)
# target_id
# abundance
# trend
# risk_level
# dominant_variant

print "location_id\ttarget\tlatest_abundance\tmean_pct_change_over_window\twindow_size\tsample_count\n";

foreach my $county (keys %calc) {
	foreach my $locid (keys %{$calc{"$county"}}) {
		foreach my $targ (keys %{$calc{"$county"}->{"$locid"}}) {
			if (scalar @{$calc{"$county"}->{"$locid"}->{"$targ"}} == 0) {
				print "$locid\t$targ\tNA\tNA\t$WINDOW\t0\n";
				next;
			}

			if (scalar @{$calc{"$county"}->{"$locid"}->{"$targ"}} == 1) {
				my $abund = $calc{"$county"}->{"$locid"}->{"$targ"}->[0]->{"target_copies_fn_per_cap"};
				$abund = $calc{"$county"}->{"$locid"}->{"$targ"}->[0]->{"target_copies_per_ld"} if "$abund" eq "NA";
				print "$locid\t$targ\t$abund\tNA\t$WINDOW\t1\n";
				next;
			}
			
			my $this_sample = $calc{"$county"}->{"$locid"}->{"$targ"}->[-1];
			if ($this_sample->{"target_copies_fn_per_cap"} eq "NA") {
				my $this_val = $this_sample->{"target_copies_per_ld"};
				my $pct_diff = 0;
				for (my $i=scalar(@{$calc{"$county"}->{"$locid"}->{"$targ"}})-1; $i>=0; $i--) {
					my $next_sample = $calc{"$county"}->{"$locid"}->{"$targ"}->[$i];
					my $next_val = $next_sample->{"target_copies_per_ld"};
					$pct_diff += 100 * ($next_val - $this_val) / ($this_val+1);
					$this_val = $next_val;
				}
				$pct_diff /= scalar(@{$calc{"$county"}->{"$locid"}->{"$targ"}});
				print "$locid\t$targ\t$this_val\t$pct_diff\t$WINDOW\t" . scalar(@{$calc{"$county"}->{"$locid"}->{"$targ"}}) . "\n";
			} else {
				my $this_val = $this_sample->{"target_copies_fn_per_cap"};
				if ("$this_val" eq "NA") {
					warn "target_copies_fn_per_cap is $this_val for the first sample from $locid $targ\n";
					$status++;
					next;
				}
				my $pct_diff = 0;
				for (my $i=1; $i<scalar(@{$calc{"$county"}->{"$locid"}->{"$targ"}}); $i++) {
					my $next_sample = $calc{"$county"}->{"$locid"}->{"$targ"}->[$i];
					my $next_val = $next_sample->{"target_copies_fn_per_cap"};
					if ("$next_val" eq "NA") {
						warn "target_copies_fn_per_cap is $next_val for the $i sample from $locid $targ\n";
						$status++;
						next;
					}
					$pct_diff += 100 * ($next_val - $this_val) / ($this_val+1);
					$this_val = $next_val;
				}
				$pct_diff /= scalar(@{$calc{"$county"}->{"$locid"}->{"$targ"}});
				print "$locid\t$targ\t$this_val\t$pct_diff\t$WINDOW\t" . scalar(@{$calc{"$county"}->{"$locid"}->{"$targ"}}) . "\n";
			}
		}
	}
}

#			$val1 = $calc{"$county"}->{"$locid"}->{"$targ"}->{"$uid"}->{"target_copies_per_ld"};

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


