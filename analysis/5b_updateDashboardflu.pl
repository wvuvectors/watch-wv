#! /usr/bin/env perl

use strict;
use warnings;
use DateTime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;
use Statistics::LineFit;

use Data::Dumper;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname [options]\n";
$usage   .=   "Build update table for the dashboard.\n";
$usage   .=   "\n";

my $debug = 0;

my %required_sample_keys   = ("Sample ID" => 1, 
															"Assay Target 1" => 1,
															"Assay Target 2" => 1,
															"Assay Target 1 Result (CN/L)" => 1,
															"Assay Target 2 Result (CN/L)" => 1,
															"Location" => 1, 
															"Sample Composite Start" => 1, 
															"Sample Composite End" => 1, 
															"Sample Flow (MGD)" => 1, 
															"Sample Received Date" => 1);

my %required_sampleQC_keys = ("Sample QC Check" => 1, 
															"Status: WaTCH" => 1, 
															"Sample Collection Method" => 1, 
															"Description" => 1);
															
my %required_site_keys     = ("status" => 1, 
															"type" => 1, 
															"category" => 1, 
															"group" => 1, 
															"level" => 1, 
															"county_population" => 1, 
															"collection_scheme" => 1, 
															"location_common_name" => 1, 
															"location_zipcode" => 1, 
															"location_info" => 1, 
															"longitude" => 1, 
															"latitude" => 1, 
															"population_served" => 1, 
															"counties_served" => 1, 
															"capacity_mgd" => 1);

my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');

my $WATCHFILE    = "updates/watchdb.LATEST.txt";
#my $WATCHFILE_MU = "/Users/tpd0001/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/mu_dashboard.LATEST.tsv";

my $FLUFILE    = "updates/watchdb_flu.LATEST.txt";

my $TABLEFILE    = "resources/WaTCH_Tables.xlsx";
my $FIELDSFILE   = "resources/fields_fluboard.txt";

my $UPDATEFILE   = "watch_fluboard.LATEST.txt";
my $UPDATEFILE_I = "watch_fluboard.$NOW.txt";

my $UPDATEPATH   = "/Users/tpd0001/github/watch-wv/data";
#my $BACKUPPATH   = "/Users/tpd0001/github/watch-wv/IGNORE";
my $BACKUPPATH   = "updates/fluboard";

my $LOD = 2855;

my %TDIST_1TAIL = ();
open (my $t1InFH, "<", "resources/tdist_1tail.txt") or die "Unable to open resources/tdist_1tail.txt for reading: $!\n";
my $tcount=0;
my @alphas = ();
while (my $line = <$t1InFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @values = split "\t", "$line", -1;
	if ($tcount == 0) {
		for (my $i=0; $i < scalar(@values); $i++) {
			push @alphas, $values[$i];
		}
	} else {
		my $df = $values[0];
		$TDIST_1TAIL{$df} = {} unless defined $TDIST_1TAIL{$df};
		for (my $i=1; $i < scalar(@values); $i++) {
			my $alpha = $alphas[$i];
			$TDIST_1TAIL{$df}->{$alpha} = $values[$i];
		}
	}
	$tcount++;
}
close $t1InFH;

my %TDIST_2TAIL = ();
open (my $t2InFH, "<", "resources/tdist_2tail.txt") or die "Unable to open resources/tdist_2tail.txt for reading: $!\n";
$tcount=0;
@alphas = ();
while (my $line = <$t2InFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @values = split "\t", "$line", -1;
	if ($tcount == 0) {
		for (my $i=0; $i < scalar(@values); $i++) {
			push @alphas, $values[$i];
		}
	} else {
		my $df = $values[0];
		$TDIST_2TAIL{$df} = {} unless defined $TDIST_2TAIL{$df};
		for (my $i=1; $i < scalar(@values); $i++) {
			my $alpha = $alphas[$i];
			$TDIST_2TAIL{$df}->{$alpha} = $values[$i];
		}
	}
	$tcount++;
}
close $t2InFH;


#
# Read data from WaTCHdb locations table
#
my %locations = ();

my $loc_wkbk = ReadData("$TABLEFILE", dtfmt => "mm/dd/yy");
my @loc_rows = Spreadsheet::Read::rows($loc_wkbk->[1]);

my %loc_headers = ();
for (my $i=0; $i < scalar(@{$loc_rows[0]}); $i++) {
	$loc_headers{$i} = $loc_rows[0][$i];
}

for (my $i=1; $i < scalar(@loc_rows); $i++) {
	my $name = $loc_rows[$i][0];
	next if !defined $name or isEmpty("$name") == 1;
	$locations{$name} = {};
	for (my $j=0; $j < scalar(@{$loc_rows[$i]}); $j++) {
		$locations{$name}->{"$loc_headers{$j}"} = $loc_rows[$i][$j];
	}
}
#print Dumper(\%locations);
#die;

# Make sure all required keys have some value for each location
foreach my $name (keys %locations) {
	my $location = $locations{"$name"};
	foreach my $key (keys %required_site_keys) {
		$location->{"$key"} = "NA" unless defined $location->{"$key"};
	}
}



#
# Input from latest WaTCH database.
# Use this to populate the keys of the master hash %asset2data
#
my %asset2data = ();

my @keysW = ();
my $count = 0;

if (-f "$WATCHFILE") {
	open (my $watchInFH, "<", "$WATCHFILE") or die "Unable to open WATCHFILE for reading: $!\n";
	while (my $line = <$watchInFH>) {
		chomp $line;
		next if $line =~ /^\s*$/;

		my @values = split "\t", "$line", -1;

		if ($count == 0) {
			# populate the keys array with values on the first line
			for (my $j=0; $j < scalar(@values); $j++) {
				$keysW[$j] = $values[$j];
			}
		} else {
			# first field is the asset ID
			my $asset_id = $values[0];
			# populate the hash entry for this asset with the remaining values
			$asset2data{$asset_id} = {} unless defined $asset2data{$asset_id};
			for (my $j=0; $j < scalar(@values); $j++) {
				$values[$j] = "NA" if isEmpty("$values[$j]") == 1;
				$asset2data{"$asset_id"}->{"$keysW[$j]"} = "$values[$j]" if defined $required_sample_keys{"$keysW[$j]"} or defined $required_sampleQC_keys{"$keysW[$j]"};
			}
			$asset2data{$asset_id}->{"responsible_lab"} = "Zoonotic Pathogen Monitoring Laboratory at West Virginia University";
			$asset2data{$asset_id}->{"responsible_lab_abbrev"} = "WVU ZPM";
		}
		$count++;
	}
	close $watchInFH;
}

#print Dumper(\%asset2data);
#die;

#
# Fetch input from latest Flu database.
#
my %flu2data = ();
my @f_field_names = ();
$count = 0;

if (-f "$FLUFILE") {
	open (my $watchInFH, "<", "$FLUFILE") or die "Unable to open $FLUFILE for reading: $!\n";
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
print "INFO: Finished flu data import. " . scalar(keys(%flu2data)) . " total samples present.\n";


#
# Remove FLU assets that do not pass QC
#
foreach my $asset_id (keys %flu2data) {
#	delete $flu2data{"$asset_id"} if $asset2data{"$asset_id"}->{"Sample QC Check"} =~ /Fail/i or $asset2data{"$asset_id"}->{"Sample QC Check"} eq "NA";
#	delete $flu2data{"$asset_id"} if $asset2data{"$asset_id"}->{"Status: WaTCH"} =~ /Rejected/i or $asset2data{"$asset_id"}->{"Status: WaTCH"} eq "NA";
	delete $flu2data{"$asset_id"} unless defined $asset2data{"$asset_id"}->{"Sample QC Check"} and $asset2data{"$asset_id"}->{"Sample QC Check"} eq "Pass";
#	delete $flu2data{"$asset_id"} unless defined $asset2data{"$asset_id"}->{"Status: WaTCH"} and $asset2data{"$asset_id"}->{"Status: WaTCH"} eq "Complete";
	delete $flu2data{"$asset_id"} if !defined $asset2data{"$asset_id"}->{"Sample Collection Method"} or $asset2data{"$asset_id"}->{"Sample Collection Method"} =~ /Grab/i;
	delete $flu2data{"$asset_id"} if defined $asset2data{"$asset_id"}->{"Description"} and $asset2data{"$asset_id"}->{"Description"} =~ /Ignore/i;
}

#print Dumper(\%asset2data);
#die;


#
# Make sure all required sample keys have some value for each asset
#
=cut
foreach my $asset_id (keys %flu2data) {
	next unless scalar keys %{$flu2data{"$asset_id"}} > 0;	# The delete function leaves the key so ignore empty assets.
	my $fluset = $flu2data{"$asset_id"};

	foreach my $key (keys %required_sample_keys) {
		$fluset->{"$key"} = "NA" unless defined $asset->{"$key"};

		# Re-format all dates
		if ($fluset->{"$key"} =~ /(\d{1,2})\/(\d{1,2})\/(\d{2,4}).*/) {
			$fluset->{"$key"} = "$1/$2/$3";
		} elsif ($fluset->{"$key"} =~ /(\d{1,2})\-(\d{1,2})\-(\d{2,4}).*/) {
			$fluset->{"$key"} = "$1/$2/$3";
		}
	}
}
=cut



foreach my $asset_id (keys %flu2data) {
	next unless scalar keys %{$flu2data{"$asset_id"}} > 0;	# The delete function leaves the key so ignore empty assets.

	my $asset = $flu2data{"$asset_id"};
	
	# Parse month, day, and year for sorting
	# Throw a warning if does not work!
	if (defined $asset->{"Sample Composite End"} and isEmpty($asset->{"Sample Composite End"}) == 0 and $asset->{"Sample Composite End"} =~ /(\d{1,2})\/(\d{1,2})\/(\d{2,4}).*/) {
		$flu2data{"$asset_id"}->{"_month"} = "$1";
		$flu2data{"$asset_id"}->{"_day"} = "$2";
		$flu2data{"$asset_id"}->{"_year"} = "$3";
	} else {
		warn "Asset $asset_id in FLUdb does not have a correctly formatted Sample Composite End field.";
	}
	
	#
	# Pre-calculate some values for faster dashboard use
	#
	
	$flu2data{"$asset_id"}->{"flua.load"} = "NA";
	$flu2data{"$asset_id"}->{"flub.load"} = "NA";

	$flu2data{"$asset_id"}->{"flua.loadcap"} = "NA";
	$flu2data{"$asset_id"}->{"flub.loadcap"} = "NA";

	# Raw values, for assets with unknown daily flow and population served
	if (defined $flu2data{"$asset_id"}->{"Assay Target 1 Result (CN/L)"}) {
		$flu2data{"$asset_id"}->{"Assay Target 1 Result (CN/L)"} =~ s/,//gi;
		
		$flu2data{"$asset_id"}->{"flua"} = $flu2data{"$asset_id"}->{"Assay Target 1 Result (CN/L)"};
		$flu2data{"$asset_id"}->{"flua"} = 1 + int rand($LOD) unless $flu2data{"$asset_id"}->{"flua"} eq "NA" or $flu2data{"$asset_id"}->{"flua"} > 1;
	}
	if (defined $flu2data{"$asset_id"}->{"Assay Target 2 Result (CN/L)"}) {
		$flu2data{"$asset_id"}->{"Assay Target 2 Result (CN/L)"} =~ s/,//gi;
		$flu2data{"$asset_id"}->{"flub"} = $flu2data{"$asset_id"}->{"Assay Target 2 Result (CN/L)"};
		$flu2data{"$asset_id"}->{"flub"} = 1 + int rand($LOD) unless $flu2data{"$asset_id"}->{"flub"} eq "NA" or $flu2data{"$asset_id"}->{"flub"} > 1;
	}
	if (defined $flu2data{"$asset_id"}->{"Sample Flow (MGD)"}) {
		$flu2data{"$asset_id"}->{"Sample Flow (MGD)"} =~ s/,//gi;
		$flu2data{"$asset_id"}->{"daily_flow"} = $flu2data{"$asset_id"}->{"Sample Flow (MGD)"};
	} else {
		$flu2data{$asset_id}->{"daily_flow"} = "NA";
	}
			
	# Pathogen load (copies normalized by flow) for assets with daily flow data
	unless ("$flu2data{$asset_id}->{flua}" eq "NA" or "$flu2data{$asset_id}->{daily_flow}" eq "NA") {
		$flu2data{"$asset_id"}->{"flua.load"} = $flu2data{"$asset_id"}->{"flua"} * 3.78541 * $flu2data{"$asset_id"}->{"daily_flow"};
	}
	unless ("$flu2data{$asset_id}->{flub}" eq "NA" or "$flu2data{$asset_id}->{daily_flow}" eq "NA") {
		$flu2data{"$asset_id"}->{"flub.load"} = $flu2data{"$asset_id"}->{"flub"} * 3.78541 * $flu2data{"$asset_id"}->{"daily_flow"};
	}

	# Pathogen load per person for assets with daily flow data and population served data
	if (defined $flu2data{"$asset_id"}->{"Location"} and defined $locations{$flu2data{"$asset_id"}->{"Location"}}) {
		my $loc_ref = $locations{$flu2data{"$asset_id"}->{"Location"}};
		unless ("$flu2data{$asset_id}->{'flua.load'}" eq "NA" or $loc_ref->{"population_served"} == -1) {
			$flu2data{"$asset_id"}->{"flua.loadcap"} = $flu2data{"$asset_id"}->{"flua.load"} / $loc_ref->{"population_served"};
		}
		unless ("$flu2data{$asset_id}->{'flub.load'}" eq "NA" or $loc_ref->{"population_served"} == -1) {
			$flu2data{"$asset_id"}->{"flub.loadcap"} = $flu2data{"$asset_id"}->{"flub.load"} / $loc_ref->{"population_served"};
		}
	}
			
}

#print Dumper(\%asset2data);
#die;


# Pre-calculate mean and 99% CI at each point over desired rolling windows (currently just 5 days)
# right end point is the most recent collection day for a particular locations
my %asset2sum = ();
my @scope_list = ("day5");
my @windows = (5);
my @target_list = ("flua", "flub", "flua.load", "flub.load", "flua.loadcap", "flub.loadcap");

foreach my $loc (keys %locations) {
	next if $debug == 1 and "$loc" ne "StarCityWWTP-01";
	my %hash_by_loc = ();
	foreach my $asset_id (keys %flu2data) {
		next unless scalar keys %{$flu2data{"$asset_id"}} > 0;	# The delete function leaves the key so ignore empty assets.

		next unless defined $flu2data{"$asset_id"}->{"_month"} and defined $flu2data{"$asset_id"}->{"_day"} and defined $flu2data{"$asset_id"}->{"_year"};

		if ("$loc" eq "$flu2data{$asset_id}->{Location}") {
			my ($m, $d, $y) = ($flu2data{"$asset_id"}->{"_month"}, $flu2data{"$asset_id"}->{"_day"}, $flu2data{"$asset_id"}->{"_year"});
			$hash_by_loc{$y} = {} unless defined $hash_by_loc{$y};
			$hash_by_loc{$y}->{$m} = {} unless defined $hash_by_loc{$y}->{$m};
			if (defined $hash_by_loc{$y}->{$m}->{$d}) {
				#print "WARNING!! Multiple samples assigned to $y-$m-$d for $loc. This is generally a BAD THING.\n";
				#print "          Asset $asset_id will be offset by 0.1 day to allow calculations but THIS SHOULD BE EXPLORED.\n";
				#print "          Overlapping IDs are $asset_id and $hash_by_loc{$y}->{$m}->{$d}.\n";
				#$d = "$d.1";
				#$asset2data{"$asset_id"}->{"_day"} = "$d";
				if ($flu2data{"$asset_id"}->{"Assay Target 1 Result (CN/L)"} > $flu2data{$hash_by_loc{$y}->{$m}->{$d}}->{"Assay Target 1 Result (CN/L)"}) {
					delete $flu2data{$hash_by_loc{$y}->{$m}->{$d}};
					$hash_by_loc{$y}->{$m}->{$d} = "$asset_id";
				} else {
					delete $flu2data{"$asset_id"};
				}
			} else {
				$hash_by_loc{$y}->{$m}->{$d} = "$asset_id";
			}
			print "$asset_id date hashed: $y-$m-$d.\n" if $debug == 1;
		}
	}
	#print Dumper(\%hash_by_loc);
	#die;
	my @loc_by_position  = ();
	#my @date_by_position = ();
	my @val_by_position  = ();
	
	my @sorted_years = sort {$a <=> $b} keys %hash_by_loc;
	foreach my $y (@sorted_years) {
		my @sorted_months = sort {$a <=> $b} keys %{$hash_by_loc{$y}};
		foreach my $m (@sorted_months) {
			my @sorted_days = sort {$a <=> $b} keys %{$hash_by_loc{$y}->{$m}};
			foreach my $d (@sorted_days) {
				my $asset_id = $hash_by_loc{$y}->{$m}->{$d};
				print "Sorting $asset_id!\n" if $debug == 1;
				push @loc_by_position, "$asset_id";
				#push @date_by_position, $asset2data{$asset_id}->{"Sample Composite End"};
				push @val_by_position, {"flua" => $flu2data{$asset_id}->{"flua"}, 
																"flub" => $flu2data{$asset_id}->{"flub"}, 
																"flua.load" => $flu2data{$asset_id}->{"flua.load"}, 
																"flub.load" => $flu2data{$asset_id}->{"flub.load"}, 
																"flua.loadcap" => $flu2data{$asset_id}->{"flua.loadcap"}, 
																"flub.loadcap" => $flu2data{$asset_id}->{"flub.loadcap"}};
			}
		}
	}
	#print Dumper(\@loc_by_position);
	#print Dumper(\@date_by_position);
	#print Dumper(\@val_by_position);
	#die;

	# loc_by_position stores asset ids for this location in calendar order (0=oldest)
	# vals_by_position stores target values for this location in calendar order (0=oldest)
	# want to calculate mean and 99% CI for each asset id for each target for each window (num of samples)
	#
	print "$loc\n" if $debug == 1;

	for (my $i=0; $i < scalar(@loc_by_position); $i++) {
		my $asset_id = $loc_by_position[$i];
		#print "$i\t$asset_id\t" . $asset2data{$asset_id}->{"Sample Composite End"} . "\n" if $debug == 1;
		
		$asset2sum{$asset_id} = {};
		foreach my $targ (@target_list) {
			foreach my $scope (@scope_list) {
				$asset2sum{$asset_id}->{"$targ.$scope.mean"} = "NA";
				$asset2sum{$asset_id}->{"$targ.$scope.ci"} = "NA";
			}
		}
		
		# extract values for each window
		for (my $j=0; $j<scalar(@windows); $j++) {
			my $win = $windows[$j];
			my $scope = $scope_list[$j];
			my ($start, $stop) = ($i-$win + 1, $i);
			unless ($start < 0) {
				foreach my $targ (@target_list) {
					extract("$targ", "$scope", [$start, $stop], $asset_id, \%asset2sum, \@val_by_position);
				}
			}
		}
		
	}	
	#print Dumper(\%asset2sum);
	#die if $debug == 1;
}


# read in the required output fields for the Dashboard
my @fields = ();
open (my $fieldsFH, "<", "$FIELDSFILE") or die "Unable to open $FIELDSFILE for reading: $!";
while (<$fieldsFH>) {
	chomp;
	push @fields, "$_";
}
close $fieldsFH;

# add the summary stats to header fields for output
foreach my $targ (sort @target_list) {
	foreach my $scope (sort @scope_list) {
		push @fields, "$targ.$scope.mean";
		push @fields, "$targ.$scope.ci";
	}
}


#
# Write the dashboard update file
#
open (my $outFH, ">", "$BACKUPPATH/$UPDATEFILE_I") or die "Unable to open $BACKUPPATH/$UPDATEFILE_I for writing: $!";
print $outFH join("\t", @fields) . "\n";
foreach my $asset_id (keys %flu2data) {
	next unless scalar keys %{$flu2data{"$asset_id"}} > 0;	# The delete function leaves the key so ignore empty assets.

	next unless defined $flu2data{$asset_id}->{"Sample Composite End"} and isEmpty($flu2data{$asset_id}->{"Sample Composite End"}) == 0;
	next unless defined $flu2data{$asset_id}->{"Location"} and isEmpty($flu2data{$asset_id}->{"Location"}) == 0;
	
	my $location = $flu2data{$asset_id}->{"Location"};
	print $outFH "$asset_id";
	for (my $i=1; $i < scalar(@fields); $i++) {
		my $field = $fields[$i];
		my $value = "";
		if (defined $flu2data{$asset_id}->{"$field"}) {
			$value = trim($flu2data{$asset_id}->{"$field"});
		} elsif (defined $locations{$location}->{"$field"}) {
			$value = trim($locations{$location}->{"$field"});
		} elsif (defined $asset2sum{$asset_id}->{"$field"}) {
			$value = trim($asset2sum{$asset_id}->{"$field"});
		} else {
			$value = "NA";
		}
		print $outFH "\t$value";
	}
	print $outFH "\n";
}
close $outFH;

`cp "$BACKUPPATH/$UPDATEFILE_I" "$UPDATEPATH/$UPDATEFILE"`;

exit 0;



sub extract {
	my $target = shift;
	my $scope  = shift;
	my $endpts = shift;
	my $asset_id = shift;
	my $asset2sum_ref = shift;
	my $val_by_position_ref = shift;
	
	my @vals = ();
	for (my $j=$endpts->[1]; $j>=$endpts->[0]; $j--) {
		unless ("$val_by_position_ref->[$j]->{$target}" eq "NA") {
			push @vals, $val_by_position_ref->[$j]->{"$target"};
		}
	}

	my $pos = $endpts->[0]-1;
	while (scalar(@vals) <= ($endpts->[1] - $endpts->[0]) and $pos >= 0) {
		unless ("$val_by_position_ref->[$pos]->{$target}" eq "NA") {
			push @vals, $val_by_position_ref->[$pos]->{"$target"};
		}
		$pos--;
	}
	
	if (scalar @vals > 1) {
		my ($m, $ci, $slope) = (0, 0, 0);
		print "$asset_id\t$target\t" if $debug == 1;
		if (scalar(@vals) > 0) {
			$m = mean(\@vals);
			$ci = ci(\@vals, 0.98);
			$slope = regress(\@vals);
		}
		print "\n" if $debug == 1;
		$asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.mean"}  = $m;
		$asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.ci"}    = $ci;
		#$asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.slope"} = $slope;
	} else {
		print "$asset_id: $target: $scope: $endpts->[0], $endpts->[1]. Less than 2 values for df calculation.\n" unless "$target" =~ /load/;
		#print Dumper($val_by_position_ref);
		#die;
	}
	
#	if ($debug == 1 and $asset2data{"$asset_id"}->{"Location"} eq "StarCityWWTP-01") {
#		print $asset2data{"$asset_id"}->{"Sample Composite End"} . "\t$scope.$target\t";
#		print "mean: " . $asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.mean"} . "\tCI: " . $asset2sum_ref->{"$asset_id"}->{"${target}.${scope}.ci"};
#		print "\t" . scalar(@vals) . " vals: " . join(",", @vals) . "\n";
#	}
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


sub mean {
	my $vals = shift;
	my $n = scalar(@$vals);
	my $sum = 0;
	foreach my $val (@$vals) {
		$sum += $val;
	}
	my $m = $sum/$n;
	return $m;
}

sub sd {
	my $vals = shift;

	my $n = scalar(@$vals);
	my $m = mean($vals);
	print "m=$m\t" if $debug == 1;
	my $sum = 0;
	
	foreach my $val (@$vals) {
		$sum += ($val - $m) * ($val - $m);
	}
	my $sdev = sqrt($sum/$n);
	return $sdev;
}

sub ci {
	my $vals = shift;
	my $level = shift;
	
	my $n = scalar(@$vals);
	my $df = $n-1;
	my $alpha = 0.5 * (1 - $level);	# 90% = 0.05, 95% = 0.025, 98% = 0.002
	my $t = $TDIST_2TAIL{$df}->{$alpha};
	if (!defined $t) {
		print "df: $df\t, alpha: $alpha\n";
		die;
	}
	my $sd = sd($vals);
	my $q = $sd/sqrt($n);
	my $conf = $t * $q;
	
	if ($debug == 1) {
		print "n=$n\t";
		print "a=$alpha\t";
		print "s=$sd\t";
		print "q=$q\t";
		print "t=$t\t";
		print "conf=$conf\t";
	}
		
	return $conf;
}

sub regress {
	my $vals = shift;
	
	my $lineFit = Statistics::LineFit->new();
	my @x = ();
	my @y = ();
	for (my $i=0; $i<scalar(@$vals); $i++) {
		my $val = $vals->[$i];
		push @x, $i;
		push @y, $val;
	}
	@y = reverse(@y);
	
	$lineFit->setData(\@x, \@y);
  my ($intercept, $slope) = $lineFit->coefficients();

	return $slope;
}



