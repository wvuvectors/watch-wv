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
$usage   .= "Usage: $progname [options]\n";
$usage   .=  "Build update table for NWSS database using WaTCHdb table on STDIN.\n";
$usage   .=   "\n";

my $LOD = 3;
my $NTC_THRESHOLD = 3;

#
# Get time stamp for now
#
my $NOW =
   DateTime
      ->now( time_zone => 'local' )
      ->set_time_zone('floating')
      ->strftime('%Y-%m-%d.%H-%M-%S');



my $MAPFILE     = "resources/WaTCH_to_NWSS.xlsx";
my $TABLEFILE   = "resources/WaTCH_Tables.xlsx";
my $SAMPLESFILE = "updates/watchdb.LATEST.txt";
my $FIPSFILE    = "resources/fips_codes.wv.txt";

my $NWSSFILE_UPDATE   = "updates/nwss/wvu_zpm.nwss_update.$NOW.csv";

# Eventually (soon!) this script will update the main WaTCHdb file directly.
#my $SAMPLEFILE_NWSSUP = "DB/watchdb.NWSS_UP.$NOW.txt";



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
	$locations{$name} = {};
	for (my $j=0; $j < scalar(@{$loc_rows[$i]}); $j++) {
		$locations{$name}->{"$loc_headers{$j}"} = $loc_rows[$i][$j];
	}
}
#print Dumper(\%locations);


#
# Read data from WaTCH-NWSS mapping table
#
my %map = ();

my $map_wkbk = ReadData("$MAPFILE", dtfmt => "mm/dd/yy");
my @map_rows = Spreadsheet::Read::rows($map_wkbk->[1]);

my %map_headers = ();
for (my $i=0; $i < scalar(@{$map_rows[0]}); $i++) {
	$map_headers{$i} = $map_rows[0][$i];
}

for (my $i=1; $i < scalar(@map_rows); $i++) {
	my $nwss_field = $map_rows[$i][1];
	$map{$nwss_field} = {};
	for (my $j=0; $j < scalar(@{$map_rows[$i]}); $j++) {
		$map{$nwss_field}->{"$map_headers{$j}"} = $map_rows[$i][$j];
	}
}
#print Dumper(\%map);


#
# Read data from table mapping county name to FIPS code
#
my %county2fips = ();

open (my $fipsFH, "<", "$FIPSFILE") or die "Unable to open $FIPSFILE for reading: $!\n";
while (my $line = <$fipsFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;

	my ($fips, $county) = split "\t", "$line", -1;
	$county2fips{$county} = $fips;
}
close $fipsFH;
#print Dumper(\%county2fips);




#
# Read sample data from the WaTCHdb database.
#
my %samples = ();

my @watchdb_fields = ();
my $count = 0;

open (my $watchInFH, "<", "$SAMPLESFILE") or die "Unable to open $SAMPLESFILE for reading: $!\n";
while (my $line = <$watchInFH>) {
	chomp $line;
	next if $line =~ /^\s*$/;

	my @fields = split "\t", "$line", -1;

if ($count == 0) {
	# populate the keys array with values on the first line
	for (my $j=0; $j < scalar(@fields); $j++) {
		$watchdb_fields[$j] = trim($fields[$j]);
	}
} else {
		# first field is the asset ID (unique sample ID)
		my $asset_id = $fields[0];

		# populate the hash entry for this sample
		$samples{$asset_id} = {} unless defined $samples{$asset_id};
		for (my $j=0; $j < scalar(@fields); $j++) {
			$samples{"$asset_id"}->{"$watchdb_fields[$j]"} = trim($fields[$j]);
		}
	}
	$count++;
}
close $watchInFH;
#print Dumper(\%samples);


#
# Fit sample data into required NWSS fields using %map as a guide. This is...messy. :-(
# The map file has WaTCH fields (and tables) that correspond to NWSS fields, where a direct mapping occurs, but some have 
# to be hacked together and there is no good way to clean that up yet. For now I just have a separate sub for the hacks.
#
my %nwss = ();

foreach my $asset (keys %samples) {
	# Skip unless WaTCH Status is "Complete" - indicates full data processing pipeline has been run successfully
	next unless $samples{$asset}->{"Status: WaTCH"} eq "Complete";

	# Skip if this sample has been rejected
	next if $samples{$asset}->{"Status: NWSS"} =~ /Rejected/;
	# Skip if we've already submitted this sample. Currently we want a full kill-and-replace submission every time.
	#next if $samples{$asset}->{"Status: NWSS"} =~ /Submitted/;
	
	# Skip if we have no assay result for this asset
	next unless defined $samples{$asset}->{"Assay Target 1 Result (CN/L)"} and $samples{$asset}->{"Assay Target 1 Result (CN/L)"} ne "";
	# Get location entry for this asset, for easier access
	my $asset_loc = $samples{$asset}->{"Location"};
	# Initialize hash entry
	$nwss{$asset} = {};
	foreach my $nwss_field (keys %map) {
		if (defined $map{$nwss_field}->{"WaTCH_table"} and $map{$nwss_field}->{"WaTCH_table"} eq "locations") {
			# Retrieve from locations hash using embedded field
			my $wfield = $map{$nwss_field}->{"WaTCH_field"};
			$nwss{$asset}->{"$nwss_field"} = $locations{$asset_loc}->{"$wfield"};
		} elsif (defined $map{$nwss_field}->{"WaTCH_table"} and $map{$nwss_field}->{"WaTCH_table"} eq "samples") {
			# Retrieve from samples hash using embedded field
			my $wfield = $map{$nwss_field}->{"WaTCH_field"};
			$nwss{$asset}->{"$nwss_field"} = $samples{$asset}->{"$wfield"};
		} elsif ($map{$nwss_field}->{"WaTCH_field"} eq "[constant]") {
			# Constant value
			$nwss{$asset}->{"$nwss_field"} = "";
			$nwss{$asset}->{"$nwss_field"} = $map{$nwss_field}->{"constant_value"} if defined $map{$nwss_field}->{"constant_value"};
		} elsif ($map{$nwss_field}->{"WaTCH_field"} eq "[empty]") {
			# Empty value
			$nwss{$asset}->{"$nwss_field"} = "";
		} elsif ($map{$nwss_field}->{"WaTCH_field"} eq "[calculated]") {
			# Calculated value. These get complicated!
			$nwss{$asset}->{"$nwss_field"} = fill_calculated("$nwss_field", $samples{$asset}, $locations{$asset_loc});
		}
	}
	
	# Format the collection date and time appropriately. I hate date hacking. DateTime::Format::Strptime doesn't work although 
	# I have not given up hope.
	#
	if (defined $nwss{$asset}->{"sample_collect_date"} and $nwss{$asset}->{"sample_collect_date"} ne "") {
		my $rawstamp = $nwss{$asset}->{"sample_collect_date"};
		$rawstamp =~ s/ AM//gi;
		if (my ($mon, $day, $yr, $hrs, $min) = $rawstamp =~ /(.+?)\/(.+?)\/(.+?) (.+?):(.+?)/) {
			$yr = "20$yr" if length $yr == 2;
			$mon = "0$mon" if length $mon == 1;
			$day = "0$day" if length $day == 1;
			my $cdate = "$yr-$mon-$day";
			$hrs = "0$hrs" if length $hrs == 1;
			$min = "0$min" if length $min == 1;
			my $ctime = "$hrs:$min";
			$nwss{$asset}->{"sample_collect_date"} = $cdate;
			$nwss{$asset}->{"sample_collect_time"} = $ctime;
		}
	}

	# Format the test result date appropriately. I hate date hacking. DateTime::Format::Strptime doesn't work although 
	# I have not given up hope.
	#
	if (defined $nwss{$asset}->{"test_result_date"} and $nwss{$asset}->{"test_result_date"} ne "") {
		my $rawstamp = $nwss{$asset}->{"test_result_date"};
		$rawstamp =~ s/ AM//gi;
		if (my ($mon, $day, $yr) = $rawstamp =~ /(.+?)\/(.+?)\/(.+)/) {
			$yr = "20$yr" if length $yr == 2;
			$mon = "0$mon" if length $mon == 1;
			$day = "0$day" if length $day == 1;
			my $cdate = "$yr-$mon-$day";
			$nwss{$asset}->{"test_result_date"} = $cdate;
		}
	}

}
#print Dumper(\%nwss);
#die;



# Output csv format to NWSS update file

my @nwss_fields = sort keys %map;

open (my $nwssoutFH, ">", "$NWSSFILE_UPDATE") or die "Unable to open $NWSSFILE_UPDATE for writing: $!";
print $nwssoutFH "\"" . join("\",\"", @nwss_fields) . "\"\n";
foreach my $sample_id (keys %nwss) {
	my %watch_sample = %{$samples{$sample_id}};
	my $valid_check = validate_nwss(\%watch_sample);
	if ($valid_check == 1) {
		$watch_sample{"Status: NWSS"} = "Validated";

		my %nwss_sample  = %{$nwss{$sample_id}};

		for (my $i=0; $i < scalar @nwss_fields; $i++) {
			my $field = $nwss_fields[$i];
			my $entry = "";
			$entry = "$nwss_sample{$field}" if defined $nwss_sample{$field};
			print $nwssoutFH "," unless $i == 0;
			print $nwssoutFH "\"$entry\"";
		}
		print $nwssoutFH "\n";
	} elsif ($valid_check == -1) {
		$watch_sample{"Status: NWSS"} = "Rejected";
	} else {
		$watch_sample{"Status: NWSS"} = "Processing";
	}
	
}
close $nwssoutFH;


# Write entire WaTCHdb file with updated Status: NWSS values. Change "Validated" to "Submitted"
#
=cut
open (my $watchoutFH, ">", "$SAMPLEFILE_NWSSUP") or die "Unable to open $SAMPLEFILE_NWSSUP for writing: $!";
print $watchoutFH join("\t", @watchdb_fields) . "\n";
foreach my $sample_id (keys %samples) {
	my %sample = %{$samples{$sample_id}};
	my $colnum = 0;
	foreach my $field (@watchdb_fields) {
		my $entry = "";
		$entry = "$sample{$field}" if defined $sample{$field};
		if ("Status: NWSS" eq "$field" and "Validated" eq "$entry") {
			$entry = "Submitted $NOW";
		}
		print $watchoutFH "\t" unless $colnum == 0;
		print $watchoutFH "$entry";
		$colnum++;
	}
	print $watchoutFH "\n";
}
close $watchoutFH;
=cut


exit 0;


sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


sub fill_calculated {
	#
	# Sub-routine to make the messy connections for [calculated] fields in the WaTCH-NWSS conversion. Eventually we will 
	# add NWSS fields to WATCH for direct mapping but lower priority than getting usable data out.
	#
	my $field = shift;
	my $sample_ref = shift;
	my $location_ref = shift;
	
	my $return_val = "";
	
	if ("$field" eq "rec_eff_target_name") {
		# samples.Assay Control, Process or empty if that is not defined
		#$return_val = $sample_ref->{"Assay Control, Process"} if defined $sample_ref->{"Assay Control, Process"} and $sample_ref->{"Assay Control, Process"} !~ /^\s*$/;
		$return_val = "";
	} elsif ("$field" eq "county_names") {
		# locations.counties_served converted to FIPS code
		#unless (defined $location_ref->{"counties_served"}) {
		# print Dumper($location_ref);
		# die;
		#}
		my @county_list = split /\s*,\s*/, $location_ref->{"counties_served"}, -1;
		my @fips_list = ();
		foreach my $county (@county_list) {
			push @fips_list, $county2fips{"$county"};
		}
		#print Dumper(\@fips_list);
		$return_val = join(",", @fips_list);
	} elsif ("$field" eq "sample_location_specify") {
		# locations.name if locations.category eq upstream; otherwise empty
		$return_val = $location_ref->{"name"} if defined $location_ref->{"category"} and $location_ref->{"category"} eq "upstream";
	} elsif ("$field" eq "rec_eff_spike_matrix") {
		# "raw sample" or empty if samples.Assay Control, Process is empty
		#$return_val = "raw sample" if defined $sample_ref->{"Assay Control, Process"} and $sample_ref->{"Assay Control, Process"} !~ /^\s*$/;
		$return_val = "";
	} elsif ("$field" eq "sars_cov2_below_lod") {
		# yes if samples.Assay Target 1 Concentration (CN/Rxn) < $LOD; otherwise no
		$return_val = "no";
		$return_val = "yes" if !defined $sample_ref->{"Assay Target 1 Concentration (CN/Rxn)"} or $sample_ref->{"Assay Target 1 Concentration (CN/Rxn)"} =~ /^\s*$/ or $sample_ref->{"Assay Target 1 Concentration (CN/Rxn)"} < $LOD;
	} elsif ("$field" eq "ntc_amplify") {
		# yes if samples.Assay Target 1 PCR NC Result (CN/Rxn) > $NTC_THRESHOLD; otherwise no
		$return_val = "yes";
		$return_val = "no" if !defined $sample_ref->{"Assay Target 1 PCR NC Result (CN/Rxn)"} or $sample_ref->{"Assay Target 1 PCR NC Result (CN/Rxn)"} =~ /^\s*$/ or $sample_ref->{"Assay Target 1 PCR NC Result (CN/Rxn)"} < $NTC_THRESHOLD;
	} elsif ("$field" eq "rec_eff_percent") {
		# -1 if if samples.Assay Control, Process is empty; 100 * (samples.Assay Control, Process Result (CN/L) / samples.Assay Control, Process Spike-In (CN/L))
		if (defined $sample_ref->{"Assay Control, Process"} and $sample_ref->{"Assay Control, Process"} !~ /^\s*$/ and $sample_ref->{"Assay Control, Process Result (CN/L)"} !~ /^\s*$/ and $sample_ref->{"Assay Control, Process Spike-In (CN/L)"} !~ /^\s*$/) {
			$return_val = 100 * ($sample_ref->{"Assay Control, Process Result (CN/L)"} / $sample_ref->{"Assay Control, Process Spike-In (CN/L)"});
		} else {
			$return_val = -1;
		}
	}	
	
	return $return_val;
}


sub validate_nwss {

	my $wsample_ref = shift;
	my $valid = 1;
	
	# Ignore any entries where Status: WaTCH is not "Complete"
	$valid = 0 unless $wsample_ref->{"Status: WaTCH"} eq "Complete";

	# Ignore any entries where Status: NWSS contains "Submitted"
	$valid = 0 if $wsample_ref->{"Status: NWSS"} =~ /Submitted/;

	# Ignore if sample has no flow data
	$valid = 0 if !defined $wsample_ref->{"Sample Flow (MGD)"} or $wsample_ref->{"Sample Flow (MGD)"} =~ /^\s*$/;

	# Ignore if sample has no target result
	$valid = 0 if !defined $wsample_ref->{"Assay Target 1 Concentration (CN/Rxn)"} or $wsample_ref->{"Assay Target 1 Concentration (CN/Rxn)"} =~ /^\s*$/;
	
	# Ignore if sample QC check failed
	$valid = -1 if !defined $wsample_ref->{"Sample QC Check"} or $wsample_ref->{"Sample QC Check"} =~ /^\s*$/ or $wsample_ref->{"Sample QC Check"} =~ /Fail/;

	return $valid;
}

