#! /usr/bin/env perl

use strict;
use warnings;

use Text::CSV qw(csv);

use Date::Format;
use DateTime qw( );
use DateTime::Format::Strptime qw( );

use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX;
use DateTime::Format::Excel;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname MU_FILE\n";
$usage   .=  "Extracts the MU data from the csv in MU_FILE, converts to watchdb result and sample ";
$usage   .=  "table format, and writes to the dashboard directory.\n";
$usage   .=   "\n";


my $status = 0;
my $muf;

while (@ARGV) {
  my $arg = shift;
  if ("$arg" eq "-h") {
		die "$usage";
	} else {
		$muf = "$arg";
	}
}

die "A file path must be provided to $progname!\n\n$usage\n" unless defined $muf;
die "$muf is not a readable file!\n\n$usage\n" unless -f "$muf";


my %mu_data = ();

my %nwss2watch = (
	"sars-cov-2" => "SARS-CoV-2",
	"n1"				 => "N1:SARS",
	"n2"				 => "N2:SARS",
	"FLUAV"			 => "Flu A",
	"FLUBV"			 => "Flu B",
	"InfA"			 => "M:FLUA",
	"InfB"			 => "NEP/NS1:FLUB"
);

my %reqd_hdrs = (
"flow_rate" => 1,
"lab_id" => 1,
"ntc_amplify" => 1,
"pcr_gene_target" => 1,
"pcr_target" => 1,
"pcr_target_avg_conc" => 1,
"pcr_target_units" => 1,
"pcr_type" => 1,
"sample_collect_date" => 1,
"sample_collect_time" => 1,
"sample_id" => 1,
"sample_type" => 1,
"site_id" => 1,
"test_result_date" => 1,
"wwtp_name" => 1
);


my $csvA  = Text::CSV->new({auto_diag => 4, binary => 1});
my $linenum = 0;
my @colheaders = ();
my $labindx = -1;
open(my $fh, "<", "$muf") or die "Unable to open $muf for reading: $!";

while (my $line = $csvA->getline($fh)) {
	my @cols = @$line;
	# Set the file colheaders if this is the first line.
	# These are used as keys of the mu_data sub-hashes.
	if ($linenum == 0) {
		foreach (@cols) {
			$labindx = scalar(@colheaders) if trim("$_") eq "lab_id";
			push @colheaders, trim("$_");
		}
	} else {
		die "Unable to parse a lab_id column from $muf in $progname!\n" unless $labindx > -1;
		if ("MUIDSL" eq "$cols[$labindx]") {
			my $uid = "mu" . "$linenum";
			$mu_data{"$uid"} = {};
			for (my $i=0; $i < scalar(@cols); $i++) {
				my $val = $cols[$i];
				$mu_data{"$uid"}->{"$colheaders[$i]"} = "$val" if defined $reqd_hdrs{"$colheaders[$i]"};
			}
		}
	}
	$linenum++;
}
close $fh;

#print Dumper(\%mu_data);
#die;



my %locations = ();
#
# Read WaTCHdb resources from resources/ dir.
#
# Read location table into hash
my $resource_wkbk = ReadData("resources/watchdb.all_tables.xlsx", dtfmt => "mm/dd/yy");
#print Dumper($resource_wkbk);
#die;

foreach my $sheet_name (keys %{$resource_wkbk->[0]->{"sheet"}}) {
#	print "$sheet_name\n";
	next unless "$sheet_name" eq "location";
	my @field_names = ();

	my $pos = $resource_wkbk->[0]->{"sheet"}->{"$sheet_name"};
	my $sheetRef = $resource_wkbk->[$pos];
	
	# The field labels are in row A (or, index 1 of each array in 'cell' ). Store these in an array.
	for (my $i=1; $i < scalar(@{$sheetRef->{"cell"}}); $i++) {
		#print "$sheetRef->{cell}->[$i]->[1]\n";
		next unless defined $sheetRef->{"cell"}->[$i]->[1] and "$sheetRef->{cell}->[$i]->[1]" ne "";
		my $col_name = trim("$sheetRef->{cell}->[$i]->[1]");
		push @field_names, "$col_name";
	}
	#print Dumper(\@field_names);
	#last;

	# The keys for each resource sheet are in column E (or, array 5 of 'cell'). 
	# Loop over these to add entries to the resources hash.
	for (my $i=2; $i < scalar(@{$sheetRef->{"cell"}->[5]}); $i++) {
		#print "$sheetRef->{cell}->[1]->[$i]\n";
		next unless defined $sheetRef->{"cell"}->[5]->[$i] and "$sheetRef->{cell}->[5]->[$i]" ne "";
		my $key = trim("$sheetRef->{cell}->[5]->[$i]");
		$locations{"$key"} = {};
	}
	#last;
	
	# To get the actual values, loop over the 'cell' array.
	# Get the key from array 1.
	# Loop over the field_names array to get the values.
	for (my $i=0; $i < scalar(@{$sheetRef->{"cell"}})-1; $i++) {
		next unless defined $sheetRef->{"cell"}->[$i+1] and defined $field_names[$i] and "$field_names[$i]" ne "";
		my $colRef = $sheetRef->{"cell"}->[$i+1];
		my $field = "$field_names[$i]";
		# Loop over the values in this column. Add each to the resources hash.
		for (my $j=1; $j < scalar(@{$colRef}); $j++) {
			next unless defined $sheetRef->{"cell"}->[5]->[$j] and "$sheetRef->{cell}->[5]->[$j]" ne "";
			my $key = "$sheetRef->{cell}->[5]->[$j]";
			my $value = "";
			$value = trim("$colRef->[$j]") if defined $colRef->[$j];
			$locations{"$key"}->{"$field"} = "$value";
		}
	}
}
#print Dumper(\%locations);
#die;

# write relevant info to sample table
#
my @sample_hdr = (
	"sample_id", 
	"sample_status", 
	"location_id", 
	"sample_event", 
	"sample_qc", 
	"sample_collection_start_datetime", 
	"sample_collection_end_datetime", 
	"sample_recovered_datetime", 
	"sample_collection_by", 
	"sample_flow", 
	"sample_received_by", 
	"sample_received_date", 
	"sample_ph_lab", 
	"sample_comment"
);

open(my $SFH, ">", "../dashboard/data/mu.sample.txt") or die "Unable to open sample file for writing: $!";

print $SFH join("\t", @sample_hdr) . "\n";
foreach my $uid (keys %mu_data) {
	my $locRef = $locations{"$mu_data{$uid}->{wwtp_name}"};
	if (defined $locRef) {
		my %thisd = %{$mu_data{"$uid"}};
		my $prnt = "";
		$prnt .= "$thisd{sample_id}\tImported\t$locRef->{location_id}\tRoutine Surveillance\tPass\t";

		# start and end datetimes as MM/DD/YYYY HH:MM in 24-hr format
		next unless ($thisd{"sample_collect_date"} =~ /^\d+?\/\d+?\/\d+$/ or $thisd{"sample_collect_date"} =~ /^\d+-\d+-\d+$/);
		next unless ($thisd{"sample_collect_time"} =~ /^\d+:\d+/);
		my $interval_hrs = $locRef->{"location_collection_window_hrs"};
		my $sdt = makeDT($thisd{"sample_collect_date"}, $thisd{"sample_collect_time"});
		$prnt .= $sdt->mdy("/") . " " . $sdt->hour() . ":" . $sdt->minute() . "\t";
		$sdt->add(hours => $interval_hrs);
		$prnt .= $sdt->mdy("/") . " " . $sdt->hour() . ":" . $sdt->minute() . "\t";
		# use end datetime as sample_recovered_datetime
		$prnt .= $sdt->mdy("/") . " " . $sdt->hour() . ":" . $sdt->minute() . "\t";

		$prnt .= "MU\t$thisd{flow_rate}\tMU\t";

		# sample received date as MM/DD/YYYY
		$prnt .= $sdt->mdy("/") . "\t";

		$prnt .= "\tImported from NWSS file\n";
		print $SFH "$prnt";
	}
}
close $SFH;


# write relevant info to result table
#
my @result_hdr = (
	"assay_id", 
	"sample_id", 
	"collection_start_datetime", 
	"collection_end_datetime", 
	"event_type", 
	"sample_flow", 
	"sample_qc", 
	"location_id", 
	"target", 
	"target_genetic_locus", 
	"lab_id", 
	"target_copies_per_l", 
	"target_copies_per_ld", 
	"target_copies_per_ldcap", 
	"target_copies_flownorm", 
	"target_copies_fn_per_cap", 
	"target_per_capita_basis", 
	"nc_copies_per_rxn", 
	"pc_copies_per_rxn", 
	"target_result_validated"
);

open(my $RFH, ">", "../dashboard/data/mu.result.txt") or die "Unable to open result file for writing: $!";

print $RFH join("\t", @result_hdr) . "\n";
foreach my $uid (keys %mu_data) {
	my $locRef = $locations{"$mu_data{$uid}->{wwtp_name}"};
	if (defined $locRef) {
		my %thisd = %{$mu_data{"$uid"}};
		my $prnt = "";
		$prnt .= "$thisd{sample_id}\t$thisd{sample_id}\t";

		# start and end datetimes as MM/DD/YYYY HH:MM in 24-hr format
		next unless ($thisd{"sample_collect_date"} =~ /^\d+?\/\d+?\/\d+$/ or $thisd{"sample_collect_date"} =~ /^\d+-\d+-\d+$/);
		next unless ($thisd{"sample_collect_time"} =~ /^\d+:\d+/);
		my $interval_hrs = $locRef->{"location_collection_window_hrs"};
		my $sdt = makeDT($thisd{"sample_collect_date"}, $thisd{"sample_collect_time"});
		$prnt .= $sdt->mdy("/") . " " . $sdt->hour() . ":" . $sdt->minute() . "\t";
		$sdt->add(hours => $interval_hrs);
		$prnt .= $sdt->mdy("/") . " " . $sdt->hour() . ":" . $sdt->minute() . "\t";

		$prnt .= "Routine Surveillance\t$thisd{flow_rate}\tPass\t$locRef->{location_id}\t";

		# target & target_genetic_locus need to be looked up
		$prnt .= "$nwss2watch{$thisd{pcr_target}}\t$nwss2watch{$thisd{pcr_gene_target}}\t";

		# copies per ld, ldcap, flownorm, fn_per_cap, and target_per_capita_basis need to be calculated
		my $copies_per_ld    = $thisd{"pcr_target_avg_conc"} / ($interval_hrs/24);
		my $copies_per_ldcap = calcCopiesPop($copies_per_ld, 100, $locRef);
		my $copies_flownorm  = calcCopiesFlowNorm($copies_per_ld, $thisd{"flow_rate"}, $locRef);
		my $copies_fn_per_cap = calcCopiesPop($copies_flownorm, 100, $locRef);
		$prnt .= "$thisd{lab_id}\t$thisd{pcr_target_avg_conc}\t";
		$prnt .= "$copies_per_ld\t$copies_per_ldcap\t$copies_flownorm\t$copies_fn_per_cap\t";
		$prnt .= "100\t0\tNA\tyes\n";

		print $RFH "$prnt";
	}
}
close $RFH;


exit $status;


sub trim {
	my $val = shift;
	
	my $trimmed = $val;
	$trimmed =~ s/ +$//;
	$trimmed =~ s/^ +//;
	
	return $trimmed;
}


sub makeDT {
	my $dstr = shift;
	my $tstr = shift;
	
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
	
	my $dtstr = "$dstr $tstr";
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


sub calcCopiesPop {
	my $cn    = shift;
	my $basis = shift;
	my $loc   = shift;

	my $val = "NA";

	my $pop = $loc->{"location_population_served"};
	if (defined $pop and isEmpty($pop) == 0) {
		$pop =~ s/,//i;
		$val = $cn / ($pop/$basis);
	}
	return $val;
}


sub calcCopiesFlowNorm {
	my $cpl  = shift;
	my $flow = shift;
	my $loc  = shift;

	my $val = "NA";
	if (defined $cpl and isEmpty($cpl) == 0) {
		$flow =~ s/,//i;
		my $hrs = $loc->{"location_collection_window_hrs"};
		# copies/liter x 3.78541 liters/gallon x flow million_gallons/day x 1 day/24 hrs x collection_time_hrs
		$val = $cpl * 3.78541 * $flow * ($hrs / 24);
	}
	return $val;
}


sub isEmpty {
	my $val = shift;
	
	my $is_empty = 0;
	$is_empty = 1 if "$val" eq "NaN" or "$val" eq "-" or "$val" eq "none" or "$val" eq "" or "$val" eq "TBD" or "$val" eq "NA";
	
	return $is_empty;
}

