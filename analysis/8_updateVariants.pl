#! /usr/bin/perl

use strict;
use warnings;

use Data::Dumper;


my @files = glob("updates/variants/*.txt");
#print Dumper(\@files);

my $watchdb_f = "updates/watchdb.LATEST.txt";

my %props = ();

foreach my $f (@files) {
	open my $fh, "<", "$f" or die "Unable to open $f: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		my ($runid, $at, $varname, $prop) = split /\t/, "$line", -1;
		next unless defined $runid and $runid ne "";
		$props{$at} = {} unless defined $props{$at};
		$props{$at}->{$varname} = $prop;
	}
	close $fh;
}

#print Dumper(\%props);
#die;

my %metadata = ();

open my $fh, "<", "$watchdb_f" or die "Unable to open $watchdb_f: $!\n";
while (my $line = <$fh>) {
	chomp $line;
	my @cols = split /\t/, "$line", -1;
	my ($at, $county, $wwtp, $start, $end, $n2_cn_per_L) = ($cols[0], $cols[5], $cols[6], $cols[11], $cols[12], $cols[26]);
	if (defined $props{$at}) {
		$metadata{$at} = {
			"at" => "$at",
			"county" => "$county",
			"start_datetime" => "$start",
			"end_datetime" => "$end",
			"facility" => "$wwtp",
			"N2_CN_per_L" => "$n2_cn_per_L"
		};
	}
}
close $fh;

#print Dumper(\%metadata);
#die;

print "AssetTag\tFacility\tCounty\tstart_datetime\tend_datetime\tVariant\tProportion\n";
foreach my $at (keys %props) {
	foreach my $varname (keys %{$props{$at}}) {
		print "$at\t$metadata{$at}->{facility}\t$metadata{$at}->{county}\t$metadata{$at}->{start_datetime}\t$metadata{$at}->{end_datetime}\t$varname\t$props{$at}->{$varname}\n";
	}
}

exit;

