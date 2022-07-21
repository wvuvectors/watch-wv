#! /usr/bin/env perl

use strict;
use warnings;

my @hdr = ();
my %indata = ();
my %outdata = ();

my $count = 0;
while (my $line = <>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @values = split "\t", "$line", -1;

	if ($count == 0) {
		for (my $i=0; $i<scalar(@values); $i++) {
			push @hdr, "$values[$i]";
		}
	} else {
		$indata{"$values[0]"} = {};
		for (my $i=1; $i<scalar(@values); $i++) {
			$indata{"$values[0]"}->{"$hdr[$i]"} = "$values[$i]";
		}
		my $subcount = 1;
		my $sid = "$values[0]";
		$sid =~ s/\.a1$//;
		for (my $i=1; $i <= 4; $i++) {
			if (defined $indata{"$values[0]"}->{"subject_$i"} and $indata{"$values[0]"}->{"subject_$i"} ne "") {
				$outdata{"$sid.a$subcount"} = {};
				$outdata{"$sid.a$subcount"}->{"sample_id"} = "$sid";
				$outdata{"$sid.a$subcount"}->{"assay_id"} = "$subcount";
				$outdata{"$sid.a$subcount"}->{"extract_id"} = $indata{"$values[0]"}->{"extract_id"};
				$outdata{"$sid.a$subcount"}->{"plate_id"} = $indata{"$values[0]"}->{"plate_id"};
				$outdata{"$sid.a$subcount"}->{"plate_cell"} = $indata{"$values[0]"}->{"plate_cell"};
				$outdata{"$sid.a$subcount"}->{"assay_input_ml"} = $indata{"$values[0]"}->{"assay_input_ml"};
				$outdata{"$sid.a$subcount"}->{"assay_type"} = $indata{"$values[0]"}->{"assay_type"};
				$outdata{"$sid.a$subcount"}->{"accepted_droplets"} = $indata{"$values[0]"}->{"accepted_droplets"};
				$outdata{"$sid.a$subcount"}->{"notes"} = $indata{"$values[0]"}->{"notes"};

				$outdata{"$sid.a$subcount"}->{"target"} = $indata{"$values[0]"}->{"subject_${i}"};
				$outdata{"$sid.a$subcount"}->{"target_category"} = $indata{"$values[0]"}->{"subject_${i}_category"};
				$outdata{"$sid.a$subcount"}->{"target_gene"} = $indata{"$values[0]"}->{"subject_${i}_target"};
				$outdata{"$sid.a$subcount"}->{"target_macromolecule"} = $indata{"$values[0]"}->{"subject_${i}_target_macromolecule"};
				$outdata{"$sid.a$subcount"}->{"target_fluorophore"} = $indata{"$values[0]"}->{"subject_${i}_fluorophore"};
				$outdata{"$sid.a$subcount"}->{"target_calculated_copies_per_reaction"} = $indata{"$values[0]"}->{"subject_${i}_calculated_copies_per_reaction"};
				$outdata{"$sid.a$subcount"}->{"target_copies_per_reaction"} = $indata{"$values[0]"}->{"subject_${i}_copies_per_reaction"};
				$outdata{"$sid.a$subcount"}->{"target_result"} = $indata{"$values[0]"}->{"subject_${i}_result"};
				
				$subcount++;
			}
		}
	}

	$count++;
}

my @outhdr = ("sample_id", "assay_id", "extract_id","plate_id","plate_cell","assay_input_ml","assay_type","accepted_droplets","target","target_category","target_gene","target_macromolecule","target_fluorophore","target_calculated_copies_per_reaction","target_copies_per_reaction","target_result","notes");

print "run_id\t" . join("\t", @outhdr) . "\n";

foreach my $key (keys %outdata) {
	print "$key";
	foreach my $h (@outhdr) {
		if (defined $outdata{$key}->{$h}) {
			print "\t$outdata{$key}->{$h}";
		} else {
			print "\t";
		}
	}
	print "\n";
}


