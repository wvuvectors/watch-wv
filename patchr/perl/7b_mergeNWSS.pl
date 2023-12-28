#! /usr/bin/env perl

use strict;
use warnings;

use Text::CSV qw(csv);
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage   .= "Usage: $progname FLIST\n";
$usage   .=  "Merges all of the files in FLIST and writes merged file to STDOUT.\n";
$usage   .=  "The first line of each file is considered a header. If a file is passed on STDIN, ";
$usage   .=  "it is used to determine the output column headers and order. Otherwise, OUTER JOIN rules apply.\n";
$usage   .=   "\n";


my $status = 0;

my @flist  = ();

while (@ARGV) {
  my $arg = shift;
  if ("$arg" eq "-h") {
		die "$usage";
	} else {
		push @flist, "$arg";
	}
}

die "No files were provided so no merge will take place.\n$usage\n" unless scalar @flist > 0;
die "Only one file was provided so no merge will take place.\n$usage\n" unless scalar @flist > 1;

my %data = ();
my @outheaders     = ();
my $headers_preset = 0;

# Read in field names and order from file, if one is provided.
while (my $line = <>) {
	chomp $line;
	push @outheaders, trim("$line") unless $line =~ /^\s*$/;
}
$headers_preset = 1 if scalar @outheaders > 0;

my $filenum = 1;
foreach my $f (@flist) {
	my $csvA  = Text::CSV->new({auto_diag => 4, binary => 1});
	my $linenum = 0;
	my @colheaders = ();
	open (my $fh, "<", "$f") or die "Unable to open $f for reading: $!\n";
	while (my $line = $csvA->getline($fh)) {
		my @cols = @$line;
		
		# Set the global headers if they are not preset.
		if (scalar @outheaders == 0) {
			foreach (@cols) {
				push @outheaders, trim("$_");
			}
		}
		
		# Set the file colheaders if this is the first line.
		# These are used as keys of the data sub-hashes.
		if ($linenum == 0) {
			foreach (@cols) {
				push @colheaders, trim("$_");
			}
		} else {
			my $uid = "f" . "$filenum" . "-" . "row" . "$linenum";
			for (my $i=0; $i < scalar(@cols); $i++) {
				my $val = $cols[$i];
				$data{"$uid"} = {} unless defined $data{"$uid"};
				$data{"$uid"}->{"$colheaders[$i]"} = "$val";
			}
		}
		$linenum++;
	}
	close $fh;
	$filenum++;
}

#print Dumper(\%data);
#die;

print join("\",\"", @outheaders) . "\"\n";
foreach my $uid (keys %data) {
	my $start = 0;
	foreach my $hdr (@outheaders) {
		my $val = "";
		$val = "$data{$uid}->{$hdr}" if defined $data{"$uid"}->{"$hdr"};
		print "," unless $start == 0;
		print "\"" . "$val" . "\"";
		$start++;
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


