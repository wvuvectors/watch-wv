#! /usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

#print "******\n";
#print "Running $progname.\n";
#print "******\n";


my $usage = "\n";
$usage   .= "Usage: $progname [options] DBDIR\n";
$usage   .=  "Extracts all the batch ids from the tables in DBDIR and writes to STDOUT.\n";
$usage   .=   "\n";

my $dbdir;

while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h") {
		die $usage;
	} else {
		$dbdir = $arg;
	}
}

my $status = 0;

die "FATAL: $progname requires a valid database directory ($dbdir).\n$usage\n" unless defined $dbdir and -d "$dbdir";

my %idcols = (
	"abatch" => 0,
	"cbatch" => 0,
	"ebatch" => 0,
	"rbatch" => 0
);


opendir(my $dbdFH, "$dbdir") or die "Although I found it, I am unable to open $dbdir: $!";
my @dbfiles  = grep { (/^watchdb\..+\.txt$/) && (!/^~/) && -f "$dbdir/$_" } readdir($dbdFH);
for my $f (@dbfiles) {
	my $table = "$f";
	$table =~ s/watchdb\.(.+?)\.txt/$1/gi;
	next unless defined $table and defined $idcols{"$table"};

	my $linenum = 0;
	my $idcol = $idcols{"$table"};

	open(my $fh, "<", "$dbdir/$f") or die "FATAL: Unable to open $dbdir/$f for reading: $!\n";
	while (my $line = <$fh>) {
		# Skip the header line.
		if ($linenum == 0) {
			$linenum++;
			next;
		}
		chomp $line;
		my @cols = split /\t/, "$line", -1;
		if (defined $cols[$idcol] and "$cols[$idcol]" ne "") {
			print "$cols[$idcol]\n";
		}
	}
	close $fh;
}

exit $status;
