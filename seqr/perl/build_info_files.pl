#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use strict;
use Data::Dumper;
use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX; 
use DateTime::Format::Excel;

#declare yo variables

# GOALS:
# Write out .txt files with info for use later on! 
#Input: build_info_files.pl <samples_date.txt


#Pull in RUNID from excel
my $sample_file = glob 'Samples*';


# READ IN EXCEL SAMPLE FILE
my $sample_data = ReadData($sample_file, dtfmt => "mm/dd/yy");

#Save runinfo for later data on page 1 (pages start with 1), row 16 (rows start with 1), column 1 (columns start with 0)	
my @sheet  = Spreadsheet::Read::rows($sample_data->[1]);
my @row_1 = Spreadsheet::Read::row($sample_data->[1], 1);
my $RUNID  = $row_1[1];
my @row_2 = Spreadsheet::Read::row($sample_data->[1], 2);
my $FLOWCELL  = $row_2[1];
my @row_3 = Spreadsheet::Read::row($sample_data->[1], 3);
my $PRIMERS  = $row_3[1];
my @row_4 = Spreadsheet::Read::row($sample_data->[1], 4);
my $INSTRUMENT  = $row_4[1];


 print Dumper ($RUNID);
 print Dumper ($FLOWCELL);
 print Dumper ($PRIMERS);
 print Dumper ($INSTRUMENT);
 


# First build the map of barcode-to-asset.
my %barcode2assetid = ();
my $line ;
my %assetid2file = ();


my @rows = Spreadsheet::Read::rows ($sample_data->[1]);
print Dumper @rows; 

for my $rowref (@rows) {
	my $num = $rowref->[0];
	my $asset_id = $rowref->[1];
	my $conc = $rowref->[2];
	my $barcode = $rowref->[3];
	
	$barcode2assetid{$barcode} = $asset_id;
	
	}


print Dumper(\%barcode2assetid);


open FH, '>', "barcode2assetid.txt";
print FH Dumper(\%barcode2assetid);


# Second, read the bam files and associate each with its asset id.


#my @files = readdir(./*.bam);
my @files = glob '*.bam';
foreach my $f (@files) {
	if ($f =~ /.+(\d{2})\.bam/) {
		# extract the barcode from the file name.
		my $file_bc = $1;
		# remove the leading zero, if nec.
		$file_bc =~ s/^0//;
		# look up the barcode to get the asset id.
		my $assetid = $barcode2assetid{$file_bc};
		# add to the hash mapping asset to file name.
		$assetid2file{$assetid} = $f;
			
	}
}

print Dumper(\%assetid2file); 

open FH, '>', "assetid2file.txt";
print FH Dumper(\%assetid2file);




#open a demix file
my @demix = glob "demix*";
foreach my $f (@demix) {
		
		my @lineages = ();
		if ($f =~ /(\d+)\_filtered_output/) {
		
		# extract the barcode from the file name.
		print "$f \n";
		my $demixfile_bc = $1;
		# remove the leading zero, if nec.
		$demixfile_bc =~ s/^0//;
		#print  "$barcode2assetid{$demixfile_bc} \n";
		
		my $outputname = $RUNID."_".$demixfile_bc."_".$barcode2assetid{$demixfile_bc}.".txt";
		print "$outputname \n"; 
		


		open(my $FH, "<", $f);
		open(my $OUTFH, ">", $outputname );


		while (my $line = <$FH>) {
		chomp $line;	
	
		if ( $line =~ /^lineages/ ){
			my @arr1 = split /\t/, "$line", -1;
			my @arr2 = split /\s/, $arr1[1], -1;	
			foreach (@arr2) {
	   			push @lineages, "$_";
			}
		
	#	print Dumper (\@lineages); die; 
		
		
		} elsif ($line =~ /^abundances/ ) {
			my @placeholder1 = split /\t/, "$line", -1;
			my @placeholder2 = split /\s/, $placeholder1[1], -1;
			for (my $i=0; $i<scalar(@placeholder2); $i++) {
				print $OUTFH "$RUNID";
				print $OUTFH "\t";
				print $OUTFH "$barcode2assetid{$demixfile_bc}";   
				print $OUTFH "\t";
				print $OUTFH "$lineages[$i]";
				print $OUTFH "\t";
				print $OUTFH "$placeholder2[$i]"; 
				print $OUTFH "\n";
			}
		}			
		


		}
 close ($FH);
 close ($OUTFH); 
 
}	
   
}





