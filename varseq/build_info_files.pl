#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#declare yo variables

# GOALS:
# Write out .txt files with info for use later on! 
#Input: samples_date.txt on stdin


# First build the map of barcode-to-asset.
my %barcode2assetid = ();
my $line ;
my $RUNID = "SARS_20240516";
my %assetid2file = ();


# For short lines, I do this:

while (my $line = (<>)) {
	chomp $line;
	#my @a = split(/\t/, "$line", -1);
	my ($num, $asset_id, $conc, $barcode) = split(/\t/, "$line", -1);

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
		
		my $outputname = "SARS_20240516_".$demixfile_bc."_".$barcode2assetid{$demixfile_bc}.".txt";
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





