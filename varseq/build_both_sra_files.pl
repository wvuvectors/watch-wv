#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX; 
use DateTime::Format::Excel;


# GOALS:
# Write out biosample and metadata files for SRAFH upload 
#Input: sample sheet from Amy, and....

my %locationid2NWSSid = ( 
	"AlpineLakeWWTP-01"    => "54077-001-01-00-00",
	"CharlestonWWTP-01"    => "54039-002-01-00-00",
	"CharlesTownWWTP-01"   => "54037-001-01-00-00",
	"CheatLakeWWTP-01"     => "54061-001-01-00-00",
	"ClarksburgWWTP-01"    => "54033-002-01-00-00",
	"ElkinsWWTP-01"        => "54083-003-01-00-00",
	"FarmingtonWWTP-01"    => "54049-002-01-00-00",
	"GraftonWWTP-01"       => "54091-001-01-00-00",
	"HuttonsvilleWWTP-01"  => "54083-002-01-00-00",
	"KeyserWWTP-01"        => "54057-001-01-00-00",
	"MadisonWWTP-01"       => "54005-001-01-00-00",
	"ManningtonWWTP-01"    => "54049-001-01-00-00",
	"MoundsvilleWWTP-01"   => "54051-003-01-00-00",
	"ParkersburgWWTP-01"   => "54107-001-01-00-00",
	"PointPleasantWWTP-01" => "54053-001-01-00-00",
	"PrincetonWWTP-01"     => "54055-002-01-00-00",
	"SalemWWTP-01"         => "54033-003-01-00-00",
	"StarCityWWTP-01"      => "54061-002-01-00-00",
	"WestonWWTP-01"        => "54041-001-01-00-00",
	"WestUnionWWTP-01"     => "54017-001-01-00-00",
	"WheelingWWTP-01"      => "54051-004-01-00-00",
	"WhiteOakWWTP-01"      => "54019-001-01-00-00"
);

my %locationid2popserved = ( 
	"AlpineLakeWWTP-01"    => "700",
	"CharlestonWWTP-01"    => "50000",
	"CharlesTownWWTP-01"   => "17000",
	"CheatLakeWWTP-01"     => "2000",
	"ClarksburgWWTP-01"    => "26498",
	"ElkinsWWTP-01"        => "13156",
	"FarmingtonWWTP-01"    => "610",
	"GraftonWWTP-01"       => "6071",
	"HuttonsvilleWWTP-01"  => "2101",
	"KeyserWWTP-01"        => "8168",
	"MadisonWWTP-01"       => "4555",
	"ManningtonWWTP-01"    => "1091",
	"MoundsvilleWWTP-01"   => "12000",
	"ParkersburgWWTP-01"   => "48050",
	"PointPleasantWWTP-01" => "5515",
	"PrincetonWWTP-01"     => "36000",
	"SalemWWTP-01"         => "1853",
	"StarCityWWTP-01"      => "48328",
	"WestonWWTP-01"        => "10364",
	"WestUnionWWTP-01"     => "564",
	"WheelingWWTP-01"      => "50000",
	"WhiteOakWWTP-01"      => "2626"
);

my %primerset2ampliconsize = ( 
	"ARTIC V3"            => "400",
	"Artic v3"            => "400",
	"VarSkip Short V2"    => "550"
);

#THESE LIKELY DO NOT CHANGE
my $sample_file = glob 'Samples*';
my %barcode2assetid = ();
my %assetid2file = ();

#hard-coded variables for biosamples sheet:
my $bioproject_accession = "PRJNA1051042";
my $organism = "wastewater metagenome";
my $geo_loc_name = "USA: West Virginia";
my $isolation_source = "Wastewater";
my $project_name = "CDC NWSS";
my $collected_by = "Not Applicable";
my $purpose_of_ww_sampling = "public health surveillance community-level";
my $ww_sample_site = "wastewater treatment plant"; 
my $ww_sample_matrix = "raw wastewater"; 
my $ww_sample_type = "composite";
my $ww_sample_duration = "24";
my $concentration_method = "ceres nanotrap";
my $extraction_method = "thermofisher magmax viral/pathogen nucleic acid isolation kit";
my $extraction_control = "Not Applicable";
my $ww_surv_target_1 = "SARS-CoV-2"; 
my $ww_surv_target_1_known_present = "Yes";

#hard-coded cariables for the (sequencing) metadata sheet:
my $title = "Wastewater metagenomic sequencing";
my $library_strategy = "AMPLICON";
my $library_source = "METAGENOMIC";
my $library_selection = "PCR"; 
my $library_layout = "single";
my $platform = "OXFORD_NANPOPORE";
my $design_description = "Methods included in custom attributes";
my $filetype = "fastq";
my $enrichment_kit = "NEBNext ARTIC SARS-CoV-2 RT-PCR Module";
my $library_preparation_kit = "NEBNext ARTIC SARS-CoV-2 Companion Kit";
my $quality_control_determination = "no quality control issues identified";
my $sequence_submitter_contact_email = 'timothy.driscoll@mail.wvu.edu';
my $raw_sequence_data_processing_method  = "Raw reads basecalled and barcodes removed by Dorado";



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

#OPEN OUTPUT FILES AND ADD HEADINGS

#open SRA file for biosamples out file and print headings
my $SRA = "SRA_BIOSAMPLES_".$RUNID.".txt";
open(my $SRAFH, ">>", $SRA);
    print $SRAFH "sample_name";
    print $SRAFH "\t";
    print $SRAFH "bioproject_accession";
    print $SRAFH "\t";
    print $SRAFH "organism";
    print $SRAFH "\t";
    print $SRAFH "collection_date";
    print $SRAFH "\t";
    print $SRAFH "geo_loc_name";
    print $SRAFH "\t";
    print $SRAFH "isolation_source";
    print $SRAFH "\t";
    print $SRAFH "collection_site_id";
    print $SRAFH "\t";
    print $SRAFH "project_name";
    print $SRAFH "\t";
    print $SRAFH "collected_by";
    print $SRAFH "\t";
    print $SRAFH "purpose_of_ww_sampling";
    print $SRAFH "\t";
    print $SRAFH "ww_sample_site";
    print $SRAFH "\t";
    print $SRAFH "ww_population";
    print $SRAFH "\t";
    print $SRAFH "ww_sample_matrix";
    print $SRAFH "\t";
    print $SRAFH "ww_sample_type";
    print $SRAFH "\t";
    print $SRAFH "ww_sample_duration";
    print $SRAFH "\t";
    print $SRAFH "concentration_method";
    print $SRAFH "\t";
    print $SRAFH "extraction_method";
    print $SRAFH "\t";
    print $SRAFH "extraction_control";
    print $SRAFH "\t";
    print $SRAFH "ww_surv_target_1";
    print $SRAFH "\t";
    print $SRAFH "ww_surv_target_1_known_present";
    print $SRAFH "\n";
    
#open SRA file for metadata out file and print headings
my $META = "SRA_METADATA_".$RUNID.".txt";
open(my $METAFH, ">>", $META);
    print $METAFH "sample_name";
    print $METAFH "\t";
    print $METAFH "library_ID";
    print $METAFH "\t";
    print $METAFH "title";
    print $METAFH "\t";
    print $METAFH "library_strategy";
    print $METAFH "\t";
    print $METAFH "library_source";
    print $METAFH "\t";
    print $METAFH "library_selection";
    print $METAFH "\t";
    print $METAFH "library_layout";
    print $METAFH "\t";
    print $METAFH "platform";
    print $METAFH "\t";
    print $METAFH "instrument_model";
    print $METAFH "\t";
    print $METAFH "design_description";
    print $METAFH "\t";
    print $METAFH "filetype";
    print $METAFH "\t";
    print $METAFH "filename";
    print $METAFH "\t";
    print $METAFH "enrichment_kit";
    print $METAFH "\t";
    print $METAFH "amplicon_PCR_primer_scheme";
    print $METAFH "\t";
    print $METAFH "library_preparation_kit";
    print $METAFH "\t";
    print $METAFH "amplicon_size";
    print $METAFH "\t";
    print $METAFH "quality_control_determination";
    print $METAFH "\t";
    print $METAFH "sequence_submitter_contact_email";
    print $METAFH "\t";
    print $METAFH "raw_sequence_data_processing_method";
    print $METAFH "\n";


 
#THIS FOR LOOP WORKS, but isnt perfect. Could be a foreach loop with a counter. 

my @rows = Spreadsheet::Read::rows ($sample_data->[1]);
print Dumper @rows; 

for my $rowref (@rows) {
	my $num = $rowref->[0];
	my $asset_id = $rowref->[1];
	my $conc = $rowref->[2];
	my $barcode = $rowref->[3];
	
	$barcode2assetid{$barcode} = $asset_id;
	}
	
#SECOND, read the bam files and associate each with its asset id.
#my @files = readdir(./*.bam);
my @files = glob '*.bam';
foreach my $bamfiles (@files) {
	if ($bamfiles =~ /.+(\d{2})\.bam/) {
		# extract the barcode from the file name.
		my $file_bc = $1;
		# remove the leading zero, if nec.
		$file_bc =~ s/^0//;
		# look up the barcode to get the asset id.
		my $assetid = $barcode2assetid{$file_bc};
		# add to the hash mapping asset to file name.
		$assetid2file{$assetid} = $bamfiles;
			
	}
}
 
my @demix = glob "demix*";
foreach my $demixfiles (@demix) {
		
		my @lineages = ();
		if ($demixfiles =~ /(\d+)\_filtered_output/) {
		
		# extract the barcode from the file name.
		print "$demixfiles \n";
		my $demixfile_bc = $1;
		# remove the leading zero, if nec.
		$demixfile_bc =~ s/^0//;


#my $watchdb = "watchdb.sample.txt";
my $watchdb = "/Users/viv0001/github/watch-wv/dashboard/data/watchdb.sample.txt";

open(my $WatchFH, "<", $watchdb);

		#READ IN EACH LINE of WATCHDB SAMPLE SHEET AND IF THE CURRENT ASSETID SAVE THE INFO
		
		while (my $line = <$WatchFH>) {
		#print "my asset ID is $barcode2assetid{$demixfile_bc} \n";
		chomp $line;	
				
			my $index = index($line, $barcode2assetid{$demixfile_bc});
			
			if ($index >= 0) {
				print "$line", "\n";
				my @working = split /\t/, $line, -1;
				#print "@working", "\n"; 
				my $location_id = $working[2];
				my $collection = $working[6];
				my @collection_date = split /\s/, $collection, -1;
				print "$barcode2assetid{$demixfile_bc}, $location_id, $collection_date[0] \n"; 
				
				print $SRAFH "$barcode2assetid{$demixfile_bc}";
				print $SRAFH "\t";
				print $SRAFH "$bioproject_accession";
				print $SRAFH "\t";
				print $SRAFH "$organism";
				print $SRAFH "\t";
				print $SRAFH "$collection_date[0]";
				#collection date corresponds to when the collection ends in the 24 hour period
				#ie if collection begins on the 6th at 9AM and ends on the 7th at 9AM the date here is the 7th
				print $SRAFH "\t";
				print $SRAFH "$geo_loc_name";
				print $SRAFH "\t";
				print $SRAFH "$isolation_source";
				print $SRAFH "\t";
				print $SRAFH "$locationid2NWSSid{$location_id}";
				print $SRAFH "\t";
				print $SRAFH "$project_name";
				print $SRAFH "\t";
				print $SRAFH "$collected_by";
				print $SRAFH "\t";
				print $SRAFH "$purpose_of_ww_sampling";
				print $SRAFH "\t";
				print $SRAFH "$ww_sample_site";
				print $SRAFH "\t";
				print $SRAFH "$locationid2popserved{$location_id}";
				print $SRAFH "\t";
				print $SRAFH "$ww_sample_matrix";
				print $SRAFH "\t";
				print $SRAFH "$ww_sample_type";
				print $SRAFH "\t";
				print $SRAFH "$ww_sample_duration";
				print $SRAFH "\t";
				print $SRAFH "$concentration_method";
				print $SRAFH "\t";
				print $SRAFH "$extraction_method";
				print $SRAFH "\t";
				print $SRAFH "$extraction_control";
				print $SRAFH "\t";
				print $SRAFH "$ww_surv_target_1";
				print $SRAFH "\t";
				print $SRAFH "$ww_surv_target_1_known_present";
				print $SRAFH "\n";
				
				print $METAFH "$barcode2assetid{$demixfile_bc}";
				print $METAFH "\t";
				print $METAFH "$barcode2assetid{$demixfile_bc}"."."."$RUNID";
				print $METAFH "\t";
				print $METAFH "$title";
				print $METAFH "\t";
				print $METAFH "$library_strategy";
				print $METAFH "\t";
				print $METAFH "$library_source";
				print $METAFH "\t";
				print $METAFH "$library_selection";
				print $METAFH "\t";
				print $METAFH "$library_layout";
				print $METAFH "\t";
				print $METAFH "$platform";
				print $METAFH "\t";
				print $METAFH "$INSTRUMENT";
				print $METAFH "\t";
				print $METAFH "$design_description";
				print $METAFH "\t";
				print $METAFH "$filetype";
				print $METAFH "\t";
				print $METAFH "$RUNID"."_"."$demixfile_bc".".fastq.gz";
				print $METAFH "\t";
				print $METAFH "$enrichment_kit";
				print $METAFH "\t";
				print $METAFH "$PRIMERS";
				print $METAFH "\t";
				print $METAFH "$library_preparation_kit";
				print $METAFH "\t";
				print $METAFH "$primerset2ampliconsize{$PRIMERS}";
				print $METAFH "\t";
				print $METAFH "$quality_control_determination";
				print $METAFH "\t";
				print $METAFH "$sequence_submitter_contact_email";
				print $METAFH "\t";
				print $METAFH "$raw_sequence_data_processing_method";
				print $METAFH "\n";
		
		
		
			}
		}	 
	}
}
