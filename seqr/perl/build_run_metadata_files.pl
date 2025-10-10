#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Spreadsheet::Read qw(ReadData);
use Spreadsheet::ParseXLSX; 
use DateTime::Format::Excel;

# GOALS:
# Write out run metadata for analysis 
#Input: sample sheet from Amy, and...

#hashes hard-coded for now, but would like to build each time from watch tables incase anything changes..

my %locationid2site = ( 
	"AlpineLakeWWTP-01"    => "Preston",
	"CharlestonWWTP-01"    => "Kanawha",
	"CharlesTownWWTP-01"   => "Jefferson",
	"CheatLakeWWTP-01"     => "Monongalia",
	"ClarksburgWWTP-01"    => "Harrison",
	"Coliseum-01"          => "Monongalia",
	"ElkinsWWTP-01"        => "Randolph",
	"FarmingtonWWTP-01"    => "Marion",
	"GraftonWWTP-01"       => "Taylor",
	"HuttonsvilleWWTP-01"  => "Randolph",
	"KeyserWWTP-01"        => "Mineral",
	"MadisonWWTP-01"       => "Boone",
	"ManningtonWWTP-01"    => "Marion",
	"MoundsvilleWWTP-01"   => "Marshall",
	"ParkersburgWWTP-01"   => "Wood",
	"PointPleasantWWTP-01" => "Mason",
	"PrincetonWWTP-01"     => "Mercer",
	"SalemWWTP-01"         => "Harrison",
	"StarCityWWTP-01"      => "Monongalia",
	"WestonWWTP-01"        => "Lewis",
	"WestUnionWWTP-01"     => "Doddridge",
	"WheelingWWTP-01"      => "Ohio",
	"WhiteOakWWTP-01"      => "Fayette",
	"3HG"                  => "Cabell",
	"BV"                   => "Cabell",
	"COM"                  => "Cabell",
	"CSB"                  => "Kanawha",
	"DUN"                  => "Kanawha",
	"HRSB"                 => "Putnam",
	"HSB"                  => "Cabell,Wayne",
	"Northern Wayne Co PSD"=> "Wayne",
	"Pea Ridge PSD LeSage" => "Cabell",
	"PRPSD"                => "Cabell",
	"SRM"                  => "Cabell",
	"SRO"                  => "Cabell",
	"TTE"                  => "Cabell",
	"TTW"                  => "Cabell",
	"AthensWWTP-01"        => "Mercer",
	"BrookhavenES-01"      => "Monongalia",
	"CheatLakeES-01"       => "Monongalia",
	"CollegePark-02"       => "Monongalia",
	"Dadisman-02"          => "Monongalia",
	"Domain-01"            => "Monongalia",
	"Evansdale-01"         => "Monongalia",
	"Evolution-01"         => "Monongalia",
	"Evolution-02"         => "Monongalia",
	"FollansbeeWWTP-01"    => "Brooke",
	"Hazelton-01"          => "Preston",
	"HonorsSummit-01"      => "Monongalia",
	"Lofts-01"             => "Monongalia",
	"Lofts-02"             => "Monongalia",
	"MorgantownHS-01"      => "Monongalia",
	"MountainviewES-01"    => "Monongalia",
	"MtValley-01"          => "Monongalia",
	"NorthES-01"           => "Monongalia",
	"NorthHighSt-02"       => "Monongalia",
	"OaklandE-02"          => "Monongalia",
	"PetersburgWWTP-01"    => "Grant",
	"Seneca-01"            => "Monongalia",
	"StadiumNE-01"         => "Monongalia",
	"StadiumNW-01"         => "Monongalia",
	"StadiumSE-01"         => "Monongalia",
	"StadiumSW-01"         => "Monongalia",
	"StMarysWWTP-01"       => "Pleasants",
	"SuncrestES-01"        => "Monongalia",
	"TowersBB-01"          => "Monongalia",
	"UniversityHS-01"      => "Monongalia",
	"UPlace-01"            => "Monongalia",
	"WestRun-01"           => "Monongalia",
	"BelingtonWWTP-01"     => "Barbour",
	"BlacksvilleWWTP-01"   => "Monongalia",
	"BridgeportWWTP-01"    => "Harrison",
	"CameronWWTP-01"       => "Marshall",
	"FairmontWWTP-01"      => "Marion",
	"PotomacState-01"      => "Mineral",
	"SistervilleWWTP-01"   => "Tyler",
	"WarmSpringsWWTP-01"   => "Morgan",
	"WeirtonWWTP-01"       => "Brooke",
	"WorthingtonWWTP-01"   => "Marion"
);


my @filtseq = "";

#ONE: GATHER THE PARAMETERS FOR VARIABLES

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
 
#THIS FOR LOOP WORKS, but isnt perfect. Could be a foreach loop with a counter. 

my @rows = Spreadsheet::Read::rows ($sample_data->[1]);
print Dumper @rows; 

for my $rowref (@rows) {
	my $num = $rowref->[0];
	my $asset_id = $rowref->[1];
	my $conc = $rowref->[2];
	my $barcode = $rowref->[3];
	
	#$barcode2assetid{$barcode} = $asset_id;
	
  
 	my $watchdb = "watchdb.sample.txt";
	open(my $WatchFH, "<", $watchdb);
		
	#open file for metadata out file
	
	my $RUN = "RUN_METADATA.".$RUNID.".txt";
	open(my $RUNMETAFH, ">>", $RUN); 
	
	while (my $line = <$WatchFH>) {
		#print "my asset ID is $barcode2assetid{$demixfile_bc} \n";
		chomp $line;	
			
			my $index = index($line, $asset_id);
			
			if ($index >= 0) {
				#print "$line", "\n";
				my @working = split /\t/, $line, -1;
				#print "@working", "\n"; 
				my $location_id = $working[2];
				my $collection = $working[6];
				my @collection_date = split /\s/, $collection, -1;
				#my $site = #print "$barcode2assetid{$demixfile_bc}, $location_id, $collection_date[0] \n"; 
				
				
				
########This bit prints out the total filtered reads
				my $mapping_metrics = "$RUNID" . "_" . "$barcode" . "_filtered_mapping_stats.txt";
				print "$RUNID" . "_" . "$barcode" . "_filtered_mapping_stats.txt";
				
				open(my $MappingMetricsFH, "<", $mapping_metrics); 

				while (my $mappingline = <$MappingMetricsFH>) {
						chomp $mappingline;	
			
						my $stats_index = index($mappingline, 'SN	raw total sequences:');					
						if ($stats_index >= 0) {
							my @filtseq = split /\t/, $mappingline, -1;
							print "$filtseq[2]  \n" ;
							
########This bit prints out the mapped reads						
						while (my $mappingline = <$MappingMetricsFH>) {
						chomp $mappingline;	

						
							my $mapped_index = index($mappingline, 'SN	reads mapped:');					
							if ($mapped_index >= 0) {
								my @mappedseq = split /\t/, $mappingline, -1;
								print "$mappedseq[2] \n" ;
								
########This bit prints out the demix percentage						
								my $demix = "demix_" . "$RUNID" . "_" . "$barcode" . "_filtered_output";
								print "demix_" . "$RUNID" . "_" . "$barcode" . "_filtered_output \n";
								
								open(my $DemixFH, "<", $demix); 

									while (my $demixline = <$DemixFH>) {
									chomp $demixline;	
			
											my $demix_index = index($demixline, 'coverage');					
											if ($demix_index >= 0) {
											my @demix_info = split /\t/, $demixline, -1;
											print "$demix_info[1]  \n" ;
			
	
########This bit prints out the total reads ....in flux						
								my $htmlfile = "$RUNID" . "_" . "$barcode" . "_fastqc.html";
								print "$RUNID" . "_" . "$barcode" . "_fastqc.html \n";
								
								open(my $htmlFH, "<", $htmlfile); 

									while (my $htmlline = <$htmlFH>) {
									chomp $htmlline;	
			
											my $html_index = index($htmlline, 'Total Sequences');					
											if ($html_index >= 0) {
											my @html_info = split 'Total Sequences</td><td>', $htmlline, -1;
											#print "$html_info[1]  \n" ;
											my @html_reads = split '</td>', $html_info[1], -1;
											print "$html_reads[0]  \n" ;
	
												
########This bit prints out the metadata file!!!!!		
				print $RUNMETAFH "$RUNID";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$barcode";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$asset_id";
				print $RUNMETAFH "\t";				
				print $RUNMETAFH "$location_id";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$locationid2site{$location_id}";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$collection_date[0]";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "Date_analyzed";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$PRIMERS";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$FLOWCELL";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "basecalling_program";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "basecalling_model";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "freyja_barcode";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "filtering_parameters";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$html_reads[0]";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$filtseq[2]";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "PERCENT READS FILTERED";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$mappedseq[2]";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "PERCENT MAPPED READS";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "$demix_info[1]";
				print $RUNMETAFH "\t";
				print $RUNMETAFH "\n";
				

                                                        }
                                                  }
					                           }
						                	}			
							
                                          }
                                       }
                                     }
                                }
			
			}	

     }
}