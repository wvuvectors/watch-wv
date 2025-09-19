#!/bin/bash
source /shared/software/conda/etc/profile.d/conda.sh

#Part 2: This file runs the freyja analysis: takes basecalled data in bam format and assigns variant info. 

#example commandline usage : ./analysis_script_testing.sh RUNID path/to/rawdata barcode_file


RUNID=$1
RAWDATA=$2
FREYJABARCODES=$3


for BAMFILE in $RAWDATA/*.bam

	#get bam file and pull out barcode
	do [[ $BAMFILE =~ ([0-9]+)\.bam$ ]] 
	    #remove leading 0
	    [[ ${BASH_REMATCH[1]} =~ ^0*(.+)$ ]] 
	    #save barcode in variable $BARCODE
	    BARCODE=${BASH_REMATCH[1]}
	    echo "$BARCODE"
	    echo "$BAMFILE"
	    echo ${BASH_REMATCH[1]}

		#These variables are static for the moment, change if you want new data structure or new reference genome. 

		WORKINGDIR=$SCRATCH/COVID_ANALYSIS/${RUNID}
		CATDATA=${WORKINGDIR}/${RUNID}_${BARCODE}.fastq.gz
		REF=$SCRATCH/COVID_ANALYSIS/SCRIPTS/NC_045512_Hu-1.fasta
		BARCODE_FILE=$SCRATCH/COVID_ANALYSIS/SCRIPTS/$FREYJABARCODES
		
		mkdir ${WORKINGDIR}
		cd ${WORKINGDIR}
		
		conda activate /scratch/viv0001/pre_process_env
		#This following works for the single bam file produced by 
		
		samtools sort $BAMFILE > ${RUNID}_${BARCODE}.sorted.bam	
		samtools fastq ${RUNID}_${BARCODE}.sorted.bam > ${RUNID}_${BARCODE}.fastq	
		gzip ${RUNID}_${BARCODE}.fastq
		
		conda deactivate
		
		#Analyze full dataset with fastqc 
		
		conda activate /scratch/viv0001/fastqc_env_scratch
		
		cd ${WORKINGDIR} 
		
		fastqc ${RUNID}_${BARCODE}.fastq.gz
		
		conda deactivate
		
		#Trim file for covid data parameters--remove first 9 bases and exclude anything more than 600bp
		
		conda activate /scratch/viv0001/pre_process_env
		trimmomatic SE -phred64 ${RUNID}_${BARCODE}.fastq.gz ${RUNID}_${BARCODE}_filtered.fastq.gz HEADCROP:9 MINLEN:250
		conda deactivate 
		
		#Run fastqc on trimmed file to check results
		
		conda activate /scratch/viv0001/fastqc_env_scratch
		fastqc ${RUNID}_${BARCODE}_filtered.fastq.gz
		conda deactivate
		
		#map trimmed covid barcode data to Hu-1 reference. 
		
		conda activate /scratch/viv0001/pre_process_env
		
		cd ${WORKINGDIR}
		
		minimap2 -ax map-ont ${REF} ${RUNID}_${BARCODE}_filtered.fastq.gz > ${RUNID}_${BARCODE}_filtered.sam
		
		samtools view -b ${RUNID}_${BARCODE}_filtered.sam > ${RUNID}_${BARCODE}_filtered.bam
		samtools sort ${RUNID}_${BARCODE}_filtered.bam > ${RUNID}_${BARCODE}_filtered.sorted.bam
		samtools index ${RUNID}_${BARCODE}_filtered.sorted.bam
		
		#gather some stats on how well the reads mapped.
		
		samtools stats ${RUNID}_${BARCODE}_filtered.sorted.bam > ${RUNID}_${BARCODE}_filtered_mapping_stats.txt
		samtools idxstats ${RUNID}_${BARCODE}_filtered.sorted.bam > ${RUNID}_${BARCODE}_filtered_idxmapping_stats.txt
		
		#deactivate pre_process_env.
		
		conda deactivate
		
		#run freyja to identify covid variants in run sample. 
		
		conda activate /scratch/viv0001/freyja_env_20240718
		
		cd ${WORKINGDIR}
		
		freyja variants --variants variants_outfile_${RUNID}_${BARCODE}_filtered --depths depths_outfile_${RUNID}_${BARCODE}_filtered --ref ${REF} ${RUNID}_${BARCODE}_filtered.sorted.bam
		freyja demix variants_outfile_${RUNID}_${BARCODE}_filtered.tsv  depths_outfile_${RUNID}_${BARCODE}_filtered --output demix_${RUNID}_${BARCODE}_filtered_output --barcodes ${BARCODE_FILE}
		
		#deactivate freyja env. 
		
		conda deactivate 

 done 

# #Make end data structure and move data and analysis files to working data folders.


mkdir "0 RUN_FILES"
mkdir "1 BASECALLING"
mkdir "2 FQ_CONCAT"
mkdir "3 QC1_RESULTS"
mkdir "4 FQ_CHECKED"
mkdir "5 QC2_RESULTS"
mkdir "6 MAPPED"
mkdir "7 DEMIX"
mkdir "8 DASHBOARD"
mkdir "9 SRA"


mv *filtered.sorted.bam "${WORKINGDIR}/6 MAPPED"
mv *filtered.sorted.bam.bai "${WORKINGDIR}/6 MAPPED"
mv *_mapping_stats.txt "${WORKINGDIR}/6 MAPPED"
mv demix* "${WORKINGDIR}/7 DEMIX"
mv depths* "${WORKINGDIR}/7 DEMIX"
mv variants* "${WORKINGDIR}/7 DEMIX"


mv *_filtered.fastq.gz "${WORKINGDIR}/4 FQ_CHECKED"
mv *.fastq.gz "${WORKINGDIR}/2 FQ_CONCAT"

mv *_filtered_fastqc.html "${WORKINGDIR}/5 QC2_RESULTS"
mv *_fastqc.html "${WORKINGDIR}/3 QC1_RESULTS"
 










