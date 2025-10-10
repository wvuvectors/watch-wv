#!/bin/bash

#SBATCH -J covid_analysis
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -p standby
#SBATCH -t 4:00:00
source /shared/software/conda/etc/profile.d/conda.sh

#Did you download the most recent barcodes from the freyja-data github for covid?


#These variables must be updated for each sample run.

RAWDATA=$SCRATCH/COVID/RUN_20230822.fq/fastq_pass/barcode13/
RUNID=RUN_20230822.fq
BARCODE=13

#Trimming parameters to change

#Keep reads longer than MINLEN
MINLEN=500

#These variables are static for the moment, change if you want new data structure or new reference genome. 

WORKINGDIR=$SCRATCH/COVID_ANALYSIS/${RUNID}
CATDATA=${WORKINGDIR}/${RUNID}_${BARCODE}.fastq.gz
REF=$SCRATCH/COVID_ANALYSIS/SCRIPTS/NC_045512_Hu-1.fasta
BARCODE_FILE=$SCRATCH/COVID_ANALYSIS/SCRIPTS/usher_barcodes09_08_2023-00-48.csv

#combine raw data files in $RAWDATA and place in analysis folder ($CATDATA).

mkdir ${WORKINGDIR}
cd ${RAWDATA}
cat *.gz > ${CATDATA}

#Analyze full dataset with fastqc 

conda activate /scratch/viv0001/fastqc_env_scratch

cd ${WORKINGDIR} 

fastqc ${RUNID}_${BARCODE}.fastq.gz

conda deactivate

#Trim file for covid data parameters--remove (HEADCROP) first 9 bases and exclude less than MINLEN. 

conda activate /scratch/viv0001/pre_process_env
trimmomatic SE -phred64 ${RUNID}_${BARCODE}.fastq.gz ${RUNID}_${BARCODE}_filtered.fastq.gz HEADCROP:9 MINLEN:${MINLEN}
conda deactivate 

#Run fastqc on trimmed file to check results

conda activate /scratch/viv0001/fastqc_env_scratch
fastqc ${RUNID}_${BARCODE}_filtered.fastq.gz
conda deactivate

#Map trimmed covid barcode data to Hu-1 reference. 

conda activate /scratch/viv0001/pre_process_env

cd ${WORKINGDIR}

minimap2 -ax map-ont ${REF} ${RUNID}_${BARCODE}_filtered.fastq.gz > ${RUNID}_${BARCODE}_filtered.sam

samtools view -b ${RUNID}_${BARCODE}_filtered.sam > ${RUNID}_${BARCODE}_filtered.bam
samtools sort ${RUNID}_${BARCODE}_filtered.bam > ${RUNID}_${BARCODE}_filtered.sorted.bam
samtools index ${RUNID}_${BARCODE}_filtered.sorted.bam

#Gather some stats on how well the reads mapped.

samtools stats ${RUNID}_${BARCODE}_filtered.sorted.bam > ${RUNID}_${BARCODE}_filtered_mapping_stats.txt
samtools idxstats ${RUNID}_${BARCODE}_filtered.sorted.bam > ${RUNID}_${BARCODE}_filtered_idxmapping_stats.txt

#Deactivate pre_process_env.

conda deactivate

#Run freyja to identify variants with supplied variant file. 

conda activate /scratch/viv0001/freyja_env_scratch

cd ${WORKINGDIR}

freyja variants --variants variants_outfile_${RUNID}_${BARCODE}_filtered --depths depths_outfile_${RUNID}_${BARCODE}_filtered --ref ${REF} ${RUNID}_${BARCODE}_filtered.sorted.bam
freyja demix variants_outfile_${RUNID}_${BARCODE}_filtered.tsv  depths_outfile_${RUNID}_${BARCODE}_filtered --output demix_${RUNID}_${BARCODE}_filtered_output --barcodes ${BARCODE_FILE}

#Deactivate freyja env. 

conda deactivate 
