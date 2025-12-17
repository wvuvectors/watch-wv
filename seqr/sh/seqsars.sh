#!/bin/bash

# SEQSARS guarded execution verion
#
# This script runs freyja analysis on barcoded COVID sequencing data from the ONT platform.
# It is generally called from an HPC scheduler script, but it can also be run separately.
# If you want to run this on your laptop, be warned that the processing requirements are not trivial.
#

#
# Fail immediately on:
# - unset variables
# - non-zero exit codes
# - failures inside pipes
set -Eeuo pipefail

# Better glob handling
shopt -s nullglob
shopt -s failglob

# Helpful error reporting
trap 'echo "[FATAL] Error on line $LINENO: $BASH_COMMAND" >&2' ERR

# Logging 
log() {
    echo "[$(date '+%F %T')] $*"
}

fail() {
    echo "[FATAL] $*" >&2
    exit 1
}

START=$(date "+%F_%H-%M")
log "Starting SEQSARS run at $START"

while getopts ":hi:r:b:g:o:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 0
			;;
		i)
			runDIR=$OPTARG
			;;
		r)
			runID=$OPTARG
			;;
		b)
			usher_barcodesFILE=$OPTARG
			;;
		g)
			refGenomeFILE=$OPTARG
			;;
		o)
			outDIR=$OPTARG
			;;
# 		c)
# 			condaDIR=$OPTARG
# 			;;
		\?)
			fail "Invalid option: -$OPTARG"
			;;
		:)
			fail "Option -$OPTARG requires an argument."
			;;
	esac
done


# Required parameters that have NO default value.
#
if [ -z "$runDIR" ]; then
	echo "Missing a valid input run directory (-i)!"
	echo "This is the parent run folder that contains the \'fastq_pass\' folder (among others)."
	echo "This is a required parameter with no default, so I'm forced to quit immediately. Sorry."
	exit 1
fi


# Required parameters for which I can try to provide a default value.
#

if [ -z "$usher_barcodesFILE" ]; then
	usher_barcodesFILE="resources/usher_barcodes.csv"
fi
if [ ! -f "$usher_barcodesFILE" ]; then
	echo "A valid USHER barcodes file for freyja was not provided (-b), so I tried the default $usher_barcodesFILE."
	echo "Sadly, $usher_barcodesFILE either does not exist or cannot be read."
	echo "The usher barcodes file is usually located in the project resources folder, and is named something akin to \'usher_barcodes.csv\'."
	echo "You can try to download it directly from https://github.com/andersen-lab/Freyja/blob/main/freyja/data/usher_barcodes.csv."
	echo "This is a required parameter, so I'm forced to quit immediately. Sorry."
	exit 1
fi

if [ -z "$refGenomeFILE" ]; then
	refGenomeFILE="resources/NC_045512_Hu-1.fasta"
fi
if [ ! -f "$refGenomeFILE" ]; then
	echo "A valid reference genome file was not provided (-g), so I tried the default $refGenomeFILE."
	echo "Sadly, $refGenomeFILE either does not exist or cannot be read."
	echo "The reference genome file is usually located in the project resources folder, and is named something akin to \'NC_045512_Hu-1.fasta\'."
	echo "This is a required parameter, so I'm forced to quit immediately. Sorry."
	exit 1
fi

if [ -z "$runID" ]; then
	echo "Missing a valid run ID (-r)!"
	echo "I'm going to use a default \'RUNON_$START\' but be aware:"
	echo "this run ID contains the date the analysis was started, not the date of sequencing."
	runID="RUNON_$START"
fi

if [ -z "$outDIR" ]; then
	echo "Missing a valid output directory (-o)!"
	echo "Writing output to the default dir of \'output/$runID/\'"
	mkdir "output"
	outDIR="output/$runID"
fi

# if [ -z "$condaDIR" ]; then
# 	echo "Missing a valid conda environment directory (-c)!"
# 	echo "I'll try to use the current directory, but if you get a conda environment error that might be why."
# 	condaDIR="."
# fi


if [[ "$outDIR" == "output/$runID" ]]; then
	echo ""
	echo "WARNING! If $outDIR is \'output/$runID\' don't forget to MOVE IT out of the github repository ASAP."
	echo "You MUST move this dir out of the git repository BEFORE committing! MUST!"
	echo "If you fail to do so, the next commit and push to github will FAIL (it's too big)."
	echo ""
	echo "You have been warned."
	echo ""
fi

# Checking tool availability
#log "Checking required tools in PATH"
#for tool in minimap2 samtools fastqc trimmomatic freyja; do
#    command -v "$tool" >/dev/null || fail "Required tool not found in PATH: $tool"
#done 

# Create the output directory if it doesn't already exist.
log "Creating output directory structure in $outDIR"
if [ ! -d "$outDIR" ]; then
	mkdir "$outDIR"
fi

mkdir -p "$outDIR/0 RUN_FILES"
mkdir -p "$outDIR/1 BASECALLING"
mkdir -p "$outDIR/2 FQ_CONCAT"
mkdir -p "$outDIR/3 QC1_RESULTS"
mkdir -p "$outDIR/4 FQ_CHECKED"
mkdir -p "$outDIR/5 QC2_RESULTS"
mkdir -p "$outDIR/6 MAPPED"
mkdir -p "$outDIR/7 DEMIX"
mkdir -p "$outDIR/8 DASHBOARD"
mkdir -p "$outDIR/9 SRA"

# Barcode discovery
log "Discovering barcode directories"
barcode_dirs=( "$runDIR"/fastq_pass/barcode* )
(( ${#barcode_dirs[@]} > 0 )) || fail "No barcode directories found in $runDIR/fastq_pass"

# Main loop
for bcDir in "${barcode_dirs[@]}"; do
	[[ -d "$bcDir" ]] || continue 
	
	log "Processing barcode directory: $bcDir" 
	bcBase=$(basename "$bcDir")
	bcBase="${bcBase//$'\r'/}"  # remove any carriage return characters
	log " Found $bcBase"
	
	# Extract the barcode id to a variable $barcode. 
	if [[ $bcBase =~ barcode([0-9]+) ]]; then
		barcode=${BASH_REMATCH[1]}
		barcode2d="$barcode"
	else 
		fail "Unexpected barcode directory name: $bcBase"
	fi 
	
	# Deal with a leading zero on the barcode number.
	if [[ $barcode =~ ^0([0-9]+)$ ]]; then
		barcode=${BASH_REMATCH[1]}
	fi
	
	# Paths for concatenated and trimmed FASTQ
	concat_fastq="$outDIR/2 FQ_CONCAT/$bcBase.fastq.gz"
	trimmed_fastq="$outDIR/4 FQ_CHECKED/${runID}_${bcBase}_tr.fastq.gz"
	
	# Check that FASTQ files exist before concatenation
	fastq_files=( "$bcDir"/*.fastq.gz )
	(( ${#fastq_files[@]} > 0 )) || fail "No FASTQ files found in $bcDir"

	# Concatenate all the fastq.gz files into one file.
	#
	log "Concatenating FASTQs for barcode $barcode2d"
	ls "$runDIR/fastq_pass/$bcDir/*.fastq.gz" | xargs cat > "$outDIR/2 FQ_CONCAT/$bcDir.fastq.gz"


	# Run fastqc on the original file.
	#
	log "Running FastQC on raw files"
	conda activate "seqr_fastqc" || fail "Failed to activate conda env: seqr_fastqc"
	fastqc -o "$outDIR/3 QC1_RESULTS/" "$outDIR/2 FQ_CONCAT/$bcDir.fastq.gz"
	conda deactivate

		
	# Trim the reads in the concatenated fastq file. These are specific for COVID sequencing:
	# Remove the first 9 bases, and exclude any reads with length > 600nt.
	#
	log "Trimming reads"
	conda activate "seqr_general" || fail "Failed to activate conda env: seqr_general"
	trimmomatic SE -phred64 "$outDIR/2 FQ_CONCAT/${runID}_$bcDir.fastq.gz" "$outDIR/4 FQ_CHECKED/${runID}_${bcDir}_tr.fastq.gz" HEADCROP:9 MINLEN:250
		conda deactivate 
		

	# Run fastqc on the trimmed file.
	#
	log "Running FASTQC on trimmed files"
	conda activate "seqr_fastqc" || fail "Failed to activate conda env: seqr_fastqc"
	fastqc -o "$outDIR/5 QC2_RESULTS/" "$outDIR/4 FQ_CHECKED/${runID}_${bcDir}_tr.fastq.gz"
	conda deactivate

	
	# Map the reads to the reference genome.
	#
	log "Mapping reads"
	conda activate "seqr_general" || fail "Failed to activate conda env: seqr_general"
	minimap2 -ax map-ont "$refGenomeFILE" "$outDIR/4 FQ_CHECKED/${runID}_${bcDir}_tr.fastq.gz" > "$outDIR/6 MAPPED/${runID}_${bcDir}_tr.sam"
	
	
	# Convert to bam (binary sam), index, and write out some mapping statistics.
	#
	samtools view -b "$outDIR/6 MAPPED/${runID}_${bcDir}_tr.sam" > "$outDIR/6 MAPPED/${runID}_${bcDir}_tr.bam"
	samtools sort "$outDIR/6 MAPPED/${runID}_${bcDir}_tr.bam" > "$outDIR/6 MAPPED/${runID}_${bcDir}_trso.bam"
	samtools index "$outDIR/6 MAPPED/${runID}_${bcDir}_trso.bam"	# Creates a file ${runID}_${bcDir}_trso.bam.bai
			
	samtools stats "$outDIR/6 MAPPED/${runID}_${bcDir}_trso.bam" > "$outDIR/6 MAPPED/${runID}_${bcDir}_trso.mapstats.txt"
	samtools idxstats "$outDIR/6 MAPPED/${runID}_${bcDir}_trso.bam" > "$outDIR/6 MAPPED/${runID}_${bcDir}_trso.idxmapstats.txt"
			
	conda deactivate


	# Run freyja to identify SARS-CoV-2 variants.
	#
	log "Running freyja demixing"
	conda activate "seqr_freyja" || fail "Failed to activate conda env: seqr_freyja"
	freyja variants --variants "$outDIR/7 DEMIX/variants_${runID}_${bcDir}" --depths "$outDIR/7 DEMIX/depths_${runID}_${bcDir}" --ref "$refGenomeFILE" "$outDIR/6 MAPPED/${runID}_${bcDir}_trso.bam"
	freyja demix "$outDIR/7 DEMIX/variants_${runID}_${bcDir}.tsv" "$outDIR/7 DEMIX/depths_${runID}_${bcDir}" --output "$outDIR/7 DEMIX/demix_${runID}_${bcDir}" --barcodes "$usher_barcodesFILE"
	conda deactivate 
	
	log "Completed barcode $barcode2d"
done

log "SEQSARS completed successfully for run $runID"
exit 0

