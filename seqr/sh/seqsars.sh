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

# Conda setup (specific to my scratch drive, should be updated before merging with main branch)
source /shared/software/conda/etc/profile.d/conda.sh || fail "Could not find conda.sh; Conda not initialized"
conda_activate() {
	set +u
	conda activate "$1"
	set -u 
}

conda_deactivate() {
	set +u 
	conda deactivate
	set -u 
}

# Create the output directory if it doesn't already exist.
log "Creating output directory structure in $outDIR"
if [ ! -d "$outDIR" ]; then
	mkdir "$outDIR"
fi

mkdir -p "$outDIR/0_RUN_FILES"
mkdir -p "$outDIR/1_BASECALLING"
mkdir -p "$outDIR/2_FQ_CONCAT"
mkdir -p "$outDIR/3_QC1_RESULTS"
mkdir -p "$outDIR/4_FQ_CHECKED"
mkdir -p "$outDIR/5_QC2_RESULTS"
mkdir -p "$outDIR/6_MAPPED"
mkdir -p "$outDIR/7_DEMIX"
mkdir -p "$outDIR/8_DASHBOARD"
mkdir -p "$outDIR/9_SRA"

# Barcode discovery
log "Discovering barcode directories"
barcode_dirs=( "$runDIR/fastq_pass"/barcode* )
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
	concat_fastq="$outDIR/2_FQ_CONCAT/$bcBase.fastq.gz"
	trimmed_fastq="$outDIR/4_FQ_CHECKED/${runID}_${bcBase}_tr.fastq.gz"
	
	# Check that FASTQ files exist before concatenation
	fastq_files=( "$bcDir"/*.fastq.gz )
	(( ${#fastq_files[@]} > 0 )) || fail "No FASTQ files found in $bcDir"

	# Concatenate all the fastq.gz files into one file.
	#
	log "Concatenating FASTQs for barcode $barcode2d"
	
	# Input FASTQs
	fastq_files=( "$bcDir"/*.fastq.gz )
	(( ${#fastq_files[@]} > 0 )) || fail "No FASTQ files found in $bcDir"
	
	# Output file 
	concat_fastq="$outDIR/2_FQ_CONCAT/$bcBase.fastq.gz"
	
	cat "${fastq_files[@]}" > "$concat_fastq"


	# Run fastqc on the original file.
	#
	log "Running FastQC on raw files"
	conda_activate "seqr_fastqc" || fail "Failed to activate conda env: seqr_fastqc"
	fastqc -o "$outDIR/3_QC1_RESULTS/" "$concat_fastq"
	conda_deactivate

		
	# Trim the reads in the concatenated fastq file. These are specific for COVID sequencing:
	# Remove the first 9 bases, and exclude any reads with length > 600nt.
	#
	log "Trimming reads"
	conda_activate "seqr_general" || fail "Failed to activate conda env: seqr_general"
	trimmed_fastq="$outDIR/4_FQ_CHECKED/${runID}_${bcBase}_tr.fastq.gz"
	trimmomatic SE -phred64 "$concat_fastq" "$trimmed_fastq" HEADCROP:9 MINLEN:250
	conda_deactivate 
		

	# Run fastqc on the trimmed file.
	#
	log "Running FASTQC on trimmed files"
	conda_activate "seqr_fastqc" || fail "Failed to activate conda env: seqr_fastqc"
	fastqc -o "$outDIR/5_QC2_RESULTS/" "$trimmed_fastq"
	conda_deactivate

	
	# Map the reads to the reference genome.
	#
	log "Mapping reads"
	conda_activate "seqr_general" || fail "Failed to activate conda env: seqr_general"
	sam="$outDIR/6_MAPPED/${runID}_${bcBase}_tr.sam"
	bam="$outDIR/6_MAPPED/${runID}_${bcBase}_tr.bam"
	bam_sorted="$outDIR/6_MAPPED/${runID}_${bcBase}_trso.bam"

	minimap2 -ax map-ont "$refGenomeFILE" "$trimmed_fastq" > "$sam"
	
	
	# Convert to bam (binary sam), index, and write out some mapping statistics.
	#
	samtools view -b "$sam" > "$bam"
	samtools sort "$bam" > "$bam_sorted"
	samtools index "$bam_sorted"	# Creates a file ${runID}_${bcDir}_trso.bam.bai
			
	samtools stats "$bam_sorted" > "$outDIR/6_MAPPED/${runID}_${bcBase}_trso.mapstats.txt"
	samtools idxstats "$bam_sorted" > "$outDIR/6_MAPPED/${runID}_${bcBase}_trso.idxmapstats.txt"
			
	conda_deactivate


	# Run freyja to identify SARS-CoV-2 variants.
	#
	log "Running freyja demixing"
	
	# Defining filenames
	demix_dir="$outDIR/7_DEMIX"
	variants_prefix="$demix_dir/variants_${runID}_${bcBase}"
	depths_file="$demix_dir/depths_${runID}_${bcBase}"
	demix_out="$demix_dir/demix_${runID}_${bcBase}"
	bam_file="$outDIR/6_MAPPED/${runID}_${bcBase}_trso.bam"
	
	# Ensuring DEMIX directory exists
	mkdir -p "$demix_dir"

	conda_activate "seqr_freyja" || fail "Failed to activate conda env: seqr_freyja"
	freyja variants \
		--variants "$variants_prefix" \
		--depths "$depths_file" \
		--ref "$refGenomeFILE" \
		"$bam_file"
		
	# Quality gate before demixing
	log "Checking coverage quality for $bcBase"

	# 1) Check mapped reads in BAM
	mapped_reads=$(samtools view -c -F 4 "$bam_file" || echo 0)

	if [[ "$mapped_reads" -lt 5000 ]]; then
		log "[WARN] Skipping freyja demix for $bcBase: low mapped reads ($mapped_reads)"
		conda_deactivate
		continue
	fi

	# 2) Check depths file exists and is non-empty
	if [[ ! -s "$depths_file" ]]; then
		log "[WARN] Skipping freyja demix for $bcBase: empty depths file"
		conda_deactivate
		continue
	fi

	# 3) Check depths file has enough lines
	min_lines=5000
	depth_lines=$(wc -l < "$depths_file" || echo 0)

	if [[ "$depth_lines" -lt $min_lines ]]; then
		log "[WARN] Skipping freyja demix for $bcBase: insufficient depth data ($depth_lines lines)"
		conda_deactivate
		continue
	fi
	
	freyja demix \
		"${variants_prefix}.tsv" \
		"$depths_file" \
		--output "$demix_out" \
		--barcodes "$usher_barcodesFILE"
	
	conda_deactivate 
	
	log "Completed barcode $barcode2d"
done

log "SEQSARS completed successfully for run $runID"
exit 0

