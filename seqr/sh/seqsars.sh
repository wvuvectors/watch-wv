#!/bin/bash

#
# This script runs freyja analysis on barcoded COVID sequencing data from the ONT platform.
# It is generally called from an HPC scheduler script, but it can also be run separately.
# If you want to run this on your laptop, be warned that the processing requirements are not trivial.
#

START=$(date "+%F_%H-%M")

while getopts ":hi:r:b:g:o:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 1
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
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
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


# Create the output directory if it doesn't already exist.
if [ ! -d "$outDIR" ]; then
	mkdir "$outDIR"
fi

mkdir "$outDIR/0 RUN_FILES"
mkdir "$outDIR/1 BASECALLING"
mkdir "$outDIR/2 FQ_CONCAT"
mkdir "$outDIR/3 QC1_RESULTS"
mkdir "$outDIR/4 FQ_CHECKED"
mkdir "$outDIR/5 QC2_RESULTS"
mkdir "$outDIR/6 MAPPED"
mkdir "$outDIR/7 DEMIX"
mkdir "$outDIR/8 DASHBOARD"
mkdir "$outDIR/9 SRA"


for bcDir in "$runDIR/fastq_pass/barcode*"; do
	# Extract the barcode id to a variable $barcode.
	if [[ $bcDir =~ barcode(\d+) ]]; then
		barcode=${BASH_REMATCH[1]}
		barcode2d="$barcode"
		
		# Deal with a leading zero on the barcode number.
		if [[ $barcode =~ 0(\d) ]]; then
			barcode=${BASH_REMATCH[1]}
		fi

		# Concatenate all the fastq.gz files into one file.
		#
		ls "$runDIR/fastq_pass/$bcDir/*.fastq.gz" | xargs cat > "$outdir/2 FQ_CONCAT/$bcDir.fastq.gz"


		# Run fastqc on the original file.
		#
		conda activate "seqr_fastqc"
		fastqc -o "$outDIR/3 QC1_RESULTS/" "$outdir/2 FQ_CONCAT/$bcDir.fastq.gz"
		conda deactivate

		
		# Trim the reads in the concatenated fastq file. These are specific for COVID sequencing:
		# Remove the first 9 bases, and exclude any reads with length > 600nt.
		#
		conda activate "seqr_general"
		trimmomatic SE -phred64 "$outdir/2 FQ_CONCAT/${runID}_$bcDir.fastq.gz" "$outdir/4 FQ_CHECKED/${runID}_${bcDir}_tr.fastq.gz" HEADCROP:9 MINLEN:250
		conda deactivate 
		

		# Run fastqc on the trimmed file.
		#
		conda activate "seqr_fastqc"
		fastqc -o "$outDIR/5 QC2_RESULTS/" "$outdir/4 FQ_CHECKED/${runID}_${bcDir}_tr.fastq.gz"
		conda deactivate

		
		# Map the reads to the reference genome.
		#
		conda activate "seqr_general"
		minimap2 -ax map-ont "$refGenomeFILE" "$outdir/4 FQ_CHECKED/${runID}_${bcDir}_tr.fastq.gz" > "$outdir/6 MAPPED/${runID}_${bcDir}_tr.sam"
		
		
		# Convert to bam (binary sam), index, and write out some mapping statistics.
		#
		samtools view -b "$outdir/6 MAPPED/${runID}_${bcDir}_tr.sam" > "$outdir/6 MAPPED/${runID}_${bcDir}_tr.bam"
		samtools sort "$outdir/6 MAPPED/${runID}_${bcDir}_tr.bam" > "$outdir/6 MAPPED/${runID}_${bcDir}_trso.bam"
		samtools index "$outdir/6 MAPPED/${runID}_${bcDir}_trso.bam"	# Creates a file ${runID}_${bcDir}_trso.bam.bai
				
		samtools stats "$outdir/6 MAPPED/${runID}_${bcDir}_trso.bam" > "$outdir/6 MAPPED/${runID}_${bcDir}_trso.mapstats.txt"
		samtools idxstats "$outdir/6 MAPPED/${runID}_${bcDir}_trso.bam" > "$outdir/6 MAPPED/${runID}_${bcDir}_trso.idxmapstats.txt"
				
		conda deactivate


		# Run freyja to identify SARS-CoV-2 variants.
		#
		conda activate "seqr_freyja"
		freyja variants --variants "$outdir/7 DEMIX/variants_${runID}_${bcDir}" --depths "$outdir/7 DEMIX/depths_${runID}_${bcDir}" --ref "$refGenomeFILE" "$outdir/6 MAPPED/${runID}_${bcDir}_trso.bam"
		freyja demix "$outdir/7 DEMIX/variants_${runID}_${bcDir}.tsv" "$outdir/7 DEMIX/depths_${runID}_${bcDir}" --output "$outdir/7 DEMIX/demix_${runID}_${bcDir}" --barcodes "$usher_barcodesFILE"
		conda deactivate 

	fi
done

exit 0

