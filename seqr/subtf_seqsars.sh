#!/bin/bash

#SBATCH -J watch_seqsars_run
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -p standby
#SBATCH -t 4:00:00

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
	echo "This is usually the parent run folder that contains the \'fastq_pass\' folder (among others)."
	echo "This is a required parameter with no default, so I'm forced to quit immediately. Sorry."
	exit 1
fi


# Required parameters for which I can try to provide a default value.
#

if [ -z "$usher_barcodesFILE" ]; then
	usher_barcodesFILE="resources/usher_barcodes.csv"
fi

if [ -z "$refGenomeFILE" ]; then
	refGenomeFILE="resources/NC_045512_Hu-1.fasta"
fi

if [ -z "$runID" ]; then
	echo "You have not provided a custom run ID (-r)."
	echo "I'm going to use the default \'RUNON_$START\' but be aware:"
	echo "this run ID contains the date the analysis was started, not the date of sequencing."
	runID="RUNON_$START"
fi

if [ -z "$outDIR" ]; then
	echo "You have not provided a custom output directory (-o)."
	echo "I will write all output to the default dir of \'output/$runID/\'"
	outDIR="output/$runID"
fi

# if [ -z "$condaDIR" ]; then
# 	echo "Missing a valid conda environment directory (-c)!"
# 	echo "I'll try to use the current directory, but if you get a conda environment error that might be why."
# 	condaDIR="."
# fi


cd $SCRATCH/watch-wv/seqr/


# This process utilizes conda environments, so we need to load the conda profile script before anything else.
# We do it here because the path is often machine-specific.
#
source /shared/software/conda/etc/profile.d/conda.sh

# We call a script here that contains all the run code for freyja analysis, etc.
# This slightly awkward structure allows us to port the actual analysis between machines easily.
# It also means the first 50 or so lines of this script are duplicated in the called script; 
# however, we don't anticipate those lines will change often (if ever).
#
./sh/seqsars.sh -i "$runDIR" -r "$runID" -b "$usher_barcodesFILE" -g "$refGenomeFILE" -o "$outDIR" #-c "$SCRATCH/conda"


