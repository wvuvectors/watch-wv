#!/bin/bash

#SBATCH -J watch_seqsars_run
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -p standby
#SBATCH -t 4:00:00

START=$(date "+%F_%H-%M")

while getopts ":hr:i:f:o:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 1
			;;
		r)
			runID=$OPTARG
			;;
		i)
			readsDIR=$OPTARG
			;;
		f)
			fr_barcodesFILE=$OPTARG
			;;
		o)
			outDIR=$OPTARG
			;;
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

if [ -z "$runID" ]; then
	echo "Missing a valid run ID (-r)!"
	echo "I'm going to use a default \'RUNON_$START\' but be aware:"
	echo "this run ID contains the date the analysis was started, not the date of sequencing."
	runID="RUNON_$START"
fi

if [ -z "$outDIR" ]; then
	echo "Missing a valid output directory (-o)!"
	echo "Writing output to the default dir of \'output/$runID/\'."
	echo "You MUST move this dir out of the git repository BEFORE committing! MUST!"
	echo "If you fail to do so, the data push to github will FAIL (it's too big)."
	outDIR="output/$runID"
fi


if [ -z "$readsDIR" ]; then
	echo "Missing a valid input directory of reads (-i)!"
	echo "This is usually the directory that *contains* the \'fastq_pass\' folder (among others)."
	echo "This is a required parameter with no default, so I'm forced to quit. Sorry."
	exit 1
fi

if [ -z "$fr_barcodesFILE" ]; then
	echo "Missing a valid freyja barcodes file (-f)!"
	echo "This is usually located in the project resources folder and has a \'.feather\' extension."
	echo "This is a required parameter with no default, so I'm forced to quit. Sorry."
	exit 1
fi

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
./sh/seqsars.sh -r "$runID" -i "$readsDIR" -f "$fr_barcodesFILE" -o "$outDIR" #-c "$SCRATCH/conda"
#usher_barcodes09_05_2024-00-47.feather


