#! /bin/bash

OUTDIR="../data"
QC=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INDIR="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -o|--output)
      OUTDIR="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -q|--noqc)
      QC=0
      shift 1 # Shift past the flag
      ;;
    -l|--logfile)
      logf="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -h|--help)
      cat << EOF
SEWER merge script.
EOF
		exit 0
		;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done


if [[ -z "$INDIR" ]];then
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : no input directory provided (-i). This is a FATAL error."
	echo "Please use the -i argument to provide a valid path to the update folder "
	echo "that you would like to merge. This is required regardless of any other settings."
	echo "Use -h or --help for more information about how to run this script. (Eg, $0 --help)."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi

if [[ ! -d "$INDIR" ]];then
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : the given input directory does not exist or is not readable. This is a FATAL error."
	echo "You have tried to use"
	echo "$INDIR"
	echo "as the path to the update folder, but either I can't find it, or I don't have permission"
	echo "to access it. Use -h or --help for more information about how to run this script."
	echo "(Eg, $0 --help)."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi

logstamp=$(date "+%Y%m%d%H%M")
UPDAY=$(date "+%B %d, %Y at %T")

# Get the run ID, which is the name of input folder.
#
RUNID=$(basename "$INDIR")


# Make the parent output directories, if necessary.
#
if [[ ! -d "$OUTDIR/" ]];then
	mkdir "$OUTDIR/"
	echo "This folder was created by sewer-merge $logstamp." > "$OUTDIR/README.txt"
fi

if [[ ! -d "$OUTDIR/merges" ]];then
	mkdir "$OUTDIR/merges"
	echo "This folder was created by sewer-merge $logstamp." > "$OUTDIR/merges/README.txt"
fi

if [[ ! -d "$OUTDIR/logs" ]];then
	mkdir "$OUTDIR/logs"
	echo "This folder was created by sewer-merge $logstamp." > "$OUTDIR/logs/README.txt"
fi


# Look for RUNID in any of the existing batch files. If it exists, do not proceed. Prompt the user to 
# make sure the current RUNID is unique.
#
# grep exits with a 0 if a match IS found, 1 if a match is NOT found, and 2 on error.
if [[ -d "$OUTDIR/latest" ]];then
	FOUND=$(grep -wi "$RUNID" "$OUTDIR/latest/watchdb.qbatch.txt")
	STATUS=$?
	if [[ $STATUS == 0 ]];then
		echo ""
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo "ERROR : a workflow labeled $RUNID has already been merged. This is a FATAL error."
		echo "Every run of the laboratory workflow is assigned a unique ID derived from the name of"
		echo "the folder that contains the raw workflow files. You have tried to merge an update using"
		echo "$RUNID"
		echo "as the workflow ID, but that ID already exists in the current database here"
		echo "$OUTDIR/latest"
		echo "Please change the workflow run folder name to something unique. If you have multiple"
		echo "workflows from the same calendar day, add a suffix to the folder name (e.g., 2026-06-22b)"
		echo "and then run sewer-make before trying the merge again."
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo ""
		exit 1
	fi
fi

# Make a directory to hold the merge that will be generated.
#
MERGEDIR="$OUTDIR/merges/$RUNID"
if [[ ! -d "$MERGEDIR" ]];then
	mkdir "$MERGEDIR"
else
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : a merge folder labeled $RUNID already exists in $OUTDIR/merges. This is a FATAL error."
	echo "Every run of the laboratory workflow is assigned a unique ID derived from the name of"
	echo "the folder that contains the raw workflow files. You have tried to make an update using"
	echo "$RUNID"
	echo "as the workflow ID, but a folder with that name already exists in the merges directory here"
	echo "$OUTDIR/merges"
	echo "Oddly, an update with this run ID does NOT appear in the existing database, so this is likely"
	echo "a result of an incomplete sewer run, where an update was merged but not replaced."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi


# Write all output to log file.
#
if [[ ! -f "$logf" ]];then
	logf="$OUTDIR/logs/sewer-merge_$RUNID.$logstamp.log"
fi

echo ""

touch "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER merge update." | tee -a "$logf"
echo "This script merges an update with the existing database." | tee -a "$logf"
echo "See $logf for general warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Options for this run of SEWER merge are as follows:" | tee -a "$logf"
echo "   Workflow ID:     $RUNID" | tee -a "$logf"
echo "   Workflow folder: $INDIR" | tee -a "$logf"
echo "   Merge folder:    $MERGEDIR" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "*****" | tee -a "$logf"
echo "1. Merging update $RUNID with the main database $OUTDIR/latest." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "#" | tee -a "$logf"
echo "Running 3A_mergeUpdate.R." | tee -a "$logf"
dup_count=$(Rscript R/3A_mergeUpdate.R "$INDIR" "$OUTDIR/latest" "$MERGEDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "3A_mergeUpdate.R exited with error code $status and caused sewer-merge to abort." | tee -a "$logf"
	echo "Arguments: $INDIR $OUTDIR/latest $MERGEDIR" | tee -a "$logf"
	echo "sewer-merge aborted during phase A (applying the update)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs and address the errors. " | tee -a "$logf"
	echo "You probably need to delete $MERGEDIR and start over."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "*****" | tee -a "$logf"
echo "2. Generating result file for merge $RUNID." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "#" | tee -a "$logf"
echo "Running 3B_generateResults.R." | tee -a "$logf"
result_count=$(Rscript R/3B_generateResults.R "$MERGEDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "3B_generateResults.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $MERGEDIR" | tee -a "$logf"
	echo "sewer-merge aborted during phase B (generating a result table)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs and address the errors. " | tee -a "$logf"
	echo "You probably need to delete $MERGEDIR and start over."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"

# 
# if [[ $QC == 1 ]];then
# 	echo "*****" | tee -a "$logf"
# 	echo "3. Generating QC outcomes for merge $RUNID." | tee -a "$logf"
# 	echo "" | tee -a "$logf"
# 	echo "These files can be found in $MERGEDIR/qc." | tee -a "$logf"
# 	mkdir "$MERGEDIR/qc"
# 
# 	echo "#" | tee -a "$logf"
# 	echo "Running 3C_qcMERGE.R." | tee -a "$logf"
# 	Rscript R/3C_qcMERGE.R "$MERGEDIR"
# 	status="${PIPESTATUS[0]}"
# 	if [[ "$status" != "0" ]]
# 	then
# 		echo "!!!!!!!!" | tee -a "$logf"
# 		echo "3C_qcMERGE.R exited with error code $status and caused sewer-merge to abort." | tee -a "$logf"
# 		echo "Arguments: $MERGEDIR" | tee -a "$logf"
# 		echo "sewer-merge aborted during QC. This is likely a problem with the QC subroutine, not " | tee -a "$logf"
# 		echo "the merge itself, but no promises. We recommend that you inspect the merge to ensure "| tee -a "$logf"
# 		echo "nothing is amiss."| tee -a "$logf"
# 		echo "!!!!!!!!" | tee -a "$logf"
# 		exit 1
# 	fi
# 	echo "Done." | tee -a "$logf"
# fi

echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER merge run $logstamp for $RUNID, started on $UPDAY, has left the chat." | tee -a "$logf"
echo "" | tee -a "$logf"
