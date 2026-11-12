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
SEWER replace script.
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
	echo "Please use the -i argument to provide a valid path to the merge folder "
	echo "that you would like to replace. This is required regardless of any other settings."
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
	echo "as the path to the merge folder, but either I can't find it, or I don't have permission"
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


# Make the log directory, if necessary. It should not be.
#
if [[ ! -d "$OUTDIR/logs" ]];then
	mkdir "$OUTDIR/logs"
	echo "This folder was created by sewer-replace $logstamp." > "$OUTDIR/logs/README.txt"
fi


# Write all output to log file.
#
if [[ ! -f "$logf" ]];then
	logf="$OUTDIR/logs/sewer-replace_$RUNID.$logstamp.log"
fi

touch "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER replace latest." | tee -a "$logf"
echo "This script replaces the existing database with a merge of new data." | tee -a "$logf"
echo "See $logf for general warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Options for this run of SEWER replace are as follows:" | tee -a "$logf"
echo "   Workflow ID:   $RUNID" | tee -a "$logf"
echo "   Merge folder:  $INDIR" | tee -a "$logf"
echo "   Output folder: $OUTDIR" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "#" | tee -a "$logf"
echo "1. Backing up main database $OUTDIR/latest to $OUTDIR/latest_bk/ before replacing it." | tee -a "$logf"
echo "" | tee -a "$logf"
cp -r "$OUTDIR/latest/" "$OUTDIR/latest_bk/"
PREVDATE=$(head -n 1 "$OUTDIR/latest/README.txt")
echo "$PREVDATE" > "$OUTDIR/latest_bk/README.txt"
echo "" >> "$OUTDIR/latest_bk/README.txt"
echo "#" >> "$OUTDIR/latest_bk/README.txt"
echo "This folder contains backups of the watchdb tables from immediately before $UPDAY." >> "$OUTDIR/latest_bk/README.txt"
echo "#" >> "$OUTDIR/latest_bk/README.txt"
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"

tables=("cbatch" "concentration" "ebatch" "extraction" "qassay" "qbatch" "qblock" "qcontrol" "sample" "result")

echo "#" | tee -a "$logf"
echo "2. Moving files from $INDIR to $OUTDIR/latest." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "#" | tee -a "$logf"
echo "Running 4A_replaceLatest.R." | tee -a "$logf"
outval=$(Rscript R/4A_replaceLatest.R "$OUTDIR/latest" "$INDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "4A_replaceLatest.R exited with error code $status and caused sewer-replace to abort." | tee -a "$logf"
	echo "Arguments: $OUTDIR/latest $INDIR" | tee -a "$logf"
	echo "sewer replace aborted during replacement. It is likely that the tables in $OUTDIR/latest " | tee -a "$logf"
	echo "have been corrupted. We recommend that you replace these files with the backup copy in "| tee -a "$logf"
	echo "$OUTIDR/latest_bk and fix the issue before trying again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
# for table in "${tables[@]}"
# do
# 	cp "$INDIR/merge.$table.txt" "$OUTDIR/latest/watchdb.$table.txt"
# 	echo "Copied $INDIR/merge.$table.txt to $OUTDIR/latest/watchdb.$table.txt"
# done

echo "" | tee -a "$logf"
echo "Updating the README file in $OUTDIR/latest." | tee -a "$logf"

echo "$UPDAY" > "$OUTDIR/latest/README.txt"
echo "" >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"
echo "This folder contains the most recent watchdb tables." >> "$OUTDIR/latest/README.txt"
echo "The latest run id included here is $RUNID." >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"

echo "Done." | tee -a "$logf"

# if [[ $QC == 1 ]];then
# 	echo "" | tee -a "$logf"
# 	echo "*****" | tee -a "$logf"
# 	echo "3. Generating QC outcomes for merge $RUNID." | tee -a "$logf"
# 	echo "" | tee -a "$logf"
# 	echo "These files can be found in $INDIR/qc." | tee -a "$logf"
# 	mkdir "$INDIR/qc"
# 	
# 	echo "#" | tee -a "$logf"
# 	echo "Running 4B_qcREPLACE.R." | tee -a "$logf"
# 	Rscript R/4B_qcREPLACE.R "$INDIR"
# 	status="${PIPESTATUS[0]}"
# 	if [[ "$status" != "0" ]];then
# 		echo "!!!!!!!!" | tee -a "$logf"
# 		echo "4B_qcREPLACE.R exited with error code $status and caused sewer-replace to abort." | tee -a "$logf"
# 		echo "Arguments: $INDIR" | tee -a "$logf"
# 		echo "sewer-merge aborted during QC. This is likely a problem with the QC subroutine, not " | tee -a "$logf"
# 		echo "the replace itself, but no promises. We recommend that you inspect the merge to ensure "| tee -a "$logf"
# 		echo "nothing is amiss."| tee -a "$logf"
# 		echo "!!!!!!!!" | tee -a "$logf"
# 		exit 1
# 	fi
# 	echo "Done." | tee -a "$logf"
# fi

echo "" | tee -a "$logf"
echo "SEWER replace run $logstamp for $RUNID, started on $UPDAY, has left the chat." | tee -a "$logf"
