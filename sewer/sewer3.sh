#! /bin/bash

while getopts ":hi:o:m:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 1
			;;
		i)
			UPDIR=$OPTARG
			;;
		o)
			OUTDIR=$OPTARG
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


if [ -z "$UPDIR" ];then
	echo "FATAL : no update directory provided to use as input (-i)."
	exit 1
fi

if [ ! -d "$UPDIR" ];then
	echo "FATAL : the update directory does not exist or is not readable."
	exit 1
fi


if [ -z "$OUTDIR" ];then
	echo "FATAL : no output directory provided (-o)."
	exit 1
fi


# Get the current date and time for logging purposes.
#
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")
TS=$(date +%s)

# Get the run ID, which is the name of input folder.
#
RUNID=$(basename "$UPDIR")
INCID="$RUNID"
#INCID="$TS.$RUNID"

# Make the output directory and all required subdirs, as necessary.
#
if [ ! -d "$OUTDIR/" ];then
	mkdir "$OUTDIR/"
	echo "This folder was created on $UPDAY." > "$OUTDIR/README.txt"
fi

if [ ! -d "$OUTDIR/latest" ];then
	mkdir "$OUTDIR/latest"
fi
if [ ! -d "$OUTDIR/incremental" ];then
	mkdir "$OUTDIR/incremental"
fi


# Look for RUNID in any of the existing batch files. If it exists, do not proceed. Prompt the user to 
# make sure the current RUNID is unique.
#
# grep exits with a 0 if a match IS found, 1 if a match is NOT found, and 2 on error.
if [ -d "$OUTDIR/latest" ];then
	FOUND=$(grep -wi "$RUNID" "$OUTDIR/latest/watchdb.qbatch.txt")
	STATUS=$?
	if [ $STATUS == 0 ];then
		echo "FATAL : Every run of the laboratory workflow is assigned an ID using the name of the folder that contains the raw files."
		echo "FATAL : It looks like a run labeled $RUNID has already been processed."
		echo "FATAL : This can occur if you forgot to change the name of the run directory."
		echo "FATAL : It might also occur if you have multiple runs on the same calendar day."
		echo "FATAL : In this case, add a suffix to the name of the run directory (e.g., 2026-06-22b) and run sewer again."
		exit 1
	fi
fi


# Make a directory to hold the update that will be generated.
#
UPDIR="$OUTDIR/updates/$RUNID"
if [ ! -d "$OUTDIR/updates" ]
then
	mkdir "$OUTDIR/updates/"
fi
mkdir "$UPDIR"
mkdir "$UPDIR/run_files"


# Write all output to log file in the update dir.
#
logf="$UPDIR/sewer.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "SEWER: Surveillance & Exploration of Wastewater to inform Epidemiological Response." | tee -a "$logf"
echo "Initiated from sewer3.sh." | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "Run ID: $RUNID" | tee -a "$logf"
echo "Update dir: $UPDIR" | tee -a "$logf"
echo "Output dir: $OUTDIR" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Now I'll start the process of applying the update from $UPDIR to $OUTDIR/latest/." | tee -a "$logf"
echo "First, I'll back up $OUTDIR/latest/ to $OUTDIR/latest_bk/ just in case we need to roll all the way back." | tee -a "$logf"
cp -r "$OUTDIR/latest/" "$OUTDIR/latest_bk/"
PREVDATE=$(head -n 1 "$OUTDIR/latest/README.txt")
echo "$PREVDATE" > "$OUTDIR/latest_bk/README.txt"
echo "" >> "$OUTDIR/latest_bk/README.txt"
echo "#" >> "$OUTDIR/latest_bk/README.txt"
echo "This folder contains backups of the watchdb tables from before $UPDAY." >> "$OUTDIR/latest_bk/README.txt"
echo "#" >> "$OUTDIR/latest_bk/README.txt"
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Now I'll append the update tables from $UPDIR to the existing data in $OUTDIR/latest." | tee -a "$logf"
echo "The merge will first be written to $OUTDIR/incremental/$RUNID." | tee -a "$logf"
echo "If successful, the new files will be copied from $OUTDIR/incremental/$RUNID to $OUTDIR/latest." | tee -a "$logf"
echo "" | tee -a "$logf"

mkdir "$OUTDIR/incremental/$RUNID"

echo "******" | tee -a "$logf"
echo "Running 3A_applyUpdate.R." | tee -a "$logf"
dup_count=$(Rscript R/3A_applyUpdate.R "$UPDIR" "$OUTDIR/latest" "$OUTDIR/incremental/$RUNID")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "3A_applyUpdate.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR $OUTDIR/latest $OUTDIR/incremental/$RUNID" | tee -a "$logf"
	echo "sewer aborted during phase 6 (applying the update)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs and address the errors. " | tee -a "$logf"
	echo "You may need to delete $UPDIR and/or $OUTDIR/incremental/$RUNID, and start over."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "It seems like the merge was successful!" | tee -a "$logf"
echo "Before we finish, though, I need to pre-calculate some processor-intensive data." | tee -a "$logf"
echo "This is primarily aimed at optimizing realtime dashboard performance and generates a file called watchdb.result.txt." | tee -a "$logf"
echo "The table in this file is an amalgam of data tables and consequently will contain a certain level of redundancy." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 3B_generateResults.R." | tee -a "$logf"
result_count=$(Rscript R/3B_generateResults.R "$OUTDIR/incremental/$RUNID")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "3B_generateResults.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $OUTDIR/incremental/$RUNID" | tee -a "$logf"
	echo "sewer aborted during phase 7 (generating the result table)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs and address the errors. " | tee -a "$logf"
	echo "You may need to delete $OUTDIR/incremental/$RUNID and start over."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Final copy of files from $OUTDIR/incremental/$RUNID to $OUTDIR/latest." | tee -a "$logf"
cp "$OUTDIR/incremental/$RUNID/watchdb."*.txt "$OUTDIR/latest/"
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Updating the README file." | tee -a "$logf"
echo "$UPDAY" > "$OUTDIR/latest/README.txt"
echo "" >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"
echo "This folder contains the most recent watchdb tables." >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "The most recent version of the watchdb can be found in $OUTDIR/latest." | tee -a "$logf"
echo "The version of the watchdb from just before this update was applied can be found in $OUTDIR/latest_bk." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "sewer run of $UPDAY has left the chat." | tee -a "$logf"
echo "" | tee -a "$logf"

