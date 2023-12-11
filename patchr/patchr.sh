#! /bin/bash

indir="$1"
WD=$(pwd)

DBDIR="data"

# Get the current date and time
START=$(date "+%F_%H-%M")
TODAY=$(date "+%F_%H-%M")

# Write all output to log file
logf="logs/patchr.$START.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "Initiated patchr.sh" | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "In the event of a catastrophic error, run"| tee -a "$logf"
echo "  ./sh/rollBack.sh $indir"| tee -a "$logf"
echo "at any point to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
echo "" | tee -a "$logf"



echo "******" | tee -a "$logf"
echo "Running 1_queryBatches.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"

upfiles=$(./perl/1_queryBatches.pl "$indir")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "1_queryBatches.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 1 (batch identification)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

if [ -z "$upfiles" ]
then
	echo "******" | tee -a "$logf"
	echo "1_queryBatches.pl found no new batches in $rundir." | tee -a "$logf"
	echo "As a result, there is nothing for patchr to do at the current time." | tee -a "$logf"
	echo "******" | tee -a "$logf"

	echo "All done. patchr will now exit."
	echo ""
	exit 1
fi


echo "1_queryBatches.pl identified unprocessed batch files, so an update will be prepared." | tee -a "$logf"
echo "Preparing an update to WaTCH now." | tee -a "$logf"

update_dir="data/updates/$START/"
echo "Creating update folder $update_dir to hold the results." | tee -a "$logf"
if [ ! -d "$update_dir" ]
then
	mkdir "$update_dir"
fi
echo "Recording batch files to process in $update_dir/update.batch_files.txt." | tee -a "$logf"
echo "$upfiles" > "$update_dir/update.batch_files.txt"


echo "******" | tee -a "$logf"
echo "Running 2_compileUpdate.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/2_compileUpdate.pl "$update_dir" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "2_compileUpdate.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 2 (run data compilation)." | tee -a "$logf"
	echo "Run "| tee -a "$logf"
	echo "    ./sh/rollBack.sh $indir"| tee -a "$logf"
	echo "to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi



