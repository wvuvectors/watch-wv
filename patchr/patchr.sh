#! /bin/bash

indir="$1"
WD=$(pwd)

DBDIR="data"

# Get the current date and time
START=$(date "+%F_%H-%M")
TODAY=$(date "+%F_%H-%M")

# Write all output to log file
logf="logs/patchr/patchr.$START.log"
if [ -f "$logf" ]; then
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
echo "Running prepRun.sh." | tee -a "$logf"
echo "******" | tee -a "$logf"


./sh/prepRun.sh "$indir" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "prepRun.sh exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 0 (run file prep)." | tee -a "$logf"
	echo "Run "| tee -a "$logf"
	echo "    ./sh/rollBack.sh $indir"| tee -a "$logf"
	echo "to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi




echo "******" | tee -a "$logf"
echo "Running 1_compileRun.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/1_compileRun.pl "$indir" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "1_compileRun.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 1 (run data compilation)." | tee -a "$logf"
	echo "Run "| tee -a "$logf"
	echo "    ./sh/rollBack.sh $indir"| tee -a "$logf"
	echo "to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi




echo "******" | tee -a "$logf"
echo "Running 2_validateRun.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"


./perl/2_validateRun.pl "$indir" | tee -a "$logf"
status="${PIPESTATUS[0]}"
#echo "Status of validation: $status" | tee -a "$logf"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "2_validateRun.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Most likely this script identifed duplicate IDs or collisions (IDs present in both the existing db and the update). Check the following files for more info:" | tee -a "$logf"
	echo "    $indir/updates/_collisions.txt" | tee -a "$logf"
	echo "    $indir/updates/_rundups.txt" | tee -a "$logf"
	echo "    $indir/updates/_watchdups.txt" | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 2 (run validation)." | tee -a "$logf"
	echo "Run "| tee -a "$logf"
	echo "    ./sh/rollBack.sh $indir"| tee -a "$logf"
	echo "to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


