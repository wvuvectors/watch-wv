#! /bin/bash

indir="$1"
WD=$(pwd)

DBDIR="data"

# Get the current date and time
START=$(date "+%F_%H-%M")
TODAY=$(date "+%F_%H-%M")

# Write all output to log file
logf="patchr_log.$START.txt"
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
	echo "sh/prepRun.sh exited with error code $status and caused patchr to abort." | tee -a "$logf"
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
	echo "perl/1_compileRun.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
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


# This script also needs to merge entries from discontiguous events.
# For example, an extraction that used a concentration from an earlier run.
# In this case, the concentration will already exist in the watchdb, but not in the update.
# The extraction will have no concntration_id value and everything downstream will fail.

./perl/2_validateRun.pl "$indir" | tee -a "$logf"
status="${PIPESTATUS[0]}"
#echo "Status of validation: $status" | tee -a "$logf"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "perl/2_validateRun.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Most likely this script identifed overlapping or duplicate IDs. Check the following files for more info:" | tee -a "$logf"
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


echo "Ok, I'm about to start working with the watchdb files in $DBDIR/watchdb/LATEST/." | tee -a "$logf"
echo "Backing up $DBDIR/watchdb/LATEST/ to $DBDIR/watchdb/LATEST_BK/." | tee -a "$logf"
cp -r "$DBDIR/watchdb/LATEST/" "$DBDIR/watchdb/LATEST_BK/"

echo "Creating a new version of watchdb in $DBDIR/watchdb/$TODAY/." | tee -a "$logf"
echo "It will be updated with data from the run in $indir." | tee -a "$logf"

if [ -d "$DBDIR/watchdb/$TODAY/" ]; then
	echo "An old version of $DBDIR/watchdb/$TODAY/ already exists, so I'm removing the old directory first." | tee -a "$logf"
	rm -r "$DBDIR/watchdb/$TODAY/"
fi
mkdir "$DBDIR/watchdb/$TODAY/"


echo "Merging the run update tables with their corresponding watchdb tables." | tee -a "$logf"
tables=("abatch" "archive" "assay" "cbatch" "concentration" "control" "ebatch" "extraction" "rbatch")
for i in ${!tables[@]}
do
	table=${tables[$i]}
	updatef="$indir/updates/update.$table.txt"
	dbinf="$DBDIR/watchdb/LATEST/watchdb.$table.txt"
	dboutf="$DBDIR/watchdb/$TODAY/watchdb.$table.txt"
	
	#echo "Appending $updatef to watchdb.$table.txt." | tee -a "$logf"
	awk '(NR>1)' "$updatef" | cat "$dbinf" - > "$dboutf"
	#cp "$dboutf" "$dbinf"
done

echo "Copying the watchdb RESULT table into $DBDIR/watchdb/$TODAY/." | tee -a "$logf"
cp "$DBDIR/watchdb/LATEST/watchdb.result.txt" "$DBDIR/watchdb/$TODAY/watchdb.result.txt"
echo "Updating the RESULT table taking into account any new assay data from $indir." | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 3_calculateResults.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"


./perl/3_calculateResults.pl "$DBDIR/watchdb/$TODAY" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "perl/3_calculateResults.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "I STRONGLY recommend deleting $DBDIR/watchdb/$TODAY after exploring this error." | tee -a "$logf"
	#echo "Removing $DBDIR/watchdb/$TODAY/" | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 3 (Calculating and updating results)." | tee -a "$logf"
	echo "Run "| tee -a "$logf"
	echo "    ./sh/rollBack.sh $indir"| tee -a "$logf"
	echo "to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "$DBDIR/watchdb/$TODAY/ now contains the most recent version of the watchdb." | tee -a "$logf"
echo "Copying from $DBDIR/watchdb/$TODAY/ to $DBDIR/watchdb/LATEST/." | tee -a "$logf"
for i in ${!tables[@]}
do
	table=${tables[$i]}
	cp "$DBDIR/watchdb/$TODAY/watchdb.$table.txt" "$DBDIR/watchdb/LATEST/watchdb.$table.txt"
done
cp "$DBDIR/watchdb/$TODAY/watchdb.result.txt" "$DBDIR/watchdb/LATEST/watchdb.result.txt"

echo "$DBDIR/watchdb/LATEST/ now contains the most recent version of the watchdb." | tee -a "$logf"
echo "$DBDIR/watchdb/LATEST_BK/ contains the watchdb right before this update was applied." | tee -a "$logf"



#echo "$DBDIR/dashboard.LATEST.txt has *not* been modified!" | tee -a "$logf"
#mv "$DBDIR/dashboard/dashboard.$START.txt" "$DBDIR/dashboard/dashboard.$START.ABORTED.txt"
#cp "$DBDIR/dashboard/dashboard.$START.txt" "$DBDIR/dashboard.LATEST.txt"

echo "" | tee -a "$logf"
echo "Finished patchr.sh" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "#############################################" | tee -a "$logf"
echo "" | tee -a "$logf"
