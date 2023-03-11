#! /bin/bash

indir="$1"
WD=$(pwd)

DBDIR="data"

# Get the current date and time
START=$(date "+%F_%H-%M")
TODAY=$(date "+%F_%H-%M")

# Write all output to log file
logf="logs/$START.patchr.txt"
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



echo "Ok, I'm about to start working with the watchdb files in $DBDIR/latest/." | tee -a "$logf"
echo "Backing up $DBDIR/latest/ to $DBDIR/latest_bk/." | tee -a "$logf"
cp -r "$DBDIR/latest/" "$DBDIR/latest_bk/"

echo "Creating a new version of watchdb in $DBDIR/incremental/$TODAY/." | tee -a "$logf"
echo "It will be updated with data from the run in $indir." | tee -a "$logf"

if [ -d "$DBDIR/incremental/$TODAY/" ]; then
	echo "An old version of $DBDIR/incremental/$TODAY/ already exists, so I'm removing the old directory first." | tee -a "$logf"
	rm -r "$DBDIR/incremental/$TODAY/"
fi
mkdir "$DBDIR/incremental/$TODAY/"



echo "Merging the run update tables with their corresponding watchdb tables." | tee -a "$logf"
tables=("abatch" "archive" "assay" "cbatch" "concentration" "control" "ebatch" "extraction" "rbatch" "sample")
for i in ${!tables[@]}
do
	table=${tables[$i]}
	updatef="$indir/updates/update.$table.txt"
	dbinf="$DBDIR/latest/watchdb.$table.txt"
	dboutf="$DBDIR/incremental/$TODAY/watchdb.$table.txt"
	
	#echo "Appending $updatef to watchdb.$table.txt." | tee -a "$logf"
	awk '(NR>1)' "$updatef" | cat "$dbinf" - > "$dboutf"
	#cp "$dboutf" "$dbinf"
done

# echo "******" | tee -a "$logf"
# echo "Running 2_mergeRun.pl." | tee -a "$logf"
# echo "******" | tee -a "$logf"
# 
# ./perl/2_mergeRun.pl "$indir" | tee -a "$logf"
# status="${PIPESTATUS[0]}"
# echo "" | tee -a "$logf"
# 
# if [[ "$status" != "0" ]]
# then
# 	echo "perl/2_mergeRun.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
# 	echo "!!!!!!!!" | tee -a "$logf"
# 	echo "patchr aborted during phase 2 (run data & watchdb merge)." | tee -a "$logf"
# 	echo "Run "| tee -a "$logf"
# 	echo "    ./sh/rollBack.sh $indir"| tee -a "$logf"
# 	echo "to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
# 	echo "!!!!!!!!" | tee -a "$logf"
# 	exit 1
# fi



echo "Copying the watchdb RESULT table into $DBDIR/incremental/$TODAY/." | tee -a "$logf"
cp "$DBDIR/latest/watchdb.result.txt" "$DBDIR/incremental/$TODAY/watchdb.result.txt"
echo "Updating the RESULT table taking into account any new assay data from $indir." | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 4_calculateResults.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"


./perl/4_calculateResults.pl "$DBDIR/incremental/$TODAY" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "4_calculateResults.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "I STRONGLY recommend deleting $DBDIR/incremental/$TODAY/ after exploring this error." | tee -a "$logf"
	#echo "Removing $DBDIR/watchdb/$TODAY/" | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 4 (calculating and updating results)." | tee -a "$logf"
	echo "Run "| tee -a "$logf"
	echo "    ./sh/rollBack.sh $indir"| tee -a "$logf"
	echo "to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "Copying from $DBDIR/incremental/$TODAY/ to $DBDIR/latest/." | tee -a "$logf"
for i in ${!tables[@]}
do
	table=${tables[$i]}
	cp "$DBDIR/incremental/$TODAY/watchdb.$table.txt" "$DBDIR/latest/watchdb.$table.txt"
done
cp "$DBDIR/incremental/$TODAY/watchdb.result.txt" "$DBDIR/latest/watchdb.result.txt"

echo "The most recent version of the watchdb can be found in two places:" | tee -a "$logf"
echo "    $DBDIR/incremental/$TODAY/" | tee -a "$logf"
echo "    and" | tee -a "$logf"
echo "    $DBDIR/latest/" | tee -a "$logf"
echo "$DBDIR/latest_bk/ contains the data from immediately before this update was applied." | tee -a "$logf"




echo "Updating the files for the Dashboard, using $DBDIR/latest/ as the source." | tee -a "$logf"

latest=("result" "sample")
for i in ${!latest[@]}
do
	table=${latest[$i]}
	cp "$DBDIR/latest/watchdb.$table.txt" "../dashboard/data/$table.txt"
done

resources=("county" "lab" "location" "wwtp")
for i in ${!resources[@]}
do
	table=${resources[$i]}
	cp "resources/watchdb.$table.txt" "../dashboard/data/$table.txt"
done






echo "" | tee -a "$logf"
echo "Finished patchr.sh" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "#############################################" | tee -a "$logf"
echo "" | tee -a "$logf"
