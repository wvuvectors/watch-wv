#! /bin/bash

indir="$1"
WD=$(pwd)


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
echo "In the event of a catastrophic error, run"| tee -a "$logf"
echo "  ./sh/rollBack.sh $indir"| tee -a "$logf"
echo "at any point to restore the original files. Then fix the error(s) and run patchr again."| tee -a "$logf"
echo "" | tee -a "$logf"


./sh/prepRun.sh "$indir" | tee -a "$logf"
status="$?"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
#	echo "sh/prepRun.sh was unable to locate a ddPCR results file." | tee -a "$logf"
#	echo "This csv file is required for run processing." | tee -a "$logf"
	echo "sh/prepRun.sh exited with error code $status  and caused patchr to abort." | tee -a "$logf"
	echo "Please see $logf for more information." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 0 (run file prep)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


./perl/1_compileRun.pl "$indir" | tee -a "$logf"
status="$?"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "perl/1_compileRun.pl reported a fatal error and caused patchr to abort." | tee -a "$logf"
	echo "Please see $logf for more information." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 1 (run data compilation)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "Backing up watchdb LATEST/ directory to LATEST_BK/." | tee -a "$logf"
cp -r "data/watchdb/LATEST/" "data/watchdb/LATEST_BK/"

if [ -d "data/watchdb/$TODAY/" ]; then
	rm -r "data/watchdb/$TODAY/"
fi
mkdir "data/watchdb/$TODAY/"

tables=("abatch" "archive" "assay" "cbatch" "concentration" "control" "ebatch" "extraction" "rbatch")
for i in ${!tables[@]}
do
	table=${tables[$i]}
	updatef="$indir/updates/update.$table.txt"
	dbinf="data/watchdb/LATEST/watchdb.$table.txt"
	dboutf="data/watchdb/$TODAY/watchdb.$table.txt"

	echo "Appending $updatef to watchdb.$table.txt." | tee -a "$logf"
	awk '(NR>1)' "$updatef" | cat "$dbinf" - > "$dboutf"
	cp "$dboutf" "$dbinf"
done

#./perl/2_validateWaTCH.pl | tee -a "$logf"
#status="$?"


