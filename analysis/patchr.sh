#! /bin/bash

indir="$1"
WD=$(pwd)


# Get the current date and time
START=$(date + '%Y-%m-%d_%H-%M')
TODAY=$(date + '%Y-%m-%d')

# Write all output to log file
logf="patchr_log.$START.txt"
touch "$logf"

echo "#############################################" | tee -a "$logf"
echo "Initiated patchr.sh" | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"


status=$(./sh/prepRun.sh "$indir" | tee -a "$logf")

if [[ "$status" == "0" ]]
then
	echo "sh/prepRun.sh was unable to locate a ddPCR results file." | tee -a "$logf"
	echo "This csv file is required for run processing." | tee -a "$logf"
	echo "Please see $logf for more information." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 0 (run file prep)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


status=$(./perl/1_compileRun.pl "$indir" | tee -a "$logf")

if [[ "$status" != "0" ]]
then
	echo "perl/1_compileRun.pl reported a fatal error." | tee -a "$logf"
	echo "patchr will now abort." | tee -a "$logf"
	echo "Please see $logf for more information." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 1 (run data compilation)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "Backing up watchdb LATEST/ directory to LATEST_BK/." | tee -a "$logf"
cp -r "data/watchdb/LATEST/" "data/watchdb/LATEST_BK/"

tables=("abatch" "archive" "assay" "cbatch" "concentration" "control" "ebatch" "extraction" "rbatch")
for i in ${!tables[@]}
do
	table=${tables[$i]}
	updatef="$indir/updates/update.$table.txt"
	dbinf="data/watchdb/LATEST/watchdb.$table.txt"
	dboutf="data/watchdb/$TODAY/watchdb.$table.txt"

	echo "Appending $updatef to watchdb.$table.txt." | tee -a "$logf"
	awk '(NR>1)' "$updatef" | cat - "$dbinf" > "$dboutf"
	cp "$dboutf" "$dbinf"
done

#status=$(./perl/2_validateWaTCH.pl | tee -a "$logf")


status=$(./perl/3_feedDashboard.pl | tee -a "$logf")

if [[ "$status" != "0" ]]
then
	echo "perl/3_feedDashboard.pl reported a fatal error." | tee -a "$logf"
	echo "patchr will now abort." | tee -a "$logf"
	echo "Please see $logf for more information." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 3 (Updating dashboard feed)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


status=$(./perl/4_feedNWSS.pl | tee -a "$logf")

if [[ "$status" != "0" ]]
then
	echo "perl/4_feedNWSS.pl reported a fatal error." | tee -a "$logf"
	echo "patchr will now abort." | tee -a "$logf"
	echo "Please see $logf for more information." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 4 (Updating NWSS feed)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


END=$(date + '%Y-%m-%d_%H-%M')

echo "" | tee -a "$logf"
echo "$END" | tee -a "$logf"
echo "Finished patchr.sh" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"

exit 0

