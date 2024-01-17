#! /bin/bash

DBDIR="data/latest"

# Get the current date and time
START=$(date "+%F_%H-%M")
TODAY=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")

# Write all output to log file
logf="logs/feeds/WVWD/feed_WVWD.$START.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "Initiated feed_WVWD.sh" | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Generating WV WaTCH Dashboard (WVWD) feed for WVU." | tee -a "$logf"

cp "$DBDIR/watchdb.result.txt" "../dashboard/data/watchdb.result.txt"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of watchdb.result.txt from $DBDIR to dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

cp "$DBDIR/watchdb.sample.txt" "../dashboard/data/watchdb.sample.txt"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of watchdb.sample.txt from $DBDIR to dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "******" | tee -a "$logf"
echo "Running 7_feedWVWD.pl '../dashboard/data'." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/7_feedWVWD.pl "../dashboard/data" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "7_feedWVWD.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during phase 3 (resource file generation)." | tee -a "$logf"
	echo "Fix the error(s) and run patchr_feed again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "All done! feed_WVWD will now exit."
echo ""

