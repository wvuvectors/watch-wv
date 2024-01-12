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

wvwdF="data/extref/wvwd/wvu_wvwd.$START.txt"
wvwdM="data/extref/wvwd/merged_wvwd.$START.txt"

wvwdP="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/READY/merged_wvwd.LATEST.txt"
wvwdA="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/ARCHIVED/WVU/wvu_wvwd.$START.txt"

echo "******" | tee -a "$logf"
echo "Running 7_feedWVWD.pl $DBDIR." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/7_feedWVWD.pl "$DBDIR" > "$wvwdF" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "7_feedWVWD.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during phase 1 (WVWD)." | tee -a "$logf"
	echo "Delete the file $wvwdF. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "All done! feed_WVWD will now exit."
echo ""

