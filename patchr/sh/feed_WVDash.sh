#! /bin/bash

DBDIR="data/latest"
RSDIR="resources"
SVDIR="../seqr/data/latest"


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
echo "Initiated feed_WVDash.sh" | tee -a "$logf"
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
	echo "Copy of watchdb.result.txt from $DBDIR to ../dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

cp "$DBDIR/watchdb.sample.txt" "../dashboard/data/watchdb.sample.txt"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of watchdb.sample.txt from $DBDIR to ../dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


cp "$RSDIR/watchdb.all_tables.xlsx" "../dashboard/data/watchdb.all_tables.xlsx"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of watchdb.all_tables.xlsx from $RSDIR to ../dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


cp "$SVDIR/seqrdb.txt" "../dashboard/data/seqrdb.txt"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of seqrdb.txt from $SVDIR to ../dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi



NWSS_F="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/merged_nwss.LATEST.csv"

if [ -f "$NWSS_F" ]
then
	echo "******" | tee -a "$logf"
	echo "Running 8a_MU2WVDash.pl." | tee -a "$logf"
	echo "******" | tee -a "$logf"

	./perl/8a_MU2WVDash.pl "$NWSS_F" | tee -a "$logf"
	status="${PIPESTATUS[0]}"
	echo "" | tee -a "$logf"

	if [[ "$status" != "0" ]]
	then
		echo "8a_MU2WVDash.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		echo "patchr_feed aborted during phase 2 (MU data merge)." | tee -a "$logf"
		echo "Delete the files mu.result.txt and mu.sample.txt in ../dashboard/data/, if they exist. "| tee -a "$logf"
		echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
	fi
else
		echo "There is no merged NWSS file available." | tee -a "$logf"
		echo "The Marshall data file in ../dashboard/data, if present, will be untouched." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		echo "MU DATA WAS NOT UPDATED. "| tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
fi


# Creating a dashboard-specific table to hold most recent risk, abundance, and trend data
# by region (state, county, or facility). Must be pre-calculated because too computationally
# expensive to do them all when the dashboard initializes.
#
echo "******" | tee -a "$logf"
echo "Running 9_buildAlertTable.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/9_buildAlertTable.pl > "../dashboard/data/wvdash.alerts.txt" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "9_buildAlertTable.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during phase 3 (WVDash optimization)." | tee -a "$logf"
	echo "Delete the file wvdash.alerts.txt in ../dashboard/data/, if it exists. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi



echo "" | tee -a "$logf"
echo "Updating the README file." | tee -a "$logf"

patchrUp=$(head -n 1 "$DBDIR/README.txt")
seqrUP=$(head -n 1 "$SVDIR/README.txt")

echo "$UPDAY" > "../dashboard/data/README.txt"
echo "" >> "../dashboard/data/README.txt"
echo "#" >> "../dashboard/data/README.txt"
echo "This folder contains the most recent dashboard data." >> "../dashboard/data/README.txt"
echo "WVU data was last updated on $patchrUp." >> "../dashboard/data/README.txt"
#echo "MU data was last updated on $DBUP_MU." >> "../dashboard/data/README.txt"
echo "SEQR data was last updated on $seqrUP." >> "../dashboard/data/README.txt"
echo "#" >> "../dashboard/data/README.txt"


echo "All done! feed_WVDash will now exit."
echo ""

