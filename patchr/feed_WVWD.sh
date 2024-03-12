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

cp "$RSDIR/watchdb.all_tables.xlsx" "../dashboard/data/watchdb.all_tables.xlsx"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of watchdb.all_tables.xlsx from $RSDIR to dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

cp "$SVDIR/seqrdb.txt" "../dashboard/data/seqrdb.txt"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of seqrdb.txt from $SVDIR to dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "" | tee -a "$logf"
echo "Updating the README file." | tee -a "$logf"
echo "$UPDAY" > "../dashboard/data/README.txt"
echo "" >> "../dashboard/data/README.txt"
echo "#" >> "../dashboard/data/README.txt"
echo "This folder contains the most recent dashboard raw data." >> "../dashboard/data/README.txt"
echo "#" >> "../dashboard/data/README.txt"



echo "All done! feed_WVWD will now exit."
echo ""

