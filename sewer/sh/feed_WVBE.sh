#! /bin/bash

DBDIR="../dashboard/data/latest"
RSDIR="resources"

# Get the current date and time
START=$(date "+%F_%H-%M")
TODAY=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")

OUTDIR="../breatheeasy/data"
OUTFILE="watch_data.breatheeasy.$START.txt"

# Write all output to log file
logf="logs/feeds/WVBE/feed_WVBE.$START.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "Initiated feed_WVBE.sh" | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Generating WV BreatheEasy (WVBE) feed for WV BPH." | tee -a "$logf"

# Date as the 1st saturday of collection week
# Target
# Abundance
# Abundance (pop norm)
# Location (county name)
# Location (FIPS code)
# Trend level
# Abundance level

echo "******" | tee -a "$logf"
echo "Running 1_feed_WVBE.R." | tee -a "$logf"
echo "******" | tee -a "$logf"

R/1_feed_WVBE.R > "$OUTDIR/incremental/$OUTFILE" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "1_feed_WVBE.R exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during phase 1 (WVBE)." | tee -a "$logf"
	echo "Delete the file $OUTF. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr_feed \"WVBE\" again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


cp "$OUTDIR/incremental/$OUTFILE" "$OUTDIR/latest/watch_data.breatheeasy.latest.txt"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of $OUTFILE from $OUTDIR/incremental/ to $OUTDIR/latest/ exited with error code $status." | tee -a "$logf"
	echo "The incremental file in $OUTDIR/incremental/$OUTFILE should be okay but for some reason I can't copy it to the latest/ directory." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "" | tee -a "$logf"
echo "Everything seems to have gone fine so I'm updating the README file." | tee -a "$logf"

touch "$OUTDIR/latest/README.txt"
echo "$UPDAY" > "$OUTDIR/latest/README.txt"
echo "" >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"
echo "This folder contains the most recent data for the breatheeasy Web site." >> "$OUTDIR/latest/README.txt"
echo "https://breatheeasy.wv.gov/." >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"


echo "All done! feed_WVBE will now exit."
echo ""

