#! /bin/bash

DBDIR="data/latest"
RSDIR="resources"
SVDIR="../seqr/data/latest"

WVD_DIR="../dashboard"

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

ALL_RESOURCE_F="$RSDIR/watchdb.all_tables.xlsx"	# This is simply copied from the patchr directory.
WVU_RESULTS_F="$DBDIR/watchdb.result.txt"	# This starts as the watchdb.result.txt file from patchr.
MU_RESULTS_F="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/mu_nwss.LATEST.csv"


echo "Ensuring all input files exist and are readable before starting." | tee -a "$logf"
if [ ! -f "$ALL_RESOURCE_F" ];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "$ALL_RESOURCE_F does not exist. This is a fatal error." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
if [ ! -f "$WVU_RESULTS_F" ];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "$WVU_RESULTS_F does not exist. This is a fatal error." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
if [ ! -f "$MU_RESULTS_F" ];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "$MU_RESULTS_F does not exist. This is a fatal error." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "All input files appear to exist, so we're moving on." | tee -a "$logf"
echo "Copying the latest WaTCH resource tables to the dashboard directory." | tee -a "$logf"

cp "$ALL_RESOURCE_F" "../dashboard/data/watchdb.all_tables.xlsx"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of $ALL_RESOURCE_F to ../dashboard/data exited with error code $status." | tee -a "$logf"
	echo "This is fatal." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "******" | tee -a "$logf"
echo "Running 8_MU2WVDash.pl." | tee -a "$logf"
echo "This subroutine takes the raw MU data file as input and converts it into watchdb format." | tee -a "$logf"
echo "The result is written into the dashboard data directory as mu.result.txt." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/8_MU2WVDash.pl < "$MU_RESULTS_F" > "$WVD_DIR/data/mu.result.txt" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "8_MU2WVDash.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during MU data conversion." | tee -a "$logf"
	echo "Delete the file mu.result.txt $WVD_DIR/data/, if it exists. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


# Identifying outliers and extreme values in the input data.
#
echo "******" | tee -a "$logf"
echo "Running 9_tagOutliers.R." | tee -a "$logf"
echo "This subroutine reads data from $WVU_RESULTS_F and $MU_RESULTS_F and adds three columns:" | tee -a "$logf"
echo "A numeric z-score for each entry, and logical columns is.outlier and is.extreme." | tee -a "$logf"
echo "The modified files are written to the dashboard data directory; the original files are unchanged." | tee -a "$logf"
echo "******" | tee -a "$logf"

echo "$WVU_RESULTS_F $WVD_DIR/data/mu.result.txt $ALL_RESOURCE_F" | R/9_tagOutliers.R | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "9_tagOutliers.R exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during the outlier tagging phase." | tee -a "$logf"
	echo "Check the logs, fix the error(s), and run patchr_feed again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

# We can remove the temporary MU table now that we've successfully added outlier data.
rm "$WVD_DIR/data/mu.result.txt"


# Creating a routine surveillance summary table to hold most recent risk, abundance, and trend data
# by region (state, county, and facility). Must be pre-calculated because too computationally
# expensive to do them all when the dashboard initializes.
#
echo "******" | tee -a "$logf"
echo "Running 10_buildRSStable.R." | tee -a "$logf"
echo "******" | tee -a "$logf"

R/10_buildRSStable.R | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "10_buildRSStable.R exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during routine surveillance summary table construction." | tee -a "$logf"
	echo "Delete the file wvdash.rsstable.txt in ../dashboard/data/, if it exists. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
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





echo "" | tee -a "$logf"
echo "Updating the README file." | tee -a "$logf"

patchrUp=$(head -n 1 "$DBDIR/README.txt")
seqrUP=$(head -n 1 "$SVDIR/README.txt")

echo "$UPDAY" > "../dashboard/data/README.txt"
echo "" >> "../dashboard/data/README.txt"
echo "#" >> "../dashboard/data/README.txt"
echo "This folder contains the most recent dashboard data." >> "../dashboard/data/README.txt"
echo "Abundance data was last updated on $patchrUp." >> "../dashboard/data/README.txt"
#echo "MU data was last updated on $DBUP_MU." >> "../dashboard/data/README.txt"
echo "SEQR data was last updated on $seqrUP." >> "../dashboard/data/README.txt"
echo "#" >> "../dashboard/data/README.txt"


echo "All done! feed_WVDash will now exit."
echo ""

