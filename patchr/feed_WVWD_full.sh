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


WVWD_MU_LTST="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/READY/mu_dashboard.LATEST.tsv"
WVWD_MU_ARCH="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/ARCHIVED/MU/mu_wvwd.$START.tsv"

if [ -f "$WVWD_MU_LTST" ]
then

	echo "Archiving original WVWD file $WVWD_MU_LTST from Marshall lab."
	cp "$WVWD_MU_LTST" "$WVWD_MU_ARCH"
	echo "Configuring dashboard file $WVWD_MU_LTST from Marshall lab."
	perl -pi -e 's/\r$//' "$WVWD_MU_LTST"
	sed $'1s/\xef\xbb\xbf//' < "$WVWD_MU_LTST" > "$WVWD_MU_LTST.tmp"
	mv "$WVWD_MU_LTST.tmp" "$WVWD_MU_LTST"
#	sed -i '' -e '$a\' "$WVWD_MU_LTST"
	echo "Done."
	
	echo "******" | tee -a "$logf"
	echo "Running 7b_mergeWVWD.pl." | tee -a "$logf"
	echo "******" | tee -a "$logf"

	./perl/7b_mergeWVWD.pl "$wvwdF" "$WVWD_MU_LTST" < "resources/WVWD_fields.txt" > "$wvwdM" | tee -a "$logf"
	status="${PIPESTATUS[0]}"
	echo "" | tee -a "$logf"

	if [[ "$status" != "0" ]]
	then
		echo "7b_mergeWVWD.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		echo "patchr_feed aborted during phase 2 (WVWD merge)." | tee -a "$logf"
		echo "Delete the file $wvwdM, if it exists. "| tee -a "$logf"
		echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
	fi
else
	echo "I could not find $WVWD_MU_LTST from Marshall lab, so there is nothing to merge." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		echo "patchr_feed aborted during phase 2 (WVWD merge)." | tee -a "$logf"
		echo "NO MERGE WAS PERFORMED. "| tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
fi


echo ""
echo "Everything seems to have gone ok so I'm going to copy the WVWD files to the WaTCH shared directory."
cp "$wvwdM" "$wvwdP"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of $wvwdM to $wvwdP exited with error code $status." | tee -a "$logf"
	echo "This is not fatal but the file copy DID NOT HAPPEN." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi

cp "$wvwdF" "$wvwdA"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of $wvwdF to $wvwdA exited with error code $status." | tee -a "$logf"
	echo "This is not fatal but the file copy DID NOT HAPPEN." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi

rm "$WVWD_MU_LTST"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Removal (rm) of $WVWD_MU_LTST exited with error code $status." | tee -a "$logf"
	echo "This is not fatal but the file was not removed." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi

echo ""
echo "Done."
echo ""


echo "All done! feed_WVWD will now exit."
echo ""

