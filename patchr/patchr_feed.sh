#! /bin/bash

WD=$(pwd)

DBDIR="data/latest"

# Get the current date and time
START=$(date "+%F_%H-%M")
TODAY=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")

# Write all output to log file
logf="logs/patchr_feed.$START.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "Initiated patchr_feed.sh" | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Generating NWSS feed for WVU." | tee -a "$logf"

nwssF="data/extref/nwss/wvu_nwss.$START.csv"
nwssM="data/extref/nwss/merged_nwss.$START.csv"

nwssP="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/merged_nwss.LATEST.csv"
nwssA="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/ARCHIVED/WVU/wvu_nwss.$START.csv"

echo "******" | tee -a "$logf"
echo "Running 7_feedNWSS.pl $DBDIR." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/7_feedNWSS.pl "$DBDIR" > "$nwssF" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "7_feedNWSS.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr_feed aborted during phase 1 (NWSS)." | tee -a "$logf"
	echo "Delete the file $nwssF. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


NWSS_MU_LTST="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/mu_nwss.LATEST.csv"
NWSS_MU_ARCH="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/ARCHIVED/MU/mu_nwss.$START.csv"

if [ -f "$NWSS_MU_LTST" ]
then

	echo "Archiving original NWSS file $NWSS_MU_LTST from Marshall lab."
	cp "$NWSS_MU_LTST" "$NWSS_MU_ARCH"
	echo "Configuring NWSS file $NWSS_MU_LTST from Marshall lab."
	perl -pi -e 's/\r$//' "$NWSS_MU_LTST"
	sed $'1s/\xef\xbb\xbf//' < "$NWSS_MU_LTST" > "$NWSS_MU_LTST.tmp"
	mv "$NWSS_MU_LTST.tmp" "$NWSS_MU_LTST"
#	sed -i '' -e '$a\' "$NWSS_MU_LTST"
	echo "Done."
	
	echo "******" | tee -a "$logf"
	echo "Running 7b_mergeNWSS.pl." | tee -a "$logf"
	echo "******" | tee -a "$logf"

	./perl/7b_mergeNWSS.pl "$nwssF" "$NWSS_MU_LTST" < "resources/NWSS_fields.txt" > "$nwssM" | tee -a "$logf"
	status="${PIPESTATUS[0]}"
	echo "" | tee -a "$logf"

	if [[ "$status" != "0" ]]
	then
		echo "7b_mergeNWSS.pl exited with error code $status and caused patchr_feed to abort." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		echo "patchr_feed aborted during phase 3 (NWSS merge)." | tee -a "$logf"
		echo "Delete the file $nwssM, if it exists. "| tee -a "$logf"
		echo "Then fix the error(s) and run patchr_feed again."| tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
	fi
else
	echo "I could not find $NWSS_MU_LTST from Marshall lab, so there is nothing to merge." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		echo "patchr_feed aborted during phase 3 (NWSS merge)." | tee -a "$logf"
		echo "NO MERGE WAS PERFORMED. "| tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
fi


echo ""
echo "Everything seems to have gone ok so I'm going to copy the NWSS files to the WaTCH shared directory."
cp "$nwssM" "$nwssP"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of $nwssM to $nwssP exited with error code $status." | tee -a "$logf"
	echo "This is not fatal but the file copy DID NOT HAPPEN." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi

cp "$nwssF" "$nwssA"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of $nwssF to $nwssA exited with error code $status." | tee -a "$logf"
	echo "This is not fatal but the file copy DID NOT HAPPEN." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi

rm "$NWSS_MU_LTST"
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Removal (rm) of $NWSS_MU_LTST exited with error code $status." | tee -a "$logf"
	echo "This is not fatal but the file was not removed." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi

echo ""
echo "Done."
echo ""


echo "All done! patchr_feed will now exit."
echo ""

