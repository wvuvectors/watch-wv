#! /bin/bash

indir="$1"
WD=$(pwd)


# Get the current date and time
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")

PATCHR_DBDIR="../patchr/data/latest"
SEQR_DBDIR="data"


# Write all output to log file
logf="logs/seqr.$START.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "Initiated seqr.sh" | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"



echo "Compiling update from run directory $indir." | tee -a "$logf"
upfile="$SEQR_DBDIR/updates/seqr_update.$START.txt"

echo "******" | tee -a "$logf"
echo "Running 1sv_compileUpdate.pl $indir > $upfile." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/1sv_compileUpdate.pl "$indir" > "$upfile"
status="${PIPESTATUS[0]}"

if [[ "$status" != "0" ]]
then
	echo "" | tee -a "$logf"
	echo "1sv_compileUpdate.pl exited with error code $status and caused seqr to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "seqr aborted during phase 1 (update compilation)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "" | tee -a "$logf"
echo "Now I'll create a new seqr folder in $SEQR_DBDIR/incremental/$START/." | tee -a "$logf"
echo "It will hold the latest seqrdb after updating with the data from $upfile." | tee -a "$logf"
echo "This is an incremental backup file." | tee -a "$logf"

if [ -d "$SEQR_DBDIR/incremental/$START/" ]
then
	echo "WARN: An old version of $SEQR_DBDIR/incremental/$START/ already exists, so I'm removing the old directory first." | tee -a "$logf"
	rm -r "$SEQR_DBDIR/incremental/$START/"
fi
mkdir "$SEQR_DBDIR/incremental/$START/"

echo "" | tee -a "$logf"
echo "Applying update $upfile to seqr database in $SEQR_DBDIR/incremental/$START/." | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 2sv_appendUpdate.pl -w $PATCHR_DBDIR -i $upfile -o $SEQR_DBDIR/incremental/$START/." | tee -a "$logf"
echo "******" | tee -a "$logf"


./perl/2sv_appendUpdate.pl -w "$PATCHR_DBDIR" -i "$upfile" -o "$SEQR_DBDIR/incremental/$START" | tee -a "$logf"
status="${PIPESTATUS[0]}"

if [[ "$status" != "0" ]]
then
	echo "" | tee -a "$logf"
	echo "2sv_appendUpdate.pl exited with error code $status and caused seqr to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "seqr aborted during phase 2 (append update to seqrdb)." | tee -a "$logf"
	echo "NO UPDATE WAS APPLIED. "| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "" | tee -a "$logf"
echo "Ok, I'm about to start modifying the seqrdb files in $SEQR_DBDIR/latest/." | tee -a "$logf"
echo "First, let's back up $SEQR_DBDIR/latest/ to $SEQR_DBDIR/latest_bk/." | tee -a "$logf"
cp -r "$SEQR_DBDIR/latest/" "$SEQR_DBDIR/latest_bk/"

status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Copy of $SEQR_DBDIR/latest/ to $SEQR_DBDIR/latest/ exited with error code $status." | tee -a "$logf"
	echo "This is fatal: I refuse to proceed without a backup of the latest database." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

PREVDATE=$(head -n 1 "$SEQR_DBDIR/latest/README.txt")
echo "$PREVDATE" > "$SEQR_DBDIR/latest_bk/README.txt"
echo "" >> "$SEQR_DBDIR/latest_bk/README.txt"
echo "#" >> "$SEQR_DBDIR/latest_bk/README.txt"
echo "This folder contains a copy of the previous seqr tables." >> "$SEQR_DBDIR/latest_bk/README.txt"
echo "#" >> "$SEQR_DBDIR/latest_bk/README.txt"


echo "" | tee -a "$logf"
echo "Backup complete!." | tee -a "$logf"
echo "Now I'm going to copy from $SEQR_DBDIR/incremental/$START to $SEQR_DBDIR/latest/." | tee -a "$logf"

tables=("seqrdb")
for i in ${!tables[@]}
do
	table=${tables[$i]}
	if [ -f "$SEQR_DBDIR/incremental/$START/$table.txt" ]
	then
		cp "$SEQR_DBDIR/incremental/$START/$table.txt" "$SEQR_DBDIR/latest/$table.txt"
		status="${PIPESTATUS[0]}"
		if [[ "$status" != "0" ]]
		then
			echo "!!!!!!!!" | tee -a "$logf"
			echo "Copy of $SEQR_DBDIR/incremental/$START/$table.txt to $SEQR_DBDIR/latest/$table.txt exited with error code $status." | tee -a "$logf"
			echo "This is not fatal, but the 'latest' folder DOES NOT CONTAIN THE LATEST FILE." | tee -a "$logf"
			echo "You may want to manually copy $SEQR_DBDIR/incremental/$START/$table.txt to $SEQR_DBDIR/latest/$table.txt" | tee -a "$logf"
			echo "but I strongly recommend that you figure out why this file copy failed." | tee -a "$logf"
			echo "!!!!!!!!" | tee -a "$logf"
		fi
	fi
done

echo "" | tee -a "$logf"
echo "Updating the README file." | tee -a "$logf"
echo "$UPDAY" > "$SEQR_DBDIR/latest/README.txt"
echo "" >> "$SEQR_DBDIR/latest/README.txt"
echo "#" >> "$SEQR_DBDIR/latest/README.txt"
echo "This folder contains the most recent seqrdb tables." >> "$SEQR_DBDIR/latest/README.txt"
echo "#" >> "$SEQR_DBDIR/latest/README.txt"


echo "" | tee -a "$logf"
echo "File copy finished."
echo "The most recent version of the seqrdb can be found in two places:" | tee -a "$logf"
echo "    $SEQR_DBDIR/incremental/$START/" | tee -a "$logf"
echo "    and" | tee -a "$logf"
echo "    $SEQR_DBDIR/latest/" | tee -a "$logf"
echo "$SEQR_DBDIR/latest_bk/ contains the data from immediately before this update was applied." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "All done! seqr will now exit."
echo "" | tee -a "$logf"


exit 0



