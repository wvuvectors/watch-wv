#! /bin/bash


# Get the current date and time
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")


while getopts ":hi:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 1
			;;
		i)
			inDIR=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

if [ -z "$inDIR" ]; then
	echo "Missing a valid input directory of demix files (-i)!"
	echo "This is a required parameter with no default, so I'm forced to quit. Sorry."
	exit 1
fi

logFILE="logs/seqr.$START.log"

WATCHDB_DIR="../patchr/data/latest"
SEQR_DIR="data"

upDIR="$SEQR_DIR/updates/$START/"

# Write all output to log file
if [ -f "$logFILE" ]
then
	rm "$logFILE"
fi

touch "$logFILE"
echo "#############################################" | tee -a "$logFILE"
echo "Initiated seqr.sh" | tee -a "$logFILE"
echo "$START" | tee -a "$logFILE"
echo "" | tee -a "$logFILE"
echo "See $logFILE for warnings, errors, and other important information." | tee -a "$logFILE"
echo "" | tee -a "$logFILE"


if [ -d "$upDIR" ]
then
	echo "WARN: An old version of $upDIR/ already exists!" | tee -a "$logFILE"
	echo "WARN: I am removing the old directory." | tee -a "$logFILE"
	rm -r "$upDIR/"
fi
mkdir "$upDIR/"




echo "Compiling update from input directory $inDIR." | tee -a "$logFILE"

echo "******" | tee -a "$logFILE"
echo "Running seqr_1_compileUpdate.pl -i $inDIR -o $upDIR." | tee -a "$logFILE"
echo "******" | tee -a "$logFILE"

./perl/seqr_1_compileUpdate.pl -i "$inDIR" -o "$upDIR"
status="${PIPESTATUS[0]}"

if [[ "$status" != "0" ]]
then
	echo "" | tee -a "$logFILE"
	echo "seqr_1_compileUpdate.pl exited with error code $status and caused seqr to abort." | tee -a "$logFILE"
	echo "!!!!!!!!" | tee -a "$logFILE"
	echo "seqr aborted during phase 1 (update compilation)." | tee -a "$logFILE"
	echo "!!!!!!!!" | tee -a "$logFILE"
	exit 1
fi


echo "" | tee -a "$logFILE"
echo "Now I'll create a new folder in $SEQR_DIR/incremental/$START/." | tee -a "$logFILE"
echo "It will hold the latest seqr data after updating with the data from $upDIR." | tee -a "$logFILE"
echo "This is an incremental backup file." | tee -a "$logFILE"

if [ -d "$SEQR_DIR/incremental/$START/" ]
then
	echo "WARN: An old version of $SEQR_DIR/incremental/$START/ already exists!" | tee -a "$logFILE"
	echo "WARN: I am removing the old directory." | tee -a "$logFILE"
	rm -r "$SEQR_DIR/incremental/$START/"
fi
mkdir "$SEQR_DIR/incremental/$START/"

echo "" | tee -a "$logFILE"
echo "Appending update from $upDIR to the latest seqr data from $SEQR_DIR/latest/." | tee -a "$logFILE"
echo "This will create an incremental version, stored in $SEQR_DIR/incremental/$START/." | tee -a "$logFILE"

echo "******" | tee -a "$logFILE"
echo "Running seqr_2_appendUpdate.pl -w $WATCHDB_DIR -d $SEQR_DIR/latest -u $upDIR -o $SEQR_DIR/incremental/$START." | tee -a "$logFILE"
echo "******" | tee -a "$logFILE"


./perl/seqr_2_appendUpdate.pl -w "$WATCHDB_DIR" -d "$SEQR_DIR/latest" -u "$upDIR" -o "$SEQR_DIR/incremental/$START" | tee -a "$logFILE"
status="${PIPESTATUS[0]}"

if [[ "$status" != "0" ]]
then
	echo "" | tee -a "$logFILE"
	echo "seqr_2_appendUpdate.pl exited with error code $status and caused seqr to abort." | tee -a "$logFILE"
	echo "!!!!!!!!" | tee -a "$logFILE"
	echo "seqr aborted during phase 2 (append update to seqrdb)." | tee -a "$logFILE"
	echo "NO UPDATE WAS APPLIED. "| tee -a "$logFILE"
	echo "!!!!!!!!" | tee -a "$logFILE"
	exit 1
fi


echo "" | tee -a "$logFILE"
echo "Ok, I'm about to start modifying the seqrdb files in $SEQR_DIR/latest/." | tee -a "$logFILE"
echo "First, let's back up $SEQR_DIR/latest/ to $SEQR_DIR/latest_bk/." | tee -a "$logFILE"
cp -r "$SEQR_DIR/latest/" "$SEQR_DIR/latest_bk/"

status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logFILE"
	echo "Copy of $SEQR_DIR/latest/ to $SEQR_DIR/latest/ exited with error code $status." | tee -a "$logFILE"
	echo "This is fatal: I refuse to proceed without a backup of the latest database." | tee -a "$logFILE"
	echo "!!!!!!!!" | tee -a "$logFILE"
	exit 1
fi

PREVDATE=$(head -n 1 "$SEQR_DIR/latest/README.txt")
echo "$PREVDATE" > "$SEQR_DIR/latest_bk/README.txt"
echo "" >> "$SEQR_DIR/latest_bk/README.txt"
echo "#" >> "$SEQR_DIR/latest_bk/README.txt"
echo "This folder contains a copy of the previous seqr tables." >> "$SEQR_DIR/latest_bk/README.txt"
echo "#" >> "$SEQR_DIR/latest_bk/README.txt"


echo "" | tee -a "$logFILE"
echo "Backup complete!." | tee -a "$logFILE"
echo "Now I'm going to copy from $SEQR_DIR/incremental/$START to $SEQR_DIR/latest/." | tee -a "$logFILE"

tables=("seqr")
for i in ${!tables[@]}
do
	table=${tables[$i]}
	if [ -f "$SEQR_DIR/incremental/$START/watchdb.$table.txt" ]
	then
		cp "$SEQR_DIR/incremental/$START/watchdb.$table.txt" "$SEQR_DIR/latest/watchdb.$table.txt"
		status="${PIPESTATUS[0]}"
		if [[ "$status" != "0" ]]
		then
			echo "!!!!!!!!" | tee -a "$logFILE"
			echo "Copy of $SEQR_DIR/incremental/$START/watchdb.$table.txt to $SEQR_DIR/latest/watchdb.$table.txt exited with error code $status." | tee -a "$logFILE"
			echo "This is not fatal, but the 'latest' folder DOES NOT CONTAIN THE LATEST FILE." | tee -a "$logFILE"
			echo "You may want to manually copy $SEQR_DIR/incremental/$START/watchdb.$table.txt to $SEQR_DIR/latest/watchdb.$table.txt" | tee -a "$logFILE"
			echo "but I strongly recommend that you figure out why this file copy failed." | tee -a "$logFILE"
			echo "!!!!!!!!" | tee -a "$logFILE"
		fi
	fi
done

echo "" | tee -a "$logFILE"
echo "Updating the README file." | tee -a "$logFILE"
echo "$UPDAY" > "$SEQR_DIR/latest/README.txt"
echo "" >> "$SEQR_DIR/latest/README.txt"
echo "#" >> "$SEQR_DIR/latest/README.txt"
echo "This folder contains the most recent seqr tables." >> "$SEQR_DIR/latest/README.txt"
echo "#" >> "$SEQR_DIR/latest/README.txt"


echo "" | tee -a "$logFILE"
echo "File copy finished." | tee -a "$logFILE"
echo "The most recent version of the seqr data can be found in two places:" | tee -a "$logFILE"
echo "    $SEQR_DIR/incremental/$START/" | tee -a "$logFILE"
echo "    and" | tee -a "$logFILE"
echo "    $SEQR_DIR/latest/" | tee -a "$logFILE"
echo "$SEQR_DIR/latest_bk/ contains the data from immediately before this update was applied." | tee -a "$logFILE"
echo "" | tee -a "$logFILE"



if [ -d "reports/$START/" ]
then
	echo "WARN: An old version of reports/$START/ already exists!" | tee -a "$logFILE"
	echo "WARN: I am removing the old directory." | tee -a "$logFILE"
	rm -r "reports/$START/"
fi
mkdir "reports/$START/"

REPORT_DIR="reports/$START"
echo "" | tee -a "$logFILE"
echo "Writing report to $REPORT_DIR."  | tee -a "$logFILE"

echo "******" | tee -a "$logFILE"
echo "Running seqr_3_generateReport.pl -w $WATCHDB_DIR -d $SEQR_DIR/latest -o $REPORT_DIR." | tee -a "$logFILE"
echo "******" | tee -a "$logFILE"

./perl/seqr_3_generateReport.pl -w "$WATCHDB_DIR" -d "$SEQR_DIR/latest" -o "$REPORT_DIR"
status="${PIPESTATUS[0]}"

if [[ "$status" != "0" ]]
then
	echo "" | tee -a "$logFILE"
	echo "seqr_3_generateReport.pl exited with error code $status and caused seqr to abort." | tee -a "$logFILE"
	echo "!!!!!!!!" | tee -a "$logFILE"
	echo "seqr aborted during phase 3 (report generation)." | tee -a "$logFILE"
	echo "The data update was successful, but a report could not be generated." | tee -a "$logFILE"
	echo "!!!!!!!!" | tee -a "$logFILE"
	exit 1
fi


echo "All done! seqr will now exit."
echo "" | tee -a "$logFILE"


exit 0



