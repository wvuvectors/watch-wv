#! /bin/bash

while getopts ":hi:o:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 1
			;;
		i)
			indir=$OPTARG
			;;
		o)
			DBDIR=$OPTARG
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

if [ -z "$indir" ]; then
	echo "WARN : no input directory provided (-i). Using the default directory:"
	echo "WARN : /Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/DRISCOLL_LAB/2 PROJECTS/WaTCH/TESTING_LAB/DATA_PCR/"
	indir="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/DRISCOLL_LAB/2 PROJECTS/WaTCH/TESTING_LAB/DATA_PCR/"
fi

if [ -z "$DBDIR" ]; then
	echo "WARN : no data output directory provided (-o). Using the default directory:"
	echo "WARN : watch-wv/patchr/data/"
	DBDIR="data"
fi

WD=$(pwd)
 
# Get the current date and time
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")

if [ ! -d "$DBDIR/" ]
then
	mkdir "$DBDIR/"
fi

if [ ! -d "$DBDIR/updates" ]
then
	mkdir "$DBDIR/updates/"
fi

if [ ! -d "$DBDIR/latest" ]
then
	mkdir "$DBDIR/latest"
	echo "Folder created  $UPDAY" > "$DBDIR/latest/README.txt"
fi

if [ ! -d "$DBDIR/latest_bk" ]
then
	mkdir "$DBDIR/latest_bk"
	echo "Folder created  $UPDAY" > "$DBDIR/latest_bk/README.txt"
fi

update_dir="$DBDIR/updates/2025-07-20_15-33/"

# Write all output to log file
logf="logs/patchr/patchr.$START.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "Initiated patchr.sh" | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "Input data dir: $indir" | tee -a "$logf"
echo "Output data dir: $DBDIR" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 2_compileUpdate.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"

./perl/2_compileUpdate.pl "$update_dir" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "2_compileUpdate.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 2 (update compilation)." | tee -a "$logf"
	echo "Delete the folder $update_dir. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi




echo "" | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "Running 3_validateUpdate.pl." | tee -a "$logf"
echo "******" | tee -a "$logf"


./perl/3_validateUpdate.pl -i "$update_dir" -d "$DBDIR/latest" | tee -a "$logf"
status="${PIPESTATUS[0]}"
#echo "Status of validation: $status" | tee -a "$logf"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "3_validateUpdate.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Most likely this script identifed duplicate IDs or collisions (IDs present in both the existing db and the update). Check the following files for more info:" | tee -a "$logf"
	echo "    $update_dir/_collisions.txt" | tee -a "$logf"
	echo "    $update_dir/_rundups.txt" | tee -a "$logf"
	echo "    $update_dir/_watchdups.txt" | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 3 (update validation)." | tee -a "$logf"
	echo "Delete the folder $update_dir. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "" | tee -a "$logf"
echo "Ok, I'm about to start modifying the watchdb files in $DBDIR/latest/." | tee -a "$logf"
echo "First, let's back up $DBDIR/latest/ to $DBDIR/latest_bk/." | tee -a "$logf"
cp -r "$DBDIR/latest/" "$DBDIR/latest_bk/"

PREVDATE=$(head -n 1 "$DBDIR/latest/README.txt")
echo "$PREVDATE" > "$DBDIR/latest_bk/README.txt"
echo "" >> "$DBDIR/latest_bk/README.txt"
echo "#" >> "$DBDIR/latest_bk/README.txt"
echo "This folder contains backups of the previous watchdb tables." >> "$DBDIR/latest_bk/README.txt"
echo "#" >> "$DBDIR/latest_bk/README.txt"


echo "" | tee -a "$logf"
echo "Now I'll create a new watchdb folder in $DBDIR/incremental/$START/." | tee -a "$logf"
echo "It will hold the latest watchdb after update with data from $update_dir." | tee -a "$logf"

if [ -d "$DBDIR/incremental/$START/" ]
then
	echo "WARN: An old version of $DBDIR/incremental/$START/ already exists, so I'm removing the old directory first." | tee -a "$logf"
	rm -r "$DBDIR/incremental/$START/"
fi
mkdir "$DBDIR/incremental/$START/"
echo "" | tee -a "$logf"

tables=("abatch" "archive" "assay" "cbatch" "concentration" "ebatch" "extraction" "rbatch" "result" "sample")
for i in ${!tables[@]}
do
	table=${tables[$i]}
	if [ -f "$DBDIR/latest/watchdb.$table.txt" ]
	then
		cp "$DBDIR/latest/watchdb.$table.txt" "$DBDIR/incremental/$START/watchdb.$table.txt"
	fi
done



echo "" | tee -a "$logf"
echo "Now I'm merging the watchdb tables from $DBDIR/latest/ with the update tables in $update_dir." | tee -a "$logf"
echo "These merged tables will be put in the folder that I just created, $DBDIR/incremental/$START/." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 4_appendUpdate.pl $DBDIR/incremental/$START $update_dir." | tee -a "$logf"
echo "******" | tee -a "$logf"


./perl/4_appendUpdate.pl "$DBDIR/incremental/$START" "$update_dir" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "4_appendUpdate.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	#echo "Removing $DBDIR/watchdb/$START/" | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 4 (appending update)." | tee -a "$logf"
	echo "I STRONGLY recommend deleting $DBDIR/incremental/$START/ after exploring this error." | tee -a "$logf"
	echo "Then delete the folder $update_dir, fix the error(s), and run patchr again. "| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi


echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Now I will generate the RESULT table using the updated watchdb tables in $DBDIR/incremental/$START." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 5_calculateResults.pl $DBDIR/incremental/$START." | tee -a "$logf"
echo "******" | tee -a "$logf"


./perl/5_calculateResults.pl "$DBDIR/incremental/$START" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "5_calculateResults.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	#echo "Removing $DBDIR/watchdb/$START/" | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 5 (calculating results)." | tee -a "$logf"
	echo "Check the file $DBDIR/incremental/$START/_result.rejections.txt for details." | tee -a "$logf"
	echo "I STRONGLY recommend deleting $DBDIR/incremental/$START/ after exploring this error." | tee -a "$logf"
	echo "Then delete the folder $update_dir, fix the error(s), and run patchr again. "| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi



echo "" | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "Running 6_generateBatchList.pl $DBDIR/incremental/$START > $DBDIR/incremental/$START/watchdb.completed_batches.txt." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"

./perl/6_generateBatchList.pl "$DBDIR/incremental/$START" > "$DBDIR/incremental/$START/watchdb.completed_batches.txt" | tee -a "$logf"
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "6_generateBatchList.pl exited with error code $status and caused patchr to abort." | tee -a "$logf"
	#echo "Removing $DBDIR/watchdb/$START/" | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "patchr aborted during phase 6 (generating list of completed batches)." | tee -a "$logf"
	echo "Sometimes this error is due to a missing or misnamed file." | tee -a "$logf"
	echo "If the fix is not simple, I STRONGLY recommend deleting $DBDIR/incremental/$START/." | tee -a "$logf"
	echo "Then fix the error(s), delete the folder $update_dir, and run patchr again. "| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi



echo "" | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "Everything appears to be ok, so on to the final step." | tee -a "$logf"
echo "I am now copying the watchdb tables from $DBDIR/incremental/$START/ to $DBDIR/latest/, overwriting the old set." | tee -a "$logf"
echo "Never fear! The previous version is still readily available in $DBDIR/latest_bk/, and older sets in the $DBDIR/incremental/ folder." | tee -a "$logf"
echo "******" | tee -a "$logf"

for i in ${!tables[@]}
do
	table=${tables[$i]}
	if [ -f "$DBDIR/incremental/$START/watchdb.$table.txt" ]
	then
		cp "$DBDIR/incremental/$START/watchdb.$table.txt" "$DBDIR/latest/watchdb.$table.txt"
	fi
done
cp "$DBDIR/incremental/$START/watchdb.completed_batches.txt" "$DBDIR/latest/watchdb.completed_batches.txt"

echo "" | tee -a "$logf"
echo "Updating the README file." | tee -a "$logf"
echo "$UPDAY" > "$DBDIR/latest/README.txt"
echo "" >> "$DBDIR/latest/README.txt"
echo "#" >> "$DBDIR/latest/README.txt"
echo "This folder contains the most recent watchdb tables." >> "$DBDIR/latest/README.txt"
echo "#" >> "$DBDIR/latest/README.txt"


echo "" | tee -a "$logf"
echo "File copy finished." | tee -a "$logf"
echo "The most recent version of the watchdb can be found in two places:" | tee -a "$logf"
echo "    $DBDIR/incremental/$START/" | tee -a "$logf"
echo "    and" | tee -a "$logf"
echo "    $DBDIR/latest/" | tee -a "$logf"
echo "$DBDIR/latest_bk/ contains the data from immediately before this update was applied." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "All done! patchr run of $START will now exit." | tee -a "$logf"
echo "" | tee -a "$logf"

