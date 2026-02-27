#! /bin/bash

while getopts ":hi:o:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 1
			;;
		i)
			INDIR=$OPTARG
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


if [ -z "$INDIR" ]; then
	echo "WARN : no input directory provided (-i). Using the default directory:"
	echo "WARN : /Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/DRISCOLL_LAB/2 PROJECTS/WaTCH/TESTING_LAB/DATA_PCR/"
	INDIR="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/DRISCOLL_LAB/2 PROJECTS/WaTCH/TESTING_LAB/DATA_PCR/"
fi

if [ -z "$DBDIR" ]; then
	echo "WARN : no data output directory provided (-o). Using the default directory:"
	echo "WARN : watch-wv/patchr/data/"
	DBDIR="data"
fi

# These paths point to the Marshall Univ lab data files.
SHARED_DIR='/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY'
MU_F='mu_nwss.LATEST.csv'


# Get the current date and time
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")

# Ensure that the database directories exist for writing output.
if [ ! -d "$DBDIR/" ]
then
	mkdir "$DBDIR/"
fi
if [ ! -d "$DBDIR/latest" ]
then
	mkdir "$DBDIR/latest"
	echo "This folder was created on $UPDAY." > "$DBDIR/latest/README.txt"
fi
if [ ! -d "$DBDIR/latest_bk" ]
then
	mkdir "$DBDIR/latest_bk"
	echo "This folder was created on $UPDAY." > "$DBDIR/latest_bk/README.txt"
fi


# Make the update directory.
if [ ! -d "$DBDIR/updates" ]
then
	mkdir "$DBDIR/updates/"
fi
if [ ! -d "$DBDIR/updates/$START" ]
then
	mkdir "$DBDIR/updates/$START"
fi
UPDIR="$DBDIR/updates/$START"
mkdir "$UPDIR/batches"


# Write all output to log file in the update dir.
logf="$UPDIR/patchr.$START.log"
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
echo "Input data dir: $INDIR" | tee -a "$logf"
echo "Output data dir: $DBDIR" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Searching for unprocessed batch files in $INDIR." | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 1_queryBatches.R $UPDIR $DBDIR/latest $INDIR." | tee -a "$logf"
num2proc=$(echo "$UPDIR" "$DBDIR/latest" "$INDIR" | ./R/1_queryBatches.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "1_queryBatches.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR $DBDIR/latest $INDIR" | tee -a "$logf"
	echo "patchr aborted during phase 1 (batch identification)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "" | tee -a "$logf"
	exit 1
fi

if [ $num2proc == 0 ]
then
	echo "1_queryBatches.R found no new batches in $INDIR." | tee -a "$logf"
	echo "As a result, there is nothing for patchr to do at this time." | tee -a "$logf"
	echo "******" | tee -a "$logf"
	echo "All done. patchr run of $START will now exit, having done nothing."
	echo "******" | tee -a "$logf"
	exit 1
fi

echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "1_queryBatches.R identified $num2proc new batch files and copied these files to: " | tee -a "$logf"
echo "$UPDIR/batches/." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "Preparing to update WaTCH from these batch files now." | tee -a "$logf"
echo "First I'll add the sample data file from AssetTiger to the update folder." | tee -a "$logf"
if [ -f "$INDIR/0 SAMPLES/AssetTagReport.csv" ]
then
	# Copy it from the input folder to the update dir.
	cp "$INDIR/0 SAMPLES/AssetTagReport.csv" "$UPDIR/batches/Samples.csv"
	echo "AT	ATDB	$UPDIR/batches/Samples.csv	Samples.csv" >> "$UPDIR/update.batch_files.txt"
else
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Unable to locate a samples file from AssetTiger in the input dir:" | tee -a "$logf"
	echo "$INDIR/0 SAMPLES/AssetTagReport.csv" | tee -a "$logf"
	echo "patchr aborted during phase 1 (no AssetTiger input file)." | tee -a "$logf"
	echo "It is possible the file was renamed or moved. "| tee -a "$logf"
	echo "It is also possible that I do not have permission to access this file. "| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"

echo "Now I'll add the MU data file to the update folder." | tee -a "$logf"
if [ -f "$SHARED_DIR/$MU_F" ]
then
	# Copy it from the shared folder to the update dir.
	cp "$SHARED_DIR/$MU_F" "$UPDIR/batches/$MU_F"
	echo "MU	MUID	$SHARED_DIR/$MU_F	$MU_F" >> "$UPDIR/update.batch_files.txt"
else
	echo "Unable to locate a data file from the MU testing lab:" | tee -a "$logf"
	echo "$SHARED_DIR/$MU_F" | tee -a "$logf"
	echo "I will use a previous file instead:" | tee -a "$logf"
	echo "$DBDIR/latest/$MU_F" | tee -a "$logf"
	cp "$DBDIR/latest/$MU_F" "$UPDIR/batches/$MU_F"
	echo "MU	MUID	$DBDIR/latest/$MU_F	$MU_F" >> "$UPDIR/update.batch_files.txt"
fi
echo "Done." | tee -a "$logf"

echo "Now I'll fix any known file issues in these input files, like Windows line endings." | tee -a "$logf"
while read -r line
do
	IFS="\t" read -r -a uprow <<< "$line"

	f="$UPDIR/batches/${uprow[3]}"

	if [[ "$f" == *.csv ]]
	then
		perl -pi -e 's/\r$//' "$f"
		sed -i '' -e '$a\' "$f"
		sed $'1s/\xef\xbb\xbf//' < "$f" > "$f.tmp"
		mv "$f.tmp" "$f"

		is_ddpcr=$(head -n 1 "$f" | grep -c "Well,")

		if [[ "$is_ddpcr" == "1" ]];then
			perl -pi -e 's/µ/u/i' "$f"
		fi
	fi
done < "$UPDIR/update.batch_files.txt"
echo "Done preparing the batch files." | tee -a "$logf"
echo "******" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "Compiling the update files." | tee -a "$logf"
echo "Note: each batch is compiled separately to allow more control over the process." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 2a_compileCPLATES.R." | tee -a "$logf"

c_count=$(echo "$UPDIR" | ./R/2a_compileCPLATES.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2a_compileCPLATES.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "patchr aborted during phase 2a (CPLATE update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"

echo "" | tee -a "$logf"

echo "Running 2b_compileEPLATES.R." | tee -a "$logf"
e_count=$(echo "$UPDIR" | ./R/2b_compileEPLATES.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2b_compileEPLATES.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "patchr aborted during phase 2b (EPLATE update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"

echo "" | tee -a "$logf"

echo "Running 2c_compileAPLATES.R." | tee -a "$logf"
a_count=$(echo "$UPDIR" | ./R/2c_compileAPLATES.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2c_compileAPLATES.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "patchr aborted during phase 2c (APLATE update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"

echo "" | tee -a "$logf"

echo "Running 2d_compileSAMPLES.R." | tee -a "$logf"
s_count=$(echo "$UPDIR" | ./R/2d_compileSAMPLES.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2d_compileSAMPLES.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "patchr aborted during phase 2d (SAMPLE update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run patchr again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "Finished compiling the update tables!" | tee -a "$logf"
echo "This update contains:" | tee -a "$logf"
echo "   $s_count new samples;" | tee -a "$logf"
echo "   $c_count new concentrations;" | tee -a "$logf"
echo "   $e_count new extractions;" | tee -a "$logf"
echo "   $a_count new assays." | tee -a "$logf"


echo "" | tee -a "$logf"
echo "Ok, I'm about to start modifying the watchdb files in $DBDIR/latest/." | tee -a "$logf"
echo "First, let's back up $DBDIR/latest/ to $DBDIR/latest_bk/ just in case." | tee -a "$logf"
cp -r "$DBDIR/latest/" "$DBDIR/latest_bk/"

PREVDATE=$(head -n 1 "$DBDIR/latest/README.txt")
echo "$PREVDATE" > "$DBDIR/latest_bk/README.txt"
echo "" >> "$DBDIR/latest_bk/README.txt"
echo "#" >> "$DBDIR/latest_bk/README.txt"
echo "This folder contains backups of the previous watchdb tables." >> "$DBDIR/latest_bk/README.txt"
echo "#" >> "$DBDIR/latest_bk/README.txt"
echo "Done." | tee -a "$logf"

echo "" | tee -a "$logf"
echo "Now I'll create a new watchdb folder in $DBDIR/incremental/$START/." | tee -a "$logf"
echo "It will hold the latest watchdb after the update from $UPDIR has been applied." | tee -a "$logf"

if [ -d "$DBDIR/incremental/$START/" ]
then
	echo "WARN: An old version of $DBDIR/incremental/$START/ already exists, so I'm removing the old directory first." | tee -a "$logf"
	rm -r "$DBDIR/incremental/$START/"
fi
mkdir "$DBDIR/incremental/$START/"

echo "" | tee -a "$logf"
echo "Copying current db files from latest to $DBDIR/incremental/$START/." | tee -a "$logf"
dbfiles=$(ls -1 "$DBDIR/latest/*.txt")
for i in ${!dbfiles[@]}
do
	f=${dbfiles[$i]}
	if [ -f "$f" ]
	then
		cp "$f" "$DBDIR/incremental/$START/"
	fi
done
echo "Done preparing the current db for the update." | tee -a "$logf"
echo "******" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "Now let's validate the update in $UPDIR against the existing data in $DBDIR/incremental/$START/." | tee -a "$logf"
echo "This checks for common issues such as spurious sample ids, missing concentration or extraction data, and more." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 3_validateUpdate.R." | tee -a "$logf"
u_count=$(echo "$UPDIR" "$DBDIR/incremental/$START" | ./R/3_validateUpdate.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "3_validateUpdate.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR $DBDIR/incremental/$START" | tee -a "$logf"
	echo "patchr aborted during phase 3 (update validation)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs, then delete $UPDIR and $DBDIR/incremental/$START, " | tee -a "$logf"
	echo "address the errors, and run patchr again from the start."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "Now I'll apply the validated update from $UPDIR to $DBDIR/incremental/$START/." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 4_applyUpdate.R." | tee -a "$logf"
u_count=$(echo "$UPDIR" "$DBDIR/incremental/$START" | ./R/4_applyUpdate.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "4_applyUpdate.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR $DBDIR/incremental/$START" | tee -a "$logf"
	echo "patchr aborted during phase 4 (applying the update to the incremental version)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs, then delete $UPDIR and $DBDIR/incremental/$START, " | tee -a "$logf"
	echo "address the errors, and run patchr again from the start."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"



echo "" | tee -a "$logf"
echo "It looks like the update was successfully applied to $DBDIR/incremental/$START." | tee -a "$logf"
echo "Before we copy it to the production folder, though, I need to regenerate some processor-intensive data." | tee -a "$logf"
echo "This is primarily aimed at optimizing data loading into the dashboard." | tee -a "$logf"

echo "First I'll generate the result table using the updated tables in $DBDIR/incremental/$START." | tee -a "$logf"
echo "This table is an amalgam of data tables and consequently will contain a certain level of redundancy." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 5_generateResults.R." | tee -a "$logf"
result_count=$(echo "$DBDIR/incremental/$START" | ./R/5_generateResults.R)
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "5_generateResults.R exited with error code $status and caused patchr to abort." | tee -a "$logf"
	echo "Arguments: $DBDIR/incremental/$START" | tee -a "$logf"
	echo "patchr aborted during phase 5 (generating the results table)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs, then delete $UPDIR and $DBDIR/incremental/$START, " | tee -a "$logf"
	echo "address the errors, and run patchr again from the start."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "Everything appears to be ok, so on to the final step." | tee -a "$logf"
echo "I am now copying the watchdb tables from $DBDIR/incremental/$START/ to $DBDIR/latest/." | tee -a "$logf"
echo "This overwrites the previous set of data, which is why we save it for last." | tee -a "$logf"
echo "Never fear! The most recent previous version is still available in $DBDIR/latest_bk/." | tee -a "$logf"
echo "Also, older sets are collected in the $DBDIR/incremental/ folder." | tee -a "$logf"
echo "******" | tee -a "$logf"


new_dbfiles=$(ls -1 "$DBDIR/incremental/$START/*.txt")
for i in ${!new_dbfiles[@]}
do
	f=${new_dbfiles[$i]}
	if [ -f "$f" ];then
		cp "$f" "$DBDIR/latest/watchdb.$table.txt"
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

