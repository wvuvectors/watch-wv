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
			OUTDIR=$OPTARG
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
	echo "FATAL : no input directory provided (-i)."
	exit 1
fi

if [ -z "$OUTDIR" ]; then
	echo "FATAL : no output directory provided (-o)."
	exit 1
fi

# Get the current date and time
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")
RUNID="$START"


# Make the output directory, if necessary.
if [ ! -d "$OUTDIR/" ]
then
	mkdir "$OUTDIR/"
	echo "This folder was created on $UPDAY." > "$OUTDIR/README.txt"
fi

# Locate and (if necessary) create the archive file that lists all the batches that have 
# been completed. This lives in the UPDATES dir.
ARCHFILE="$OUTDIR/UPDATES/completed_batches.txt"
if [ ! -f "$ARCHFILE/" ]
then
	echo "batch_id" > "$ARCHFILE"
fi


# Make a directory to hold this run, in the UPDATES sub-folder.
UPDIR="$OUTDIR/UPDATES/$RUNID"
if [ ! -d "$OUTDIR/UPDATES" ]
then
	mkdir "$OUTDIR/UPDATES/"
fi
mkdir "$UPDIR"
mkdir "$UPDIR/BATCH_FILES"


# Write all output to log file in the update dir.
logf="$UPDIR/sewer.$RUNID.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

touch "$logf"
echo "#############################################" | tee -a "$logf"
echo "SEWER: Surveillance & Exploration of Wastewater to inform Epidemiological Response." | tee -a "$logf"
echo "Initiated from sewer.sh." | tee -a "$logf"
echo "$START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "Input data dir: $INDIR" | tee -a "$logf"
echo "Output run dir: $UPDIR/$RUNID" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Searching for unprocessed batch files in $INDIR." | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 1_queryBatches.R $INDIR $UPDIR $ARCHFILE." | tee -a "$logf"
num2proc=$(Rscript R/1_queryBatches.R "$INDIR" "$UPDIR" "$ARCHFILE")

status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
echo "" | tee -a "$logf"

if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "1_queryBatches.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $INDIR $UPDIR $ARCHFILE" | tee -a "$logf"
	echo "sewer aborted during phase 1 (batch identification)." | tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	echo "" | tee -a "$logf"
	exit 1
fi

if [ $num2proc == 0 ]
then
	echo "1_queryBatches.R found no new batches in $INDIR." | tee -a "$logf"
	echo "As a result, there is nothing for sewer to do at this time." | tee -a "$logf"
	echo "******" | tee -a "$logf"
	echo "All done. sewer run of $START will now exit, having done nothing."
	echo "******" | tee -a "$logf"
	exit 1
fi

echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"



echo "1_queryBatches.R identified $num2proc new batch files and copied these files to: " | tee -a "$logf"
echo "   $UPDIR/" | tee -a "$logf"
echo "The original batch files in $INDIR will remain untouched." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "Preparing to update WaTCH from these batch files now." | tee -a "$logf"
echo "First I'll add the sample data file from AssetTiger to the update folder." | tee -a "$logf"
if [ -f "$INDIR/AssetTagReport.csv" ]
then
	# Copy it from the input folder to the update folder.
	cp "$INDIR/AssetTagReport.csv" "$UPDIR/BATCH_FILES/Samples.csv"
	echo "Sample	AT	$INDIR/AssetTagReport.csv	Samples.csv" >> "$UPDIR/update_files.txt"
else
	echo "!!!!!!!!" | tee -a "$logf"
	echo "Unable to locate a Samples file from AssetTiger in the input dir:" | tee -a "$logf"
	echo "$INDIR/AssetTagReport.csv" | tee -a "$logf"
	echo "sewer aborted during phase 1 (no AssetTiger input file)." | tee -a "$logf"
	echo "It is possible the file was renamed or moved. "| tee -a "$logf"
	echo "It is also possible that I do not have permission to access this file. "| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"


echo "Now I'll fix any known file issues in these input files:" | tee -a "$logf"
echo "   Windows line endings are replaced with UNIX line endings." | tee -a "$logf"
echo "   Special characters, including micron, are replaced with standard alphanumerics." | tee -a "$logf"
while read -r line
do
	IFS="\t" read -r -a uprow <<< "$line"

	f="$UPDIR/BATCH_FILES/${uprow[4]}"

	if [[ "$f" == *.csv ]]
	then
		perl -pi -e 's/\r$//' "$f"
		sed -i '' -e '$a\' "$f"
		sed $'1s/\xef\xbb\xbf//' < "$f" > "$f.tmp"
		mv "$f.tmp" "$f"
		
		is_qia=$(head -n 1 "$f" | grep -c "sep=,")
		if [[ "$is_qia" == "1" ]];then
			sed -i '1d' "$f"
		fi
		
		is_result=$(head -n 1 "$f" | grep -c "Well,")
		if [[ "$is_result" == "1" ]];then
			perl -pi -e 's/µ/u/i' "$f"
		fi
	fi
done < "$UPDIR/update_files.txt"

echo "Done preparing the batch files." | tee -a "$logf"
echo "******" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "Compiling the update files." | tee -a "$logf"
echo "Note: each batch is compiled separately to allow more control over the process." | tee -a "$logf"
echo "" | tee -a "$logf"



echo "******" | tee -a "$logf"
echo "Running 2a_compileCBATCHES.R." | tee -a "$logf"

c_count=$(Rscript R/2a_compileCBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2a_compileCBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 2a (CBATCH update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "" | tee -a "$logf"
if [[ $c_count < 0 ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "FATAL: At least one row of concentration data is missing required data." | tee -a "$logf"
	echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
	echo "FATAL: Check the validation_key column in $UPDIR/update.concentration.txt and $UPDIR/update.cbatch.txt." | tee -a "$logf"
	echo "FATAL: Fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Running 2b_compileEBATCHES.R." | tee -a "$logf"
e_count=$(Rscript R/2b_compileEBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2b_compileEBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 2b (EBATCH update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "" | tee -a "$logf"
if [[ $e_count < 0 ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "FATAL: At least one row of extraction data is missing required data." | tee -a "$logf"
	echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
	echo "FATAL: Check the validation_key column in $UPDIR/update.extraction.txt and $UPDIR/update.ebatch.txt." | tee -a "$logf"
	echo "FATAL: Fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Running 2c_compileQBATCHES.R." | tee -a "$logf"
q_count=$(Rscript R/2c_compileQBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2c_compileQBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 2c (QBATCH update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "" | tee -a "$logf"
if [[ $q_count < 0 ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "FATAL: At least one row of quantification data is missing required data." | tee -a "$logf"
	echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
	echo "FATAL: Check the validation_key columns in $UPDIR/update.quant_run.txt, $UPDIR/update.quant_control.txt, $UPDIR/update.quant_block.txt, and $UPDIR/update.qbatch.txt." | tee -a "$logf"
	echo "FATAL: Fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Running 2d_compileQRESULTS.R." | tee -a "$logf"
r_count=$(Rscript R/2d_compileQRESULTS.R "$UPDIR")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2d_compileQRESULTS.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 2d (QRESULT update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "" | tee -a "$logf"
if [[ $r_count < 0 ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "FATAL: At least one row of quantification data is missing results." | tee -a "$logf"
	echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
	echo "FATAL: Check the validation_key columns in $UPDIR/update.quant_run_merged.txt and $UPDIR/update.quant_control_merged.txt." | tee -a "$logf"
	echo "FATAL: Fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Running 3_addSAMPLES.R." | tee -a "$logf"
s_count=$(Rscript R/3_addSAMPLES.R "$UPDIR")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "3_addSAMPLES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 3 (SAMPLE addition)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "" | tee -a "$logf"
if [[ $s_count < 0 ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "FATAL: At least one row of sample data is missing required information." | tee -a "$logf"
	echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
	echo "FATAL: Check the validation_key columns in $UPDIR/update.sample.txt." | tee -a "$logf"
	echo "FATAL: Fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Finished compiling the update tables!" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "******" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "This update contains:" | tee -a "$logf"
echo "   $s_count new samples;" | tee -a "$logf"
echo "   $c_count new concentrations;" | tee -a "$logf"
echo "   $e_count new extractions;" | tee -a "$logf"
echo "   $q_count new quantification runs;" | tee -a "$logf"
echo "   $r_count new assay results (including controls)." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Now I'm going to validate the update in $UPDIR against the existing data in $OUTDIR/latest/." | tee -a "$logf"
echo "This checks for common issues such as spurious sample ids, missing concentration or extraction data, and more." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Running 4_validateUpdate.R." | tee -a "$logf"
v_count=$(Rscript R/4_validateUpdate.R "$UPDIR" "$OUTDIR/latest")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "4_validateUpdate.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR $OUTDIR/latest" | tee -a "$logf"
	echo "sewer aborted during phase 4 (update validation)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs, then delete $UPDIR, address the errors, " | tee -a "$logf"
	echo "and run sewer again from the start."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"




echo "" | tee -a "$logf"
echo "Adding a validate update is a simple matter of concatenating the update to the existing tables in $OUTDIR/latest/." | tee -a "$logf"
echo "So that's what I'll do now." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "First, I'll back up $OUTDIR/latest/ to $OUTDIR/latest_bk/ just in case we need to roll all the way back." | tee -a "$logf"
cp -r "$OUTDIR/latest/" "$OUTDIR/latest_bk/"
PREVDATE=$(head -n 1 "$OUTDIR/latest/README.txt")
echo "$PREVDATE" > "$OUTDIR/latest_bk/README.txt"
echo "" >> "$OUTDIR/latest_bk/README.txt"
echo "#" >> "$OUTDIR/latest_bk/README.txt"
echo "This folder contains backups of the watchdb tables from before $UPDAY." >> "$OUTDIR/latest_bk/README.txt"
echo "#" >> "$OUTDIR/latest_bk/README.txt"
echo "Done." | tee -a "$logf"



echo "" | tee -a "$logf"
echo "Now I'll apply the validated update from $UPDIR to $OUTDIR/latest/." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 5_applyUpdate.R." | tee -a "$logf"
u_count=$(Rscript R/5_applyUpdate.R "$UPDIR" "$OUTDIR/latest")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "5_applyUpdate.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR $OUTDIR/latest" | tee -a "$logf"
	echo "sewer aborted during phase 5 (applying the update)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs and address the errors. " | tee -a "$logf"
	echo "You may need to rollback $OUTDIR/latest to $OUTDIR/latest_bk, and/or delete $UPDIR and start over."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"



echo "" | tee -a "$logf"
echo "It seems like the update was applied successfully!" | tee -a "$logf"
echo "Before we finish, though, I need to regenerate some processor-intensive data." | tee -a "$logf"
echo "This is primarily aimed at optimizing data loading into the dashboard." | tee -a "$logf"

echo "First I'll generate a new result table for $OUTDIR/latest." | tee -a "$logf"
echo "This table is an amalgam of data tables and consequently will contain a certain level of redundancy." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "******" | tee -a "$logf"
echo "Running 6_generateResults.R." | tee -a "$logf"
result_count=$(Rscript R/6_generateResults.R echo "$OUTDIR/latest")
status="${PIPESTATUS[0]}"
echo "" | tee -a "$logf"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "6_generateResults.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $OUTDIR/latest" | tee -a "$logf"
	echo "sewer aborted during phase 6 (generating the result table)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs, then delete $UPDIR and $OUTDIR/incremental/$START, " | tee -a "$logf"
	echo "address the errors, and run sewer again from the start."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


cp "$OUTDIR/incremental/$START/watchdb.completed_batches.txt" "$OUTDIR/latest/watchdb.completed_batches.txt"

echo "" | tee -a "$logf"
echo "Updating the README file." | tee -a "$logf"
echo "$UPDAY" > "$OUTDIR/latest/README.txt"
echo "" >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"
echo "This folder contains the most recent watchdb tables." >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"


echo "" | tee -a "$logf"
echo "File copy finished." | tee -a "$logf"
echo "The most recent version of the watchdb can be found in $OUTDIR/latest." | tee -a "$logf"
echo "The version of the watchdb from just before this update was applied can be found in $OUTDIR/latest_bk." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "All done! sewer run of $START will now exit." | tee -a "$logf"
echo "" | tee -a "$logf"

