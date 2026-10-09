#! /bin/bash

while getopts ":hi:o:m:" opt; do
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
	echo "FATAL : no input run directory provided (-i)."
	exit 1
fi

if [ ! -d "$INDIR" ]; then
	echo "FATAL : the input directory does not exist or is not readable."
	exit 1
fi


if [ -z "$OUTDIR" ]; then
	echo "FATAL : no output directory provided (-o)."
	exit 1
fi


# Get the current date and time for logging purposes.
#
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")
TS=$(date +%s)

# Get the run ID, which is the name of input folder.
#
RUNID=$(basename $INDIR)
INCID="$RUNID"
#INCID="$TS.$RUNID"

# Make the output directory and all required subdirs, as necessary.
#
if [ ! -d "$OUTDIR/" ]
then
	mkdir "$OUTDIR/"
	echo "This folder was created on $UPDAY." > "$OUTDIR/README.txt"
fi

if [ ! -d "$OUTDIR/latest" ]
then
	mkdir "$OUTDIR/latest"
fi
if [ ! -d "$OUTDIR/incremental" ]
then
	mkdir "$OUTDIR/incremental"
fi
if [ ! -d "$OUTDIR/updates" ]
then
	mkdir "$OUTDIR/updates"
fi


# Look for RUNID in any of the existing batch files. If it exists, do not proceed. Prompt the user to 
# make sure the current RUNID is unique.
#
FOUND=$(grep -wi "$RUNID" "$OUTDIR/latest/watchdb.qbatch.txt")
if [ $FOUND == 0 ]
then
	echo "FATAL : Every run of the laboratory workflow is assigned an ID using the name of the folder that contains the raw files."
	echo "FATAL : It looks like a run labeled $RUNID has already been processed."
	echo "FATAL : This can occur if you forgot to change the name of the run directory."
	echo "FATAL : It might also occur if you have multiple runs on the same calendar day."
	echo "FATAL : In this case, add a suffix to the name of the run directory (e.g., 2026-06-22b) and run sewer again."
	exit 1
fi


# Make a directory to hold the update that will be generated.
#
UPDIR="$OUTDIR/updates/$RUNID"
if [ ! -d "$OUTDIR/updates" ]
then
	mkdir "$OUTDIR/updates/"
fi
mkdir "$UPDIR"
mkdir "$UPDIR/run_files"


# Write all output to log file in the update dir.
#
logf="$UPDIR/sewer.log"
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
echo "Run ID: $RUNID" | tee -a "$logf"
echo "Input dir: $INDIR" | tee -a "$logf"
echo "Update dir: $UPDIR" | tee -a "$logf"
echo "" | tee -a "$logf"


# Move all of the input files to the update dir and "fix" them.
#
echo "Copying all of the input files for run $RUNID to the run_files/ directory in $UPDIR." | tee -a "$logf"
echo "The original run files in $INDIR will remain untouched." | tee -a "$logf"
for fin in "$INDIR"/*
do
	f=$(basename "$fin")
	fup="$UPDIR/run_files/$f"
	cp "$fin" "$fup"
	echo "Copied $fin to $fup"
done

echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"

# echo "Also copying the sample metadata file to the run_files/ directory in $UPDIR." | tee -a "$logf"
# cp "$METADATA_F" "$UPDIR/run_files/all_samples.csv"
# echo "Copied $METADATA_F to $UPDIR/run_files/all_samples.csv"
# echo "Done." | tee -a "$logf"
# echo "" | tee -a "$logf"

echo "Now I'll fix any known file issues in these input files:" | tee -a "$logf"
echo "   Windows line endings are replaced with UNIX line endings in csv or tab-delim text files." | tee -a "$logf"
echo "   Special characters, including micron, are replaced with standard alphanumerics." | tee -a "$logf"
for f in "$UPDIR/run_files"/*
do

	if [[ "$f" == *.csv ]]
	then
		perl -pi -e 's/\r$//' "$f"
		sed -i '' -e '$a\' "$f"
		sed $'1s/\xef\xbb\xbf//' < "$f" > "$f.tmp"
		mv "$f.tmp" "$f"
		
# 		is_sample=$(head -n 1 "$f" | grep -c "Asset Tag ID")
# 		if [[ "$is_sample" == "1" ]];then
# 			mv "$f" "Samples.csv"
# 		fi

		is_qia=$(head -n 1 "$f" | grep -c "sep=,")
		if [[ "$is_qia" == "1" ]];then
			sed -i '1d' "$f"
		fi
		
		is_qx=$(head -n 1 "$f" | grep -c "Well,")
		if [[ "$is_qx" == "1" ]];then
			perl -pi -e 's/µ/u/i' "$f"
		fi
	fi
done
echo "Done preparing the run files." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"

echo "Compiling the update files." | tee -a "$logf"
echo "Note: each batch is compiled separately to allow more control over the process." | tee -a "$logf"
echo "" | tee -a "$logf"



echo "******" | tee -a "$logf"
echo "Running 1_compileCBATCHES.R." | tee -a "$logf"
c_count=$(Rscript R/1_compileCBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "1_compileCBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 1 (CBATCH update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
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


echo "******" | tee -a "$logf"
echo "Running 2_compileEBATCHES.R." | tee -a "$logf"
e_count=$(Rscript R/2_compileEBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2_compileEBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 2 (EBATCH update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
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


echo "******" | tee -a "$logf"
echo "Running 3_compileQBATCHES.R." | tee -a "$logf"
q_count=$(Rscript R/3_compileQBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "3_compileQBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 3 (QBATCH update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
if [[ $q_count < 0 ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "FATAL: At least one row of quantification data is missing required data." | tee -a "$logf"
	echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
	echo "FATAL: Check the validation_key columns in $UPDIR/update.qassay.txt, $UPDIR/update.qcontrol.txt, $UPDIR/update.qblock.txt, and $UPDIR/update.qbatch.txt." | tee -a "$logf"
	echo "FATAL: Fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 4_compileQRESULTS.R." | tee -a "$logf"
r_count=$(Rscript R/4_compileQRESULTS.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "4_compileQRESULTS.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 4 (QRESULT update compilation)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
if [[ $r_count < 0 ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "FATAL: At least one row of quantification data is missing results." | tee -a "$logf"
	echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
	echo "FATAL: Check the validation_key columns in $UPDIR/update.qassay_merged.txt and $UPDIR/update.qcontrol_merged.txt." | tee -a "$logf"
	echo "FATAL: Fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "Finishing the merge of assay and result data." | tee -a "$logf"
mv "$UPDIR/update.qassay_merged.txt" "$UPDIR/update.qassay.txt"
mv "$UPDIR/update.qcontrol_merged.txt" "$UPDIR/update.qcontrol.txt"
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 5_addSAMPLES.R." | tee -a "$logf"
s_count=$(Rscript R/5_addSAMPLES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "5_addSAMPLES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR" | tee -a "$logf"
	echo "sewer aborted during phase 5 (SAMPLE addition)." | tee -a "$logf"
	echo "Delete the folder $UPDIR. "| tee -a "$logf"
	echo "Then fix the error(s) and run sewer again."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
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


echo "******" | tee -a "$logf"
echo "Finished compiling all update tables!" | tee -a "$logf"
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
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Now I'll start the process of applying the update from $UPDIR to $OUTDIR/latest/." | tee -a "$logf"
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


echo "Now I'll append the update tables from $UPDIR to the existing data in $OUTDIR/latest." | tee -a "$logf"
echo "The merge will first be written to $OUTDIR/incremental/$RUNID." | tee -a "$logf"
echo "If successful, the new files will be copied from $OUTDIR/incremental/$RUNID to $OUTDIR/latest." | tee -a "$logf"
echo "" | tee -a "$logf"

mkdir "$OUTDIR/incremental/$RUNID"

echo "******" | tee -a "$logf"
echo "Running 6_applyUpdate.R." | tee -a "$logf"
dup_count=$(Rscript R/6_applyUpdate.R "$UPDIR" "$OUTDIR/latest" "$OUTDIR/incremental/$RUNID")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "6_applyUpdate.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $UPDIR $OUTDIR/latest $OUTDIR/incremental/$RUNID" | tee -a "$logf"
	echo "sewer aborted during phase 6 (applying the update)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs and address the errors. " | tee -a "$logf"
	echo "You may need to delete $UPDIR and/or $OUTDIR/incremental/$RUNID, and start over."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "It seems like the merge was successful!" | tee -a "$logf"
echo "Before we finish, though, I need to pre-calculate some processor-intensive data." | tee -a "$logf"
echo "This is primarily aimed at optimizing realtime dashboard performance and generates a file called watchdb.result.txt." | tee -a "$logf"
echo "The table in this file is an amalgam of data tables and consequently will contain a certain level of redundancy." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Running 7_generateResults.R." | tee -a "$logf"
result_count=$(Rscript R/7_generateResults.R "$OUTDIR/incremental/$RUNID")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "7_generateResults.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "Arguments: $OUTDIR/incremental/$RUNID" | tee -a "$logf"
	echo "sewer aborted during phase 7 (generating the result table)." | tee -a "$logf"
	echo "The safest route to recovery is to check the run logs and address the errors. " | tee -a "$logf"
	echo "You may need to delete $UPDIR and/or $OUTDIR/incremental/$RUNID, and start over."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "Final copy of files from $OUTDIR/incremental/$RUNID to $OUTDIR/latest." | tee -a "$logf"
cp "$OUTDIR/incremental/$RUNID/watchdb."*.txt "$OUTDIR/latest/"
echo "Done." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "Updating the README file." | tee -a "$logf"
echo "$UPDAY" > "$OUTDIR/latest/README.txt"
echo "" >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"
echo "This folder contains the most recent watchdb tables." >> "$OUTDIR/latest/README.txt"
echo "#" >> "$OUTDIR/latest/README.txt"
echo "" | tee -a "$logf"


echo "******" | tee -a "$logf"
echo "The most recent version of the watchdb can be found in $OUTDIR/latest." | tee -a "$logf"
echo "The version of the watchdb from just before this update was applied can be found in $OUTDIR/latest_bk." | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "sewer run of $UPDAY has left the chat." | tee -a "$logf"
echo "" | tee -a "$logf"

