#! /bin/bash

OUTDIR="../data"
QC=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INDIR="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -o|--output)
      OUTDIR="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -q|--noqc)
      QC=0
      shift 1 # Shift past the flag
      ;;
    -l|--logfile)
      logf="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -h|--help)
      cat << EOF
SEWER make update script.
EOF
		exit 0
		;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done


if [[ -z "$INDIR" ]];then
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : no input directory provided (-i). This is a FATAL error."
	echo "Please use the -i argument to provide a valid path to the workflow run folder "
	echo "that you would like to process. This is required regardless of any other settings."
	echo "Use -h or --help for more information about how to run this script. (Eg, $0 --help)."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi

if [[ ! -d "$INDIR" ]];then
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : the given input directory does not exist or is not readable. This is a FATAL error."
	echo "You have tried to use"
	echo "$INDIR"
	echo "as the path to the input run folder, but either I can't find it, or I don't have permission"
	echo "to access it. Use -h or --help for more information about how to run this script."
	echo "(Eg, $0 --help)."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi

logstamp=$(date "+%Y%m%d%H%M")
UPDAY=$(date "+%B %d, %Y at %T")

# Get the run ID, which is the name of input folder.
#
RUNID=$(basename "$INDIR")


# Make the parent output directories, if necessary.
#
if [[ ! -d "$OUTDIR/" ]];then
	mkdir "$OUTDIR/"
	echo "This folder was created by sewer-make $logstamp." > "$OUTDIR/README.txt"
fi

if [[ ! -d "$OUTDIR/updates" ]];then
	mkdir "$OUTDIR/updates"
	echo "This folder was created by sewer-make $logstamp." > "$OUTDIR/updates/README.txt"
fi

if [[ ! -d "$OUTDIR/logs" ]];then
	mkdir "$OUTDIR/logs"
	echo "This folder was created by sewer-make $logstamp." > "$OUTDIR/logs/README.txt"
fi


# Look for RUNID in any of the existing batch files. If it exists, do not proceed. Prompt the user to 
# make sure the current RUNID is unique.
#
# grep exits with a 0 if a match IS found, 1 if a match is NOT found, and 2 on error.
if [[ -d "$OUTDIR/latest" ]];then
	FOUND=$(grep -wi "$RUNID" "$OUTDIR/latest/watchdb.qbatch.txt")
	STATUS=$?
	if [[ $STATUS == 0 ]];then
		echo ""
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo "ERROR : a workflow labeled $RUNID has already been processed. This is a FATAL error."
		echo "Every run of the laboratory workflow is assigned a unique ID derived from the name of"
		echo "the folder that contains the raw workflow files. You have tried to make an update using"
		echo "$RUNID"
		echo "as the workflow ID, but that ID already exists in the current database here"
		echo "$OUTDIR/latest"
		echo "Please change the workflow run folder name to something unique. If you have multiple"
		echo "workflows from the same calendar day, add a suffix to the folder name (e.g., 2026-06-22b)"
		echo "and then run sewer-make again."
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo ""
		exit 1
	fi
fi

# Make a directory to hold the update that will be generated.
#
UPDIR="$OUTDIR/updates/$RUNID"
if [[ ! -d "$UPDIR" ]];then
	mkdir "$UPDIR"
	mkdir "$UPDIR/run_files"
else
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : an update folder labeled $RUNID already exists in $OUTDIR/updates. This is a FATAL error."
	echo "Every run of the laboratory workflow is assigned a unique ID derived from the name of"
	echo "the folder that contains the raw workflow files. You have tried to make an update using"
	echo "$RUNID"
	echo "as the workflow ID, but a folder with that name already exists in the updates directory here"
	echo "$OUTDIR/updates"
	echo "Oddly, an udate with this run ID does not appear in the existing database, so this is likely"
	echo "a result of an incomplete sewer run, where an update was made but not merged. Please rectify"
	echo "this issue and then run sewer-make again."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi


# Write all output to log file.
#
if [[ ! -f "$logf" ]];then
	logf="$OUTDIR/logs/sewer-make_$RUNID.$logstamp.log"
fi

echo ""

touch "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER make update." | tee -a "$logf"
echo "This script makes an update from the laboratory workflow files in $INDIR." | tee -a "$logf"
echo "See $logf for general warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Options for this run of SEWER make are as follows:" | tee -a "$logf"
echo "   Workflow ID:     $RUNID" | tee -a "$logf"
echo "   Workflow folder: $INDIR" | tee -a "$logf"
echo "   Update folder:   $UPDIR" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"
echo "" | tee -a "$logf"


echo "*****" | tee -a "$logf"
echo "1. Preparing the input files." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Copying all of the relevant input files for workflow $RUNID to the run_files/ folder in $UPDIR." | tee -a "$logf"
echo "The original files in $INDIR will remain untouched." | tee -a "$logf"
for fin in "$INDIR"/*
do
	f=$(basename "$fin")
	if [[ "$f" == *".ddpcr" || "$f" == *".png" ]];then
		echo "Skipping copy of $fin"
	else
		fup="$UPDIR/run_files/$f"
		cp "$fin" "$fup"
		echo "Copied $fin to $fup"
	fi
done

echo "Also copying the collection metadata file to the run_files/ directory in $UPDIR." | tee -a "$logf"
echo "This file is stored in resources/COLLECTIONS and is updated weekly with new collection metadata." | tee -a "$logf"
for f in "resources/COLLECTIONS"/*.xlsx
do
	cp "$f" "$UPDIR/run_files/collections.xlsx"
	echo "Copied $f to $UPDIR/run_files/collections.xlsx"
done

echo "Done with run file copy." | tee -a "$logf"
echo "" | tee -a "$logf"


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
echo "" | tee -a "$logf"

echo "*****" | tee -a "$logf"
echo "2. Building the update files." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "#" | tee -a "$logf"
echo "Running 2A_compileCBATCHES.R." | tee -a "$logf"
c_count=$(Rscript R/2A_compileCBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2A_compileCBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
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


echo "#" | tee -a "$logf"
echo "Running 2B_compileEBATCHES.R." | tee -a "$logf"
e_count=$(Rscript R/2B_compileEBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2B_compileEBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
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


echo "#" | tee -a "$logf"
echo "Running 2C_compileQBATCHES.R." | tee -a "$logf"
q_count=$(Rscript R/2C_compileQBATCHES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2C_compileQBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
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


echo "#" | tee -a "$logf"
echo "Running 2D_compileQRESULTS.R." | tee -a "$logf"
r_count=$(Rscript R/2D_compileQRESULTS.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2D_compileQRESULTS.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
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

echo "Finishing the merge of assay batch and result data." | tee -a "$logf"
mv "$UPDIR/update.qassay_merged.txt" "$UPDIR/update.qassay.txt"
mv "$UPDIR/update.qcontrol_merged.txt" "$UPDIR/update.qcontrol.txt"
echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"


echo "#" | tee -a "$logf"
echo "Running 2E_addSAMPLES.R." | tee -a "$logf"
s_count=$(Rscript R/2E_addSAMPLES.R "$UPDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "2E_addSAMPLES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
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
echo "   $r_count new assays (including controls)." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "##########" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "******" | tee -a "$logf"
echo "" | tee -a "$logf"

if [[ $QC == 1 ]];then
	echo "*****" | tee -a "$logf"
	echo "3. Generating QC outcomes for update $RUNID." | tee -a "$logf"
	echo "" | tee -a "$logf"
	echo "These files can be found in $UPDIR/qc." | tee -a "$logf"
	mkdir "$UPDIR/qc"

	echo "#" | tee -a "$logf"
	echo "Running 2F_qcUPDATE.R." | tee -a "$logf"
	Rscript R/2F_qcUPDATE.R "$UPDIR"
	status="${PIPESTATUS[0]}"
	if [[ "$status" != "0" ]]
	then
		echo "!!!!!!!!" | tee -a "$logf"
		echo "2F_qcUPDATE.R exited with error code $status and caused sewer-make to abort." | tee -a "$logf"
		echo "Arguments: $UPDIR" | tee -a "$logf"
		echo "sewer-make aborted during QC. This is likely a problem with the QC subroutine, not " | tee -a "$logf"
		echo "the update itself, but no promises. We recommend that you inspect the update to ensure "| tee -a "$logf"
		echo "nothing is amiss."| tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
	fi
	echo "Done." | tee -a "$logf"
fi

echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER make run $logstamp for $RUNID, started on $UPDAY, has left the chat." | tee -a "$logf"
echo "" | tee -a "$logf"

