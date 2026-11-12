#! /bin/bash

# Default option settings.
#
OUTDIR="../data"
PRECHECK=0
ALL=0
QC=1
MAKE_UPDATE=0
MERGE_UPDATE=0
REPLACE_DB=0
#

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
    -p|--precheck)
      PRECHECK=1
      shift 1 # Shift past the flag
      ;;
    -a|--all)
      ALL=1
      shift 1 # Shift past the flag
      ;;
    -q|--noqc)
      QC=0
      shift 1 # Shift past the flag
      ;;
    -u|--update)
      MAKE_UPDATE=1
      shift 1 # Shift past the flag
      ;;
    -m|--merge)
      MERGE_UPDATE=1
      shift 1 # Shift past the flag
      ;;
    -r|--replace)
      REPLACE_DB=1
      shift 1 # Shift past the flag
      ;;
    -h|--help)
      cat << EOF
###################################################################################################
SEWER: Surveillance & Exploration of Wastewater to inform Epidemiological Response.
(Timothy Driscoll; github: wvuvectors)

This is the master control program for SEWER. SEWER takes a set of laboratory workflow 
files, converts them into a set of update tables, merges the update with the existing 
tables, and then replaces the primary database files with the merge.

A minimal run of this program looks like this:
  $0 -i <MY_INPUT_DIR>
It is equivalent to this:
  $0 -i <MY_INPUT_DIR> -o ../data/ --all

The above will process the given workflow run in MY_INPUT_DIR all the way through to replacing 
the current database tables with the updated version. All QC output will be generated EXCEPT 
the precheck of the original run files. Output files will be written to the default location 
of ../data/.

If you want to run each major step separately, with full QC, this will do the trick:
  $0 -i <MY_INPUT_DIR> -o ../data/ --precheck
  $0 -i <MY_INPUT_DIR> -o ../data/ --update
  $0 -i <MY_INPUT_DIR> -o ../data/ --merge
  $0 -i <MY_INPUT_DIR> -o ../data/ --replace

Use the settings below to tweak the functionality of this script as desired.

	-i|--input <path_to_workflow_run_folder> 
	Folder containing the original run files.
	REQUIRED.
	
	-o|--output <path_to_output_folder>
	DEFAULT: ../data
	
	-p|--precheck
	Run QC on the given original workflow run and stop. Supercedes all settings, including --all.
	
	-a|--all
	Run all of qc, update, merge, and replace. Supercedes all settings for these individual 
	options, if present. Functionally equivalent to 
	$0 --update --merge --replace
	
	-q|--noqc
	Do NOT include the QC routines for any of the selected options. QC is turned on by default. 
	This setting supercedes the --all flag.
	
	-u|--update
	Generate the update files and stop. Superceded by the --all flag, if set.
	
	-m|--merge
	Merge an existing update with the current database and stop. Superceded by the --all flag, if set.
	Requires an update folder in OUTDIR.
	
	-r|--replace
	Replace the current database with an existing merge and stop. Superceded by the --all flag, if set.
	Requires a merge folder in OUTDIR.
	
	-h|--help
	Print this helpful block of text.

###################################################################################################

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

# Get the run ID, which is the name of input folder. This is also used as the folder name 
# for the update and incremental merge steps.
#
RUNID=$(basename "$INDIR")


# Not all combination of settings are compatible, as some override others. Make sure that 
# we have a working combination of settings. Before moving forward.
#
if [[ $PRECHECK == 1 ]];then
	# If a precheck has been requested, set all other options to FALSE.
	ALL=0
	QC=0
	MAKE_UPDATE=0
	MERGE_UPDATE=0
	REPLACE_DB=0
elif [[ $MAKE_UPDATE == 0 && $MERGE_UPDATE == 0 && $REPLACE_DB == 0 ]];then
	ALL=1
fi

if [[ $ALL == 0 && $MAKE_UPDATE == 1 && $MERGE_UPDATE == 0 && $REPLACE_DB == 1 ]];then
	# This combination is just not allowed.
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : You are trying to skip the incremental merge step after making a new update. This is a FATAL error."
	echo "By skipping the -m|--merge flag, you have tried to skip a required step in this pipeline."
	echo "Use -h or --help for more information about how to run this script (Eg, $0 --help)."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi


if [[ $ALL == 0 && $MAKE_UPDATE == 0 && $MERGE_UPDATE == 1 ]];then
	# For this combination, an update must exist in $OUTDIR.
	if [[ ! -d "$OUTDIR/updates/" || ! -d "$OUTDIR/updates/$RUNID" ]];then
		echo ""
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo "ERROR : You are trying to merge an existing update, but no update exists. This is a FATAL error."
		echo "By setting the -m|--merge flag but NOT the -u|--update flag, You have tried to use"
		echo "an existing update $RUNID located here"
		echo "$OUTDIR/updates/$RUNID"
		echo "Either I can't find this update, or I don't have permission to access it. Is it possible"
		echo "that you still need to make this update? If so, re-run with both the --update and --merge"
		echo "flags set. If you already made the update, it's possible the folder has been moved; "
		echo "Verify that the above path is correct and the folder exists and I can access it."
		echo "Use -h or --help for more information about how to run this script (Eg, $0 --help)."
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo ""
		exit 1
	fi
fi

if [[ $ALL == 0 && $MAKE_UPDATE == 0 && $MERGE_UPDATE == 0 && $REPLACE_DB == 1 ]];then
	# For this combination, a merge set must exist in $OUTDIR.
	if [[ ! -d "$OUTDIR/merges/" || ! -d "$OUTDIR/merges/$RUNID" ]];then
		echo ""
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo "ERROR : You are trying to replace the existing tables with an incremental merge,"
		echo "but no such merge exists. This is a FATAL error."
		echo "By setting the -r|--replace flag by itself, You have tried to use an existing merge of $RUNID located here"
		echo "$OUTDIR/merges/$RUNID"
		echo "Either I can't find this merge, or I don't have permission to access it. Is it possible"
		echo "that you still need to make this merge? If so, re-run with both the --merge and --replace"
		echo "flags set. If you already made the merge, it's possible the folder has been moved; "
		echo "Verify that the above path is correct and the folder exists and I can access it."
		echo "Use -h or --help for more information about how to run this script (Eg, $0 --help)."
		echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		echo ""
		exit 1
	fi
fi


# Get the current date and time for logging purposes.
#
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")
TS=$(date +%s)
logstamp=$(date "+%Y%m%d%H%M")


# Make sure the main output folder exist, or create it.
# This folder will hold the updates, merges, and the latest database tables.
#
if [[ ! -d "$OUTDIR" ]];then
	mkdir "$OUTDIR"
	echo "This folder was created by sewer.sh $logstamp." > "$OUTDIR/README.txt"
fi

if [[ ! -d "$OUTDIR/logs" ]];then
	mkdir "$OUTDIR/logs"
	echo "This folder was created by sewer.sh $logstamp." > "$OUTDIR/README.txt"
fi


# Write all output for this script to a common log file.
#
logf="$OUTDIR/logs/sewer_$logstamp.$RUNID.log"
if [[ -f "$logf" ]]
then
	rm "$logf"
fi

echo ""

touch "$logf"
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER: Surveillance & Exploration of Wastewater to inform Epidemiological Response." | tee -a "$logf"
echo "Initiated $UPDAY" | tee -a "$logf"
echo "See $logf for general warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Options for this run of SEWER are as follows:" | tee -a "$logf"
echo "   Workflow ID:     $RUNID" | tee -a "$logf"
echo "   Workflow folder: $INDIR" | tee -a "$logf"
echo "   Output folder:   $OUTDIR" | tee -a "$logf"
echo "" | tee -a "$logf"
if [[ $PRECHECK == 1 ]];then
	echo "   QC the workflow:  YES" | tee -a "$logf"
	echo "   QC steps:          NO" | tee -a "$logf"
	echo "   Make update:       NO" | tee -a "$logf"
	echo "   Merge update:      NO" | tee -a "$logf"
	echo "   Replace database:  NO" | tee -a "$logf"
elif [[ $ALL == 1 ]];then
	echo "   QC the workflow:   NO" | tee -a "$logf"
	echo "   QC steps:         YES" | tee -a "$logf"
	echo "   Make update:      YES" | tee -a "$logf"
	echo "   Merge update:     YES" | tee -a "$logf"
	echo "   Replace database: YES" | tee -a "$logf"
else
	echo "   QC the workflow:   NO" | tee -a "$logf"
fi
if [[ $PRECHECK == 0 && $ALL == 0 ]];then
	if [[ $QC == 1 ]];then
		echo "   QC the steps:     YES" | tee -a "$logf"
	else
		echo "   QC the steps:      NO" | tee -a "$logf"
	fi
	if [[ $MAKE_UPDATE == 1 ]];then
		echo "   Make update:      YES" | tee -a "$logf"
	else
		echo "   Make update:       NO" | tee -a "$logf"
	fi
	if [[ $MERGE_UPDATE == 1 ]];then
		echo "   Merge update:     YES" | tee -a "$logf"
	else
		echo "   Merge update:      NO" | tee -a "$logf"
	fi
	if [[ $REPLACE_DB == 1 ]];then
		echo "   Replace database: YES" | tee -a "$logf"
	else
		echo "   Replace database:  NO" | tee -a "$logf"
	fi
fi
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"

echo "" | tee -a "$logf"

if [[ $PRECHECK == 1 ]];then
	echo "Running a pre-check QC on the original (input) workflow files." | tee -a "$logf"
	echo "The output from this QC will appear in the input dir ($INDIR)." | tee -a "$logf"
	echo "Note: No other tasks can be run in conjunction with the pre-check." | tee -a "$logf"
	./sewer_check.sh --input "$INDIR" --logfile "$logf"
	status="${PIPESTATUS[0]}"
	if [[ "$status" != "0" ]];then
		echo "pre-check FAIL!"
	else
		echo "" | tee -a "$logf"
		echo "SEWER pre-check run of $UPDAY has finished. QC results can be found with the workflow files in:" | tee -a "$logf"
		echo "$INDIR" | tee -a "$logf"
		echo "The pre-check completed without error, though that does not necessarily mean a successful QC." | tee -a "$logf"
		echo "Inspect the QC result files to determine the next steps." | tee -a "$logf"
		echo "" | tee -a "$logf"
	fi
fi


if [[ $ALL == 1 || $MAKE_UPDATE == 1 ]];then
	# Make the update.
	echo "Making an update for workflow $RUNID." | tee -a "$logf"
	echo "The output from this make will appear in $OUTDIR/updates/$RUNID." | tee -a "$logf"
	if [[ $QC == 0 ]];then
		./sewer_make.sh --input "$INDIR" --output "$OUTDIR" --logfile "$logf" --noqc
	else
		./sewer_make.sh --input "$INDIR" --output "$OUTDIR" --logfile "$logf"
	fi
	status="${PIPESTATUS[0]}"
	if [[ "$status" != "0" ]];then
		echo "make FAIL!"
	else
		echo "" | tee -a "$logf"
		echo "SEWER make run of $UPDAY has finished. Update files can be found in:" | tee -a "$logf"
		echo "$OUTDIR/updates/$RUNID" | tee -a "$logf"
		echo "" | tee -a "$logf"
	fi
fi

if [[ $ALL == 1 || $MERGE_UPDATE == 1 ]];then
	# Merge the update.
	echo "Merging the update for workflow $RUNID with the primary database." | tee -a "$logf"
	echo "The output from this merge will appear in $OUTDIR/merges/$RUNID." | tee -a "$logf"
	if [[ $QC == 0 ]];then
		./sewer_merge.sh --input "$OUTDIR/updates/$RUNID" --output "$OUTDIR" --logfile "$logf" --noqc
	else
		./sewer_merge.sh --input "$OUTDIR/updates/$RUNID" --output "$OUTDIR" --logfile "$logf"
	fi
fi

if [[ $ALL == 1 || $REPLACE_DB == 1 ]];then
	# Replace the existing database.
	echo "Replacing the existing primary database with the merge of $RUNID." | tee -a "$logf"
	echo "The output from this merge will appear in $OUTDIR/latest." | tee -a "$logf"
	echo "The previous version of the primary database, before $RUNID, will appear in $OUTDIR/latest-bk." | tee -a "$logf"
	if [[ $QC == 0 ]];then
		./sewer_replace.sh --input "$OUTDIR/merges/$RUNID" --output "$OUTDIR" --logfile "$logf" --noqc
	else
		./sewer_replace.sh --input "$OUTDIR/merges/$RUNID" --output "$OUTDIR" --logfile "$logf"
	fi
fi

echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER run $logstamp for $RUNID, started on $UPDAY, has left the chat." | tee -a "$logf"
echo "" | tee -a "$logf"
