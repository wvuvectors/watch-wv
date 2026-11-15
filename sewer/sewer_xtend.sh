#! /bin/bash

# Default option settings.
#
DATADIR="../data/latest"
OUTDIR="../data/xtend"
ALL=0
QC=1
MAKE_NWSS=0
MAKE_WVDASH=0
MAKE_WVBEASY=0
MU_NWSS_F="/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/mu_nwss.LATEST.csv"
#

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      DATADIR="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -o|--output)
      OUTDIR="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -m|--marshall)
      MU_NWSS_F="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -a|--all)
      ALL=1
      shift 1 # Shift past the flag
      ;;
    -q|--noqc)
      QC=0
      shift 1 # Shift past the flag
      ;;
    -n|--nwss)
      MAKE_NWSS=1
      shift 1 # Shift past the flag
      ;;
    -w|--watch)
      MAKE_WVDASH=1
      shift 1 # Shift past the flag
      ;;
    -b|--breatheeasy)
      MAKE_WVBEASY=1
      shift 1 # Shift past the flag
      ;;
    -h|--help)
      cat << EOF
###################################################################################################
SEWER: Surveillance & Exploration of Wastewater to inform Epidemiological Response.
(Timothy Driscoll; github: wvuvectors)

This is the XTEND component of the SEWER program. SEWER takes a set of laboratory workflow 
files, converts them into a set of update tables, merges the update with the existing 
tables, and then replaces the primary database files with the merge. The XTEND component takes 
the SEWER database and creates files for external references including NWSS, WV BreatheEasy, 
and the WV Dashboard.

A minimal run of this program looks like this:
  $0 -i <SEWER_DATA_DIR>
It is equivalent to this:
  $0 -i <SEWER_DATA_DIR> --all

The above will use the data in SEWER_DATA_DIR to generate files that can be used as input 
into the NWSS site, WV BreathEasy site, and the WV state WaTCH dashboard. The output files 
are written to the extref/ directory in SEWER_DATA_DIR.

If you want to run each step separately, with full QC, this will do the trick:
  $0 -i <SEWER_DATA_DIR> --nwss
  $0 -i <SEWER_DATA_DIR> --watch
  $0 -i <SEWER_DATA_DIR> --breatheeasy

Use the settings below to tweak the functionality of this script as desired.

	-i|--input <path_to_watch_data_folder> 
	Folder containing the watch data tables to use.
	DEFAULT: ../data/latest
	
	-o|--output <path_to_output_folder> 
	Path to the parent destination folder for the chosen xtend processes.
	DEFAULT: ../data/xtend

	-m|--marshall <path_to_mu_data_file> 
	File with data from MUIDSL to include with the WaTCH data. Required NWSS csv format.
	DEFAULT: /Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/mu_nwss.LATEST.csv
	
	-a|--all
	Run all of qc, nwss, watch, and breatheeasy. Supercedes all settings for these individual 
	options, if present. Functionally equivalent to 
	$0 --nwss --watch --breatheeasy
	
	-q|--noqc
	Do NOT include the QC routines for any of the selected options. QC is turned on by default. 
	This setting supercedes the --all flag.
	
	-u|--nwss
	Generate the file for transfer to the NWSS site and stop. Superceded by the --all flag, if set.
	
	-w|--watch
	Generate the file for use by the WaTCH dashboard and stop. Superceded by the --all flag, if set.
	
	-b|--breatheeasy
	Generate the file for use by the WV breatheeasy dashboard and stop. Superceded by the --all flag, if set.
	
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

if [[ ! -d "$DATADIR" ]];then
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "ERROR : the given input directory does not exist or is not readable. This is a FATAL error."
	echo "You have tried to use"
	echo "$DATADIR"
	echo "as the path to the WaTCH data folder, but either I can't find it, or I don't have permission"
	echo "to access it. Use -h or --help for more information about how to run this script."
	echo "(Eg, $0 --help)."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit 1
fi

# Get the current date and time for logging purposes.
#
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")
TS=$(date +%s)
logstamp=$(date "+%Y%m%d%H%M")


# Make sure the main output folder (extref) exist, or create it.
# This folder will hold the updates, merges, and the latest database tables.
#
if [[ ! -d "$OUTDIR" ]];then
	mkdir "$OUTDIR"
	echo "This folder was created by sewer_xtend.sh $logstamp." > "$OUTDIR/README.txt"
fi
if [[ ! -d "$OUTDIR/logs" ]];then
	mkdir "$OUTDIR/logs"
	echo "This folder was created by sewer.sh $logstamp." > "$OUTDIR/logs/README.txt"
fi


# Write all output for this script to a common log file.
#
logf="$OUTDIR/logs/sewer_xtend_$logstamp.log"
if [[ -f "$logf" ]]
then
	rm "$logf"
fi

echo ""

touch "$logf"
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER-XTEND: Surveillance & Exploration of Wastewater to inform Epidemiological Response." | tee -a "$logf"
echo "Initiated $UPDAY" | tee -a "$logf"
echo "See $logf for general warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Options for this run of SEWER-XTEND are as follows:" | tee -a "$logf"
echo "   Data folder:   $DATADIR" | tee -a "$logf"
echo "   Output folder: $OUTDIR" | tee -a "$logf"
echo "----------------------------------------------" | tee -a "$logf"
if [[ $ALL == 1 ]];then
	echo "   QC steps:                 YES" | tee -a "$logf"
	echo "   Make NWSS file:           YES" | tee -a "$logf"
	echo "   Make WaTCH dash file:     YES" | tee -a "$logf"
	echo "   Make WV BreatheEasy file: YES" | tee -a "$logf"
else
	if [[ $QC == 1 ]];then
		echo "   QC steps:                 YES" | tee -a "$logf"
	else
		echo "   QC steps:                  NO" | tee -a "$logf"
	fi
	if [[ $MAKE_NWSS == 1 ]];then
		echo "   Make NWSS file:           YES" | tee -a "$logf"
	else
		echo "   Make NWSS file:            NO" | tee -a "$logf"
	fi
	if [[ $MAKE_WVDASH == 1 ]];then
		echo "   Make WaTCH dash file:     YES" | tee -a "$logf"
	else
		echo "   Make WaTCH dash file:      NO" | tee -a "$logf"
	fi
	if [[ $MAKE_WVBEASY == 1 ]];then
		echo "   Make WV BreatheEasy file: YES" | tee -a "$logf"
	else
		echo "   Make WV BreatheEasy file:  NO" | tee -a "$logf"
	fi
fi
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"

echo "" | tee -a "$logf"

if [[ $ALL == 1 || $MAKE_NWSS == 1 ]];then

	echo "Making an update file for the NWSS site." | tee -a "$logf"
	echo "The output from this make will appear in $OUTDIR/nwss." | tee -a "$logf"
	echo "An archive will also be created in $OUTDIR/nwss/archives." | tee -a "$logf"
	
	# Make sure the proper NWSS output directories exist.
	#
	OUTDIR_NWSS="$OUTDIR/nwss"
	if [[ ! -d "$OUTDIR_NWSS" ]];then
		mkdir "$OUTDIR_NWSS"
		echo "This folder was created by sewer_xtend.sh $logstamp." > "$OUTDIR_NWSS/README.txt"
	fi
	if [[ ! -d "$OUTDIR_NWSS/archive" ]];then
		mkdir "$OUTDIR_NWSS/archive"
		echo "This folder was created by sewer_xtend.sh $logstamp." > "$OUTDIR_NWSS/archive/README.txt"
	fi
	
	echo "#" | tee -a "$logf"
	echo "Running 9A_xtendNWSS.R." | tee -a "$logf"
	nwss_count=$(Rscript R/9A_xtendNWSS.R "$DATADIR" "$OUTDIR_NWSS" "$MU_NWSS_F" "logstamp")
	status="${PIPESTATUS[0]}"
	if [[ "$status" != "0" ]]
	then
		echo "!!!!!!!!" | tee -a "$logf"
		echo "9A_xtendNWSS.R exited with error code $status and caused sewer-xtend to abort." | tee -a "$logf"
		echo "Options passed to this script:" | tee -a "$logf"
		echo "$DATADIR" | tee -a "$logf"
		echo "$OUTDIR_NWSS" | tee -a "$logf"
		echo "$MU_NWSS_F" | tee -a "$logf"
		echo "$logstamp" | tee -a "$logf"
		echo "sewer-xtend aborted during phase A (NWSS file make)." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
	fi
	if [[ $nwss_count < 1 ]]
	then
		echo "!!!!!!!!" | tee -a "$logf"
		echo "FATAL: There is no NWSS data in the output." | tee -a "$logf"
		echo "FATAL: This has caused sewer to abort." | tee -a "$logf"
		echo "!!!!!!!!" | tee -a "$logf"
		exit 1
	fi
	echo "Done." | tee -a "$logf"
	echo "" | tee -a "$logf"
	
# 	if [[ $QC == 0 ]];then
# 		./sewer_make.sh --input "$INDIR" --output "$OUTDIR" --logfile "$logf" --noqc
# 	else
# 		./sewer_make.sh --input "$INDIR" --output "$OUTDIR" --logfile "$logf"
# 	fi
# 	status="${PIPESTATUS[0]}"
# 	if [[ "$status" != "0" ]];then
# 		echo "make FAIL!"
# 	else
# 		echo "" | tee -a "$logf"
# 		echo "SEWER make run of $UPDAY has finished. Update files can be found in:" | tee -a "$logf"
# 		echo "$OUTDIR/updates/$RUNID" | tee -a "$logf"
# 		echo "" | tee -a "$logf"
# 	fi
fi

echo "" | tee -a "$logf"
echo "SEWER-XTEND has left the chat." | tee -a "$logf"
echo "Run $logstamp, started on $UPDAY." | tee -a "$logf"
echo "" | tee -a "$logf"

