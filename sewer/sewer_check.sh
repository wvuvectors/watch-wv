#! /bin/bash

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INDIR="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -l|--logfile)
      logf="$2"
      shift 2 # Shift past the flag and its value
      ;;
    -h|--help)
      cat << EOF
SEWER precheck script.
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

# Get the run ID, which is the name of input folder.
#
RUNID=$(basename "$INDIR")

mkdir "$INDIR/qc"

# Write all output to log file in the qc folder, unless passed a logfile.
#
if [[ ! -f "$logf" ]];then
	logstamp=$(date "+%Y%m%d%H%M")
	logf="$INDIR/qc/precheck_$RUNID.$logstamp.log"
fi

echo ""

touch "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER precheck." | tee -a "$logf"
echo "Initial data quality check of laboratory workflow files." | tee -a "$logf"
echo "See $logf for general warnings, errors, and other important information." | tee -a "$logf"
echo "" | tee -a "$logf"
echo "Options for this run of SEWER precheck are as follows:" | tee -a "$logf"
echo "   Workflow ID:     $RUNID" | tee -a "$logf"
echo "   Workflow folder: $INDIR" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "########################" | tee -a "$logf"

echo "" | tee -a "$logf"

echo "Running 1A_qcBATCHES.R." | tee -a "$logf"
batch_qc_res=$(Rscript R/1A_qcBATCHES.R "$INDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "1A_qcBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "This is due to a syntax error, and not necessarily a problem with the data. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "##############################" | tee -a "$logf"
echo "Results from the quality control check of the workflow batch files:" | tee -a "$logf"
echo "" | tee -a "$logf"
if [[ "$batch_qc_res" == "0" ]];then
	echo "1A_qcBATCHES.R returned result code $batch_qc_res:" | tee -a "$logf"
	echo "All QC checks were successful, and no major problems were detected. "| tee -a "$logf"
elif [[ "$batch_qc_res" == "1" ]];then
	echo "1A_qcBATCHES.R returned error code $batch_qc_res:" | tee -a "$logf"
	echo "This run is missing at least one required batch file. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
elif [[ "$batch_qc_res" == "2" ]];then
	echo "1A_qcBATCHES.R returned error code $batch_qc_res:" | tee -a "$logf"
	echo "One or more of the batch files is missing a valid batch ID. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
elif [[ "$batch_qc_res" == "3" ]];then
	echo "1A_qcBATCHES.R returned error code $batch_qc_res:" | tee -a "$logf"
	echo "This run is missing at least one required batch file, AND one or more of the existing batch files is missing a valid batch ID."| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
fi
echo "##############################" | tee -a "$logf"

echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"

echo "Running 1B_qcDATA.R." | tee -a "$logf"
data_qc_res=$(Rscript R/1B_qcDATA.R "$INDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "1B_qcDATA.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "This is due to a syntax error, and not necessarily a problem with the data. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

echo "##############################" | tee -a "$logf"
echo "Results from the quality control check of the workflow data:" | tee -a "$logf"
echo "" | tee -a "$logf"
if [[ "$data_qc_res" == "0" ]];then
	echo "1B_qcDATA.R returned result code $data_qc_res:" | tee -a "$logf"
	echo "All QC checks were successful, and no major problems were detected. "| tee -a "$logf"
elif [[ "$data_qc_res" == "1" ]];then
	echo "1B_qcDATA.R returned result code $data_qc_res:" | tee -a "$logf"
	echo "Too many wells failed to meet minimum accepted partitions. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
elif [[ "$data_qc_res" == "2" ]];then
	echo "1B_qcDATA.R returned result code $data_qc_res:" | tee -a "$logf"
	echo "One or more of the data files is missing a valid batch ID. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
elif [[ "$data_qc_res" == "3" ]];then
	echo "1B_qcDATA.R returned result code $data_qc_res:" | tee -a "$logf"
	echo "Too many wells failed to meet minimum accepted partitions, AND one or more of the data files is missing a valid batch ID."| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
fi
echo "##############################" | tee -a "$logf"

echo "Done." | tee -a "$logf"
echo "" | tee -a "$logf"

exit 0
