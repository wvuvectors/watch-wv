#! /bin/bash

while getopts ":hi:" opt; do
	case $opt in
		h)
			echo "help not available."
			exit 1
			;;
		i)
			INDIR=$OPTARG
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


# Get the current date and time for logging purposes.
#
START=$(date "+%F_%H-%M")
UPDAY=$(date "+%B %d, %Y at %T")
TS=$(date +%s)

# Get the run ID, which is the name of input folder.
#
RUNID=$(basename "$INDIR")


# Write all output to log file.
#
logf="$INDIR/sewerQC.log"
if [ -f "$logf" ]
then
	rm "$logf"
fi

echo ""

touch "$logf"
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "SEWER: Surveillance & Exploration of Wastewater to inform Epidemiological Response." | tee -a "$logf"
echo "Initial data quality check workflow, initiated from sewerQC.sh." | tee -a "$logf"
echo "Script run datetime stamp: $START" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "See $logf for warnings, errors, and other important information." | tee -a "$logf"
echo "Run ID: $RUNID" | tee -a "$logf"
echo "Input directory: $INDIR" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"


echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"
echo "Running 0_qcBATCHES.R." | tee -a "$logf"
batch_qc_res=$(Rscript R/0_qcBATCHES.R "$INDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcBATCHES.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "This is due to a syntax error, and not necessarily a problem with the data. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

if [[ "$batch_qc_res" == "1" ]];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcBATCHES.R returned error code $batch_qc_res:" | tee -a "$logf"
	echo "This run is missing at least one required batch file. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
elif [[ "$batch_qc_res" == "2" ]];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcBATCHES.R returned error code $batch_qc_res:" | tee -a "$logf"
	echo "One or more of the batch files is missing a valid batch ID. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
elif [[ "$batch_qc_res" == "3" ]];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcBATCHES.R returned error code $batch_qc_res:" | tee -a "$logf"
	echo "This run is missing at least one required batch file, AND one or more of the existing batch files is missing a valid batch ID."| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi
echo "Done." | tee -a "$logf"
echo "#############################################" | tee -a "$logf"



echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "#############################################" | tee -a "$logf"
echo "Running 0_qcDATA.R." | tee -a "$logf"
data_qc_res=$(Rscript R/0_qcDATA.R "$INDIR")
status="${PIPESTATUS[0]}"
if [[ "$status" != "0" ]]
then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcDATA.R exited with error code $status and caused sewer to abort." | tee -a "$logf"
	echo "This is due to a syntax error, and not necessarily a problem with the data. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
	exit 1
fi

if [[ "$data_qc_res" == "1" ]];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcDATA.R returned error code $data_qc_res:" | tee -a "$logf"
	echo "Too many wells failed to meet minimum accepted partitions. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
elif [[ "$data_qc_res" == "2" ]];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcDATA.R returned error code $data_qc_res:" | tee -a "$logf"
	echo "One or more of the data files is missing a valid batch ID. "| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
elif [[ "$data_qc_res" == "3" ]];then
	echo "!!!!!!!!" | tee -a "$logf"
	echo "0_qcDATA.R returned error code $data_qc_res:" | tee -a "$logf"
	echo "Too many wells failed to meet minimum accepted partitions, AND one or more of the data files is missing a valid batch ID."| tee -a "$logf"
	echo "Check any output files that were generated, and the log file, for hints as to the problem."| tee -a "$logf"
	echo "!!!!!!!!" | tee -a "$logf"
fi
echo "Done." | tee -a "$logf"
echo "#############################################" | tee -a "$logf"



echo "" | tee -a "$logf"
echo "" | tee -a "$logf"
echo "sewerQC run of $UPDAY has left the chat." | tee -a "$logf"
echo "" | tee -a "$logf"

