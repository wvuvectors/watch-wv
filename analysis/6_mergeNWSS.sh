#! /bin/bash

indir=$1

cd "$indir"

# get the current datetime
now=$(date +'%Y-%m-%d.%k_%M')


mu_len=$(wc -l < READY/mu_nwss.LATEST.csv)
mu_len=$((mu_len-1))

wvu_len=$(wc -l < READY/wvu_nwss.LATEST.csv)
wvu_len=$((wvu_len-1))

if [ $mu_len -lt 1 ]
then
	echo "The MU data file appears to be empty."
	echo "No merge will be performed."
	exit 1
fi

if [ $wvu_len -lt 1 ]
then
	echo "The WVU data file appears to be empty."
	echo "No merge will be performed."
	exit 1
fi


if [ -f "READY/merged_nwss.LATEST.csv" ]
then
	echo "A merged file already exists! The new merge will OVERWRITE this file."
	echo "In most cases, this is desired behavior as each new NWSS update contains all data."
	echo "Just in case the existing merge file is required for some reason, I am saving it as 'merged_nwss.OLD.txt' for now."
	echo "I recommend deleting or archiving 'merged_nwss.OLD.csv' properly ASAP to avoid confusion."
	mv READY/merged_nwss.LATEST.csv READY/merged_nwss.OLD.csv
fi


# Convert Windows line endings to unix because of course
perl -pi -e 's/\r$//' READY/mu_nwss.LATEST.csv
perl -pi -e 's/\r$//' READY/wvu_nwss.LATEST.csv

# Remove any trailing empty columns
# This can occur when exporting from damn Excel
perl -pi -e 's/,+$//' READY/mu_nwss.LATEST.csv
perl -pi -e 's/,+$//' READY/wvu_nwss.LATEST.csv

# Remove any empty rows
# This can occur when exporting from damn Excel
sed -i '' '/^[[:space:]]*$/d' READY/mu_nwss.LATEST.csv
sed -i '' '/^[[:space:]]*$/d' READY/wvu_nwss.LATEST.csv

# Ensure both files contain a line ending at the end of the file
tail -c1 < "READY/mu_nwss.LATEST.csv" | read -r _ || echo >> "READY/mu_nwss.LATEST.csv"
tail -c1 < "READY/wvu_nwss.LATEST.csv" | read -r _ || echo >> "READY/wvu_nwss.LATEST.csv"

#sed -i '' -e '$a\' READY/mu_nwss.LATEST.csv
#sed -i '' -e '$a\' READY/wvu_nwss.LATEST.csv


# remove the header line from the WVU data file
# since this is the second file in the concat we don't need the header
tail -n +2 READY/wvu_nwss.LATEST.csv > READY/wvu_tmp.csv


# concatenate MU and WVU data files
cat READY/mu_nwss.LATEST.csv READY/wvu_tmp.csv > READY/merged_nwss.LATEST.csv

# Recalculate the length of each individul file since processing might have changed this
mu_len=$(wc -l < READY/mu_nwss.LATEST.csv)
mu_len=$((mu_len-1))
wvu_len=$(wc -l < READY/wvu_nwss.LATEST.csv)
wvu_len=$((wvu_len-1))
sum_len=$((mu_len + wvu_len))

# Calculate the length of the merged file
merged_len=$(wc -l < READY/merged_nwss.LATEST.csv)
merged_len=$((merged_len-1))



if [ $merged_len -lt 1 ]
then
	echo "******"
	echo "The merged data file appears to be empty."
	echo "The merge has occurred but all intermediate files will remain so you can find the problem."
	echo "No files have been archived!"
	echo "******"
	exit 1
fi

echo "Data rows in original MU file : $mu_len"
echo "Data rows in original WVU file: $wvu_len"
echo "Data rows in merged file      : $merged_len"

if [ $sum_len -ne $merged_len ]
then
	echo "******"
	echo "The merged file ($merged_len) does not contain the expected number of data rows ($sum_len)."
	echo "The merge has occurred but all intermediate files will remain so you can find the problem."
	echo "No files have been archived!"
	echo "******"
	exit 1
else
	echo "******"
	echo "Merge complete and validated."
	echo "File 'READY/merged_nwss.LATEST.csv' is now ready for submission to NWSS."
	echo "Archiving original files."
	echo "******"
fi

# archive the individual data files
mv READY/mu_nwss.LATEST.csv UNMERGED/mu_nwss.$now.csv
mv READY/wvu_nwss.LATEST.csv UNMERGED/wvu_nwss.$now.csv

cp READY/merged_nwss.LATEST.csv MERGED/merged_nwss.$now.csv

# clean up tmp files
rm READY/wvu_tmp.csv

exit 0
