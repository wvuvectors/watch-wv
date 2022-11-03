#! /bin/bash

indir=$1

cd "$indir"

if [ -f "merged_nwss.LATEST.txt" ]
then
	echo "A merged file already exists! The new merge will OVERWRITE this file."
	echo "In most cases, this is desired behavior as each new NWSS update contains all data."
	echo "Just in case the existing merge file is required for some reason, I am saving it as 'merged_nwss.OLD.txt' for now."
	echo "I recommend deleting or archiving this file properly ASAP to avoid confusion."
	mv merged_nwss.LATEST.txt merged_nwss.OLD.txt
fi

# get the current datetime
now=$(date +'%Y-%m-%d.%k_%M')

# remove the header line from the ZPM data file
# since this is the second file in the concat we don't need the header
tail -n +2 zpm_nwss.LATEST.txt > zpm_tmp.txt

# concatenate MU and WVU data files
concat mu_nwss.LATEST.txt zpm_tmp.txt > merged_nwss.LATEST.txt

# archive the individual data files
mv mu_nwss.LATEST.txt UNMERGED/mu_nwss.$now.txt
mv zpm_nwss.LATEST.txt UNMERGED/zpm_nwss.$now.txt

cp merged_nwss.LATEST.txt MERGED/merged_nwss.$now.txt

# clean up tmp files
rm zpm_tmp.txt
