#! /bin/bash

indir=$1

cd "$indir"

if [ -f "READY/merged_nwss.LATEST.txt" ]
then
	echo "A merged file already exists! The new merge will OVERWRITE this file."
	echo "In most cases, this is desired behavior as each new NWSS update contains all data."
	echo "Just in case the existing merge file is required for some reason, I am saving it as 'merged_nwss.OLD.txt' for now."
	echo "I recommend deleting or archiving 'merged_nwss.OLD.txt' properly ASAP to avoid confusion."
	mv READY/merged_nwss.LATEST.txt READY/merged_nwss.OLD.txt
fi

# get the current datetime
now=$(date +'%Y-%m-%d.%k_%M')

# remove the header line from the ZPM data file
# since this is the second file in the concat we don't need the header
tail -n +2 READY/wvu_nwss.LATEST.txt > READY/wvu_tmp.txt

# concatenate MU and WVU data files
concat READY/mu_nwss.LATEST.txt READY/wvu_tmp.txt > READY/merged_nwss.LATEST.txt

# archive the individual data files
mv READY/mu_nwss.LATEST.txt UNMERGED/mu_nwss.$now.txt
mv READY/wvu_nwss.LATEST.txt UNMERGED/wvu_nwss.$now.txt

cp READY/merged_nwss.LATEST.txt MERGED/merged_nwss.$now.txt

# clean up tmp files
rm READY/zpm_tmp.txt
