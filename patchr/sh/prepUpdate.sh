#! /bin/bash


updir="$1"
WD=$(pwd)

status=0

WATCHFILE_MU='/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/mu_dashboard.LATEST.tsv'
WATCHFILE_MUBK='/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/ARCHIVED/mu_dashboard.$now.tsv'

if [ -f "$WATCHFILE_MU" ]
then
	echo "Configuring MU dashboard file $WATCHFILE_MU."
	cp "$WATCHFILE_MU" "$WATCHFILE_MUBK"
	perl -pi -e 's/\r$//' "$WATCHFILE_MU"
	sed -i '' -e '$a\' "$WATCHFILE_MU"
	echo "Done. The original, unprocessed MU file has been archived as $WATCHFILE_MUBK."
fi


echo "Switching to update directory $updir."
cd "$updir"

echo "Creating directory batches/ to hold the validated batch files."
if [ -d "batches/" ]
then
	echo "Directory already exists so I'll just use it."
else
	mkdir "batches/"
fi

echo "Preparing unprocessed batch files listed in update.batch_files.txt."
echo "Validated batch files will be stored in batches/ and update.batch_files.txt will be revised."
echo "Further steps in patchr will use the files in batches/."

mv update.batch_files.txt update.batch_files.original.txt

while read -r f
do

	# Copy the batch files to the batches directory for validation.
	rfile=$(basename "$f")
	echo "Copying file $rfile to batches/$rfile."
	cp "$f" "batches/$rfile"
	
	echo "$updir/batches/$rfile" >> update.batch_files.txt
	
	# Fix known file issues, like Windows line endings.
	echo "Repairing known issues in $rfile."
	if [[ "$rfile" == *.csv ]]
	then
		perl -pi -e 's/\r$//' "batches/$rfile"
		sed $'1s/\xef\xbb\xbf//' < "batches/$rfile" > "batches/$rfile.tmp"
		mv "batches/$rfile.tmp" "batches/$rfile"
		echo "Done."

		is_ddpcr=$(head -n 1 "batches/$rfile" | grep -c "Well,")
		is_asset=$(grep -c "Asset Tag ID," "batches/$rfile")

		if [[ "$is_ddpcr" == "1" ]]
		then
			echo "File $rfile appears to be a ddPCR results file."
			echo "  Fixing known incompatibilities in $rfile."
			perl -pi -e 's/,Taget,/,Target,/i' "batches/$rfile"
			perl -pi -e 's/,N1,/,SARS-CoV-2 N1,/i' "batches/$rfile"
			perl -pi -e 's/,N2,/,SARS-CoV-2 N2,/i' "batches/$rfile"
			perl -pi -e 's/,RP,/,Human RNase P,/i' "batches/$rfile"
			perl -pi -e 's/µ/u/i' "batches/$rfile"
			echo "Done."
		elif [[ "$is_asset" == "1" ]]
		then
			echo "File $rfile appears to be an asset file."
			echo "It has been copied to batches/."
			#echo "You may want to process this file separately, using perl/2_updateSamples.pl."
		else
			echo "File $rfile is a csv file of unknown data format."
			echo "It has been copied to batches/ but will be ignored."
		fi
	fi
done <update.batch_files.original.txt

#mv update.batch_files.txt update.batch_files.original.txt
#mv update.batch_files.txt.tmp update.batch_files.txt


cd "$WD"

exit $status

