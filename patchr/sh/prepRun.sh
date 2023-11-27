#! /bin/bash


WATCHFILE_MU='/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/mu_dashboard.LATEST.tsv'
WATCHFILE_MUBK='/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/ARCHIVED/mu_dashboard.$now.tsv'

if [ -f "$WATCHFILE_MU" ]
then
	echo "Configuring MU dashboard file $WATCHFILE_MU."
	cp "$WATCHFILE_MU" "$WATCHFILE_MUBK"
	perl -pi -e 's/\r$//' "$WATCHFILE_MU"
	sed -i '' -e '$a\' "$WATCHFILE_MU"
	echo "Done. The original, unprocessed file has been archived as $WATCHFILE_MUBK."
fi


indir="$1"
WD=$(pwd)

status=0

echo "Switching to run directory $indir."
cd "$indir"

echo "Creating directory $rundir/restore/ to hold the original data files."
if [ -d "restore/" ]
then
	echo "Directory already exists."
else
	mkdir "restore/"
fi

echo "Preparing run files in $indir."

csv_files=*.csv
for f in $csv_files
do
	echo "Backing up csv file $f to $indir/restore/."
	cp "$f" "restore/$f"
	echo "Repairing known issues in $f."
	perl -pi -e 's/\r$//' "$f"
	sed $'1s/\xef\xbb\xbf//' < "$f" > "$f.tmp"
	mv "$f.tmp" "$f"
	echo "Done."
	
	#echo "Checking data type of file $f"
	is_ddpcr=$(head -n 1 "$f" | grep -c "Well,")
	is_asset=$(grep -c "Asset Tag ID," "$f")

	if [[ "$is_ddpcr" == "1" ]]
	then
		echo "File $f appears to be a ddPCR results file."
		echo "  Fixing known incompatibilities in $f."
		perl -pi -e 's/,Taget,/,Target,/i' "$f"
		perl -pi -e 's/,N1,/,SARS-CoV-2 N1,/i' "$f"
		perl -pi -e 's/,N2,/,SARS-CoV-2 N2,/i' "$f"
		perl -pi -e 's/,RP,/,Human RNase P,/i' "$f"
		perl -pi -e 's/µ/u/i' "$f"
		echo "Done."
	elif [[ "$is_asset" == "1" ]]
	then
		echo "File $f appears to be an asset file."
		echo "It has been backed up to restore/ but will be ignored."
		echo "You may want to process this file separately, using perl/2_updateSamples.pl."
	else
		echo "File $f is a csv file of unknown data format."
		echo "It has been backed up to restore/ but will be ignored."
	fi
done

xlsx_files=*.xlsx
for f in $xlsx_files
do
	echo "Backing up xlsx file $f to $indir/restore/."
	cp "$f" "restore/$f"
#	echo "File $f has been backed up but not modified in any other way."
#	echo "(Excel files from the testing laboratory are not processed further by this script.)"
done


cd "$WD"

exit $status

