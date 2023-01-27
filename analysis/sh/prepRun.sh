#! /bin/bash


indir="$1"
WD=$(pwd)

status=0

echo "Switching to run directory $indir."
cd "$indir"

echo "Creating directory restore/ to hold original data files."
if [ -d "restore/" ]
then
	echo "Directory already exists."
else
	mkdir "restore/"
fi

echo "Processing files in $indir."

csv_files=*.csv
for f in $csv_files
do
	echo "Backing up csv file $f to restore/"
	cp "$f" "restore/$f"
	echo "Repairing known issues in $f"
	perl -pi -e 's/\r$//' "$f"
	sed $'1s/\xef\xbb\xbf//' < "$f" > "$f.tmp"
	mv "$f.tmp" "$f"

	#echo "Checking data type of file $f"
	is_ddpcr=$(head -n 1 "$f" | grep -c "Well,")
	is_asset=$(grep -c "Asset Tag ID," "$f")

	if [[ "$is_ddpcr" == "1" ]]
	then
		echo "File $f appears to be a ddPCR results file."
		echo "  Fixing known incompatibilities in $f..."
#		perl -pi -e 's/,Taget,/,Target,/i' "$f"
#		perl -pi -e 's/,N1,/,SARS-CoV-2 N1,/i' "$f"
#		perl -pi -e 's/,N2,/,SARS-CoV-2 N2,/i' "$f"
#		perl -pi -e 's/,RP,/,Human RNase P,/i' "$f"
		perl -pi -e 's/µ/u/i' "$f"
	elif [[ "$is_asset" == "1" ]]
	then
		echo "File $f appears to be an asset file. It will be ignored. You may want to process this separately using 2_updateSamples.pl."
	else
		echo "File $f is a csv file of unknown data format. It will be ignored."
	fi
	echo "Done."
done

xlsx_files=*.xlsx
for f in $xlsx_files
do
	echo "Backing up xlsx file  $f to restore/"
	cp "$f" "restore/$f"
#	echo "File $f has been backed up but not modified in any other way."
#	echo "(Excel files from the testing laboratory are not processed further by this script.)"
done


cd "$WD"

exit $status

