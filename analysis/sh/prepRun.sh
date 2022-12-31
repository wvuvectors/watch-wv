#! /bin/bash

echo "******"
echo "Running prepRun.sh."
echo "******"


indir="$1"
WD=$(pwd)

status=0

echo "Switching to run directory $indir."
cd "$indir"

echo "Creating directory restore/ to hold original data files."
if [ -d "restore/" ]
	echo "Directory already exists."
else
	mkdir "restore/"
fi

echo "Processing csv format files for run compilation."

csv_files=*.csv
for f in $csv_files
do
	echo "Backing up $f to restore/ folder..."
	cp "$f" "restore/$f"
	echo "Repairing known file issues..."
	perl -pi -e 's/\r$//' "$f"
	sed $'1s/\xef\xbb\xbf//' < "$f" > "$f.tmp"
	mv "$f.tmp" "$f"
	echo "Done."

	echo "Checking data type of file $f..."
	is_ddpcr=$(head -n 1 "$f" | grep -c "Well,")
	is_asset=$(grep -c "Asset Tag ID," "$f")

	if [[ "$is_ddpcr" == "1" ]]
	then
		echo "File $f appears to be a ddPCR results file."
		echo "  Fixing known incompatibilities..."
#		perl -pi -e 's/,Taget,/,Target,/i' "$f"
#		perl -pi -e 's/,N1,/,SARS-CoV-2 N1,/i' "$f"
#		perl -pi -e 's/,N2,/,SARS-CoV-2 N2,/i' "$f"
#		perl -pi -e 's/,RP,/,Human RNase P,/i' "$f"
		perl -pi -e 's/µ/u/i' "$f"
		$status=1
		echo "Done."
	elif [[ "$is_asset" == "1" ]]
	then
		echo "File $f appears to be an asset file."
		echo "  This should be processed for new sample data using 2_updateSamples.pl."
		echo "Done."
	else
		echo "File $f is a csv file of unknown data format. It has been backed up but will be ignored otherwise."
		echo "Done."		
	fi
done

echo "Processing xlsx format files for run compilation."

xlsx_files=*.xlsx
for f in $xlsx_files
do
	echo "Backing up $f to restore/ folder..."
	cp "$f" "restore/$f"
	echo "Done."
	echo "File $f has been backed up but not modified in any other way."
	echo "(Excel files from the testing laboratory are not processed further by this script.)"
	echo "Done."
done


cd "$WD"

echo "******"
echo "Finished prepRun.sh."
echo "In the event of a catastrophic error, run:"
echo "  ./sh/rollBack.sh $indir"
echo "to restore the original state. Then fix the error(s) and run this script again."
echo "******"

exit $status

