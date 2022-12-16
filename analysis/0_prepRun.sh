#! /bin/bash

echo "******"

indir="$1"
WD=$(pwd)

runf="run_data.csv"
assf="assets.csv"

echo "Switching to run directory $indir."
cd "$indir"

echo "Creating directory restore/ to hold original data files."
if [ -d "restore/" ]
	echo "Directory already exists."
else
	mkdir "restore/"
fi

echo "Looking for the ddPCR results file, in csv format."

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

	echo "Checking csv data type..."
	is_ddpcr=$(head -n 1 "$f" | grep -c "Well,")
	is_asset=$(grep -c "Asset Tag ID," "$f")

	if [[ "$is_ddpcr" == "1" ]]
	then
		echo "Found ddPCR data file $f."
		echo "Renaming $f to '$runf' for easier processing..."
		mv "$f" "$runf"
		echo "Fixing known incompatibilities..."
#		perl -pi -e 's/,Taget,/,Target,/i' "$runf"
#		perl -pi -e 's/,N1,/,SARS-CoV-2 N1,/i' "$runf"
#		perl -pi -e 's/,N2,/,SARS-CoV-2 N2,/i' "$runf"
#		perl -pi -e 's/,RP,/,Human RNase P,/i' "$runf"
		perl -pi -e 's/µ/u/i' "$runf"
		echo "Done."
	elif [[ "$is_asset" == "1" ]]
	then
		echo "Found ddPCR data file $f."
		echo "Renaming $f to '$assf' for easier processing..."
		mv "$f" "$assf"
		echo "Extracting asset data..."
		grep "," "$assf" > "$assf"
		echo "Done."
	fi
done

cd "$WD"

echo "Finished run file processing."
echo "In the event of a catastrophic error, run the 'restore.sh' script to restore the original state. Then fix the error(s) and run this script again."
echo "******"

