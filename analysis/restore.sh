#! /bin/bash

indir="$1"
WD=$(pwd)

echo "Restoring the run in '$indir' to its original state."

# Remove the processed files
cd "$indir"
files2rm=*.csv

for f in $files2rm
do
	rm "$f"
done


# Restore the original files
cd "restore/"
files2rs=*.csv

for f in $files2rs
do
	cp "$f" "../$f"
done


# Remove the restore folder
cd ../
rm -r "restore/"
cd "$WD"


echo "Done."

