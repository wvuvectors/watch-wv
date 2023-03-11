#! /bin/bash

indir="$1"
WD=$(pwd)

echo "Restoring the run in '$indir' to its original state."

cd "$indir"

# Remove the processed files
csv2rm=*.csv
for f in $csv2rm
do
	rm "$f"
done
xlsx2rm=*.xlsx
for f in $xlsx2rm
do
	rm "$f"
done


cd "restore/"

# Restore the original files
csv2rs=*.csv
for f in $csv2rs
do
	cp "$f" "../$f"
done
xlsx2rs=*.xlsx
for f in $xlsx2rs
do
	cp "$f" "../$f"
done


cd ../

# Remove the updates folder, if it exists
if [ -d "updates/" ]; then
	rm -r "updates/"
fi


# Remove the restore folder
#rm -r "restore/"

cd "$WD"
echo "Done."
echo "You may delete the restore folder after verifying the restore was successful."

exit 0
