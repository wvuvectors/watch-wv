#! /bin/bash

indir="$1"
WD=$(pwd)

now=$(date +'%Y-%m-%d.%k_%M')

WATCHFILE_MU="/Users/tpd0001/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/mu_dashboard.LATEST.tsv"
WATCHFILE_MUBK="/Users/tpd0001/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR DASHBOARD/ARCHIVED/mu_dashboard.$now.tsv"

if [ -f "$WATCHFILE_MU" ]
then
	cp "$WATCHFILE_MU" "$WATCHFILE_MUBK"
	perl -pi -e 's/\r$//' "$WATCHFILE_MU"
	sed -i '' -e '$a\' "$WATCHFILE_MU"
fi


cd "$indir"
csv_files=*.csv

for f in $csv_files
do
	cp "$f" "$f.bak"
	perl -pi -e 's/\r$//' "$f"
	sed $'1s/\xef\xbb\xbf//' < "$f" > "$f.tmp"
	mv "$f.tmp" "$f"
	is_ddpcr=$(head -n 1 "$f" | grep -c "Well,")
	is_asset=$(grep -c "Asset Tag ID," "$f")
	if [[ "$is_ddpcr" == "1" ]]
	then
		echo "Processing data file $f"
		perl -pi -e 's/,Taget,/,Target,/i' "$f"
		perl -pi -e 's/,N1,/,SARS-CoV-2 N1,/i' "$f"
		perl -pi -e 's/,N2,/,SARS-CoV-2 N2,/i' "$f"
		perl -pi -e 's/,RP,/,Human RNase P,/i' "$f"
		perl -pi -e 's/µ/u/i' "$f"
		mv "$f" "assay_plate.csv"
	elif [[ "$is_asset" == "1" ]]
	then
		echo "Processing asset file $f"
		grep "," "$f" > "assets_all.csv"
	fi
done

cd "$WD"
./1_compile_run.pl "$indir"
