#! /bin/bash

indir=$1

cd "$indir"

# Convert Windows line endings to unix because of course
my_files=*.txt

for f in $my_files
do
	perl -pi -e 's/\r$//' "$f"
	tail -c1 < "$f" | read -r _ || echo >> "$f"
done

exit 0
