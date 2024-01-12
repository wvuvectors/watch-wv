#! /bin/bash

target=$1

if [[ "$target" == "NWSS" ]]
then
	./feed_NWSS.sh
elif [[ "$target" == "WVWD" ]]
then
	./feed_WVWD.sh
else
	echo "Unknown feed target selected!"
	exit 1
fi


echo "All done! patchr_feed will now exit."
echo ""

exit 0
