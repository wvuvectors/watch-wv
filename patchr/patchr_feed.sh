#! /bin/bash

target=$1

if [[ "$target" == "NWSS" ]]
then
	./sh/feed_NWSS.sh
elif [[ "$target" == "WVDash" ]]
then
	./sh/feed_WVDash.sh
elif [[ "$target" == "WVBE" ]]
then
	./sh/feed_WVBE.sh
else
	echo "Unknown feed target selected! Please choose one of the following:"
	echo "NWSS (for Natl Wastewater Surveillance System)"
	echo "WVDash (for WV dashboard run by WVU)"
	echo "WVBE (for WV BPH breatheeasy Web site)"
	exit 1
fi


echo "All done! patchr_feed will now exit."
echo ""

exit 0
