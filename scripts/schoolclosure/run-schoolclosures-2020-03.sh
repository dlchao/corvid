#!/bin/bash
#script for generating config files for liberal leave, work from home, total school closure (non-reactive) sweep
#runs combinations of duration of closure, ascertainment of cases, delay of closure after first ascertained case
# After the runs are completed, the model outputs are moved to "results".
# This script will take a long time to run everything on a single core.
# It just runs one simulation at a time.
# Users might want to modify it to take advantage of multi-core systems 
# or clusters.

# run main school closure/liberal leave/work from home combinations
for r0 in 2.0 2.6; do
    printr0=${r0//./}
    for responseday in 60 74 100; do
	for printcompliance in 0 25 50; do
	    compliance="0."${printcompliance}
	    if test $compliance = "0.100"
	    then
		compliance=1
	    fi
	    # "0-day" school closures mean "no" school closure
	    cp template-config-seattle-liberalleave-workfromhome-schoolclosures config
	    sed -i "s/PARAMR0/$r0/g" config
	    sed -i "s/PRINTR0/$printr0/g" config
	    sed -i "s/DURATION/0/g" config
	    sed -i "s/RESPONSEDAY/$responseday/g" config
	    sed -i "s/SCHOOLPOLICY/none/g" config
	    sed -i "s/PRINTCOMPLIANCE/$printcompliance/g" config
	    sed -i "s/PARAMCOMPLIANCE/$compliance/g" config
	    ./corvid config
	    mv config tempconfig-seattle-$printr0-$responseday-school-none-0-llwfh-$printcompliance # just in case you want to see the config file later

	    for schoolpolicy in "all"; do
		for duration in 14 42 180; do
		    cp template-config-seattle-liberalleave-workfromhome-schoolclosures config
		    sed -i "s/PARAMR0/$r0/g" config
		    sed -i "s/PRINTR0/$printr0/g" config
		    sed -i "s/DURATION/$duration/g" config
		    sed -i "s/RESPONSEDAY/$responseday/g" config
		    sed -i "s/SCHOOLPOLICY/$schoolpolicy/g" config
		    sed -i "s/PRINTCOMPLIANCE/$printcompliance/g" config
		    sed -i "s/PARAMCOMPLIANCE/$compliance/g" config
		    ./corvid config
		    mv config tempconfig-seattle-$printr0-$responseday-school-$schoolpolicy-$duration-llwfh-$printcompliance # just in case you want to see the config file later
		done
	    done
	done
    done
done

# run a few more pre-determined scenarios
./corvid config-seattle20
./corvid config-seattle26
./corvid config-seattle-26-100-school-all-42-llwfh-0-heavylogging
./corvid config-seattle-26-1-school-all-365-llwfh-0-heavylogging
./corvid config-seattle-26-60-school-all-365-llwfh-0-heavylogging
./corvid config-seattle-26-60-school-all-42-llwfh-0-heavylogging
./corvid config-seattle-26-90-school-all-42-llwfh-0-heavylogging

# move all results to a new directory
mkdir results
mv Summary* Log* Tracts* Individuals* tempconfig* results

exit 0
