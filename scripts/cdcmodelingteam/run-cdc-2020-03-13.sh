#!/bin/bash
# script for running CDC Modeling Team scenarios on March 13, 2020
# We test combinations of R0, duration of closure, and threshold number
# of cases before schools close.
# After the runs are completed, the model outputs are moved to "results".
# This script will take a long time to run everything on a single core.
# It just runs one simulation at a time.
# Users might want to modify it to take advantage of multi-core systems 
# or clusters.

./corvid config-seattle20-heavylogging-cdc-2020-03-13
./corvid config-seattle25-heavylogging-cdc-2020-03-13
./corvid config-seattle25-libleavewfh50-trigger0.001-cdc-2020-03-13
./corvid config-seattle25-school-libleavewfh50-trigger0.001-cdc-2020-03-13

# no intervention
for r0 in 2.5 2.0 3.0; do
    printr0=${r0//./}
    cp template-config-seattle-cdc-2020-03-13 config
    sed -i "s/PARAMR0/$r0/g" config
    sed -i "s/PRINTR0/$printr0/g" config
    sed -i "s/PRINTCOMPLIANCE/$printcompliance/g" config
    sed -i "s/PARAMCOMPLIANCE/$compliance/g" config
    ./corvid config
done

# school closures
for r0 in 2.5 2.0 3.0; do
    printr0=${r0//./}
    for trigger in 0.0001 0.001 0.01 0.1; do
	for duration in 14 28 56 140; do
	    cp template-config-seattleschool-cdc-2020-03-13 config
	    sed -i "s/PARAMR0/$r0/g" config
	    sed -i "s/PRINTR0/$printr0/g" config
	    sed -i "s/DURATION/$duration/g" config
	    sed -i "s/RESPONSETHRESHOLD/$trigger/g" config
	    ./corvid config
	done
    done
done

# move results to a subdirectory
#mkdir results
#mv Summary* Tract* Log* Individual* results

exit 0
