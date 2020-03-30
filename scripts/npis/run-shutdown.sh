#!/bin/bash
# script for running shutdown scenarios
# We test combinations of duration of shutdowns, and threshold number
# of cases before shutdowns.
# A "shutdown" is the simultaneous closure of schools, most workers stay
# home, and community contacts are reduced.
# This script will take a long time to run everything on a single core.
# It just runs one simulation at a time.
# Users might want to modify it to take advantage of multi-core systems 
# or clusters.

#  shutdowns
for r0 in 2.6; do
    for contactreduction in 0.75; do
	printr0=${r0//./}
	for threshold in 0.001 0.0001; do
#	    for duration in 28 56 90; do
	    for duration in 56 90; do
#		for compliance in 0.75 0.9; do
		for compliance in 0.9; do
		    for randomseed in 1 2 3 4 5 6 7 8 9 10; do
			cp template-config-seattle-shutdown configshutdown
			sed -i "s/PARAMR0/$r0/g" configshutdown
			sed -i "s/PRINTR0/$printr0/g" configshutdown
			sed -i "s/DURATION/$duration/g" configshutdown
			sed -i "s/RESPONSETHRESHOLD/$threshold/g" configshutdown
			sed -i "s/RANDOMSEED/$randomseed/g" configshutdown
			sed -i "s/COMPLIANCE/$compliance/g" configshutdown
			sed -i "s/CONTACTREDUCTION/$contactreduction/g" configshutdown
			./corvid configshutdown
		    done
		done
	    done
	done
    done
done

exit 0
