#!/bin/bash
# script for running "test and isolate" strategy
# We test combinations of isolation compliance and policy start time.
# After the runs are completed, the model outputs are moved to "results".
# This script will take a long time to run everything on a single core.
# It just runs one simulation at a time.
# Users might want to modify it to take advantage of multi-core systems 
# or clusters.

# test-and-isolate runs
for r0 in 2.6; do
    printr0=${r0//./}
    for responseday in 30 90; do
	for testdelay in 0; do
	    for testfraction in 0.75; do
                # test and isolate only
		cp template-config-seattle-testandisolate configt
		sed -i "s/PARAMR0/$r0/g" configt
		sed -i "s/PRINTR0/$printr0/g" configt
		sed -i "s/RESPONSEDAY/$responseday/g" configt
		sed -i "s/TESTFRACTION/$testfraction/g" configt
		sed -i "s/TESTDELAY/$testdelay/g" configt
		sed -i "s/QUARANTINECOMPLIANCE/0/g" configt # no quarantine
		sed -i "s/QUARANTINEDAYS/0/g" configt # no quarantine
		sed -i "s/AVPOLICY/none/g" configt # no antivirals
 		./corvid configt
		
                # test and isolate + AV
		cp template-config-seattle-testandisolate configav
		sed -i "s/PARAMR0/$r0/g" configav
		sed -i "s/PRINTR0/$printr0/g" configav
		sed -i "s/RESPONSEDAY/$responseday/g" configav
		sed -i "s/TESTFRACTION/$testfraction/g" configav
		sed -i "s/TESTDELAY/$testdelay/g" configav
		sed -i "s/QUARANTINECOMPLIANCE/0/g" configav # no quarantine
		sed -i "s/QUARANTINEDAYS/0/g" configav # no quarantine
		sed -i "s/AVPOLICY/treatmentonly/g" configav # antiviral for case
  		./corvid configav

                # test and isolate + home quarantine
		cp template-config-seattle-testandisolate configq
		sed -i "s/PARAMR0/$r0/g" configq
		sed -i "s/PRINTR0/$printr0/g" configq
		sed -i "s/RESPONSEDAY/$responseday/g" configq
		sed -i "s/TESTFRACTION/$testfraction/g" configq
		sed -i "s/TESTDELAY/$testdelay/g" configq
		sed -i "s/QUARANTINECOMPLIANCE/1.0/g" configq # great home quarantine
		sed -i "s/QUARANTINEDAYS/14/g" configq # home quarantine for 14 days
		sed -i "s/AVPOLICY/none/g" configq # antiviral for case
		./corvid configq
	    done
	done
    done
done

exit 0

# move results to a subdirectory
#mkdir results
#mv Summary* Tract* Log* Individual* results

exit 0
