#!/bin/bash
# script for "no intervention" runs
# "individual" files are output

# no intervention
for r0 in 2.0 2.6; do
    printr0=${r0//./}
    cp template-config-seattle configbaseline
    sed -i "s/PARAMR0/$r0/g" configbaseline
    sed -i "s/PRINTR0/$printr0/g" configbaseline
    ./corvid configbaseline
done

# move results to a subdirectory
#mkdir results
#mv Summary* Tract* Log* Individual* results

exit 0
